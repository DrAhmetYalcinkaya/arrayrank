#' Detect Significant Hits using an Isolation Forest Anomaly Pipeline
#'
#' This function takes a tidy data frame and identifies significant hits using a
#' two-pass anomaly detection pipeline. It is the core analysis engine designed
#' to work with pre-formatted data where samples are rows and analytes are columns.
#'
#' @param data A tidy data frame with samples as rows. It must contain a sample
#'   identifier column, a cohort column, and numeric columns for each analyte.
#' @param sample_id_col A character string specifying the name of the column in
#'   `data` that contains the unique sample identifiers.
#' @param cohort_col A character string specifying the name of the column in
#'   `data` that contains the group/cohort labels.
#' @param controls A character string specifying the label of the control group.
#' defaults to NULL if working without a control group.
#' @param examine An optional character vector specifying the non-control
#'   group labels to be analyzed. If NULL, all groups except the control group
#'   are analyzed.
#' @param anomaly_scale A numeric value for the composite anomaly score threshold, used
#'   to define the initial set of anomalies for threshold setting. Default is 0.6.
#' @param fold_threshold A numeric value specifying the fold change threshold for
#'   detecting significant hits in the final classification. Default is 10.
#' @param penalize A numeric factor used to adjust the score of hits found in the
#'   control group. A value of 0 (default) nullifies their contribution, while a
#'   negative value (e.g., -0.5) actively penalizes the analyte's rank score.
#' @param detailed_results A logical value. If FALSE (default), the function runs in
#'   a high-performance mode and returns only a single, ranked data frame of hits.
#'   If TRUE, it returns a full list object containing detailed, sample-level
#'   results for every analyte.
#' @param mad_multiplier A numeric multiplier for the Median Absolute Deviation (MAD)
#'   used in setting the sigmoid inflection point. Default is 2.
#' @param sigmoid_k A numeric value for the slope of the sigmoid curve. Default is 0.0005.
#' @param ntrees A numeric value for the number of trees to build in the
#'   isolation forest. Default is 200.
#'
#' @return By default (`detailed_results = FALSE`), a single data frame of ranked
#'   significant hits. If `detailed_results = TRUE`, a list containing two summary
#'   data frames (`summary_hits_only`, `summary_all_analytes`) and detailed results
#'   for each analyte. Both summary types will include columns for each sample's
#'   raw value, with examine groups appearing before control groups.
#' @export
#'
#' @importFrom dplyr %>% filter mutate select left_join arrange
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom isotree isolation.forest
#'
#' @examples
#' \dontrun{
#' # First, prepare data using the updated assemble.output function
#' all_data <- assemble.output(key = sample_key, norm_df = normalized_data)
#'
#' # Run isolate.hits on the pre-tidied data frame
#' ranked_hits <- isolate.hits(data = all_data$tidied_analysis_set,
#'                             sample_id_col = "Array",
#'                             cohort_col = "Group",
#'                             controls = "Control",
#'                             examine = "Treatment",
#'                             penalize = -0.5)
#' }
isolate.hits <- function(data,
                         sample_id_col,
                         cohort_col,
                         controls = NULL,
                         examine = NULL,
                         mad_multiplier = 2,
                         sigmoid_k = 0.0005,
                         anomaly_scale = 0.6,
                         fold_threshold = 10,
                         ntrees = 200,
                         penalize = 0,
                         detailed_results = FALSE) {
  # main check and warning:
  if(is.null(controls) & is.null(examine)) {
    message("WARNING: Having BOTH the examine and controls unspecified is ill-advised, but the analysis will still _try to_ run...")
  }

  # assessing if a targeted examination group is defined
  if (is.null(examine)) {
    examine <- unique(data[[cohort_col]])
    message("Target examination group unspecified. Will attempt to 'examine' all groups except controls.")
  }
  #checking presence of control group
  if (is.null(controls)) {
    baseline_group <- examine
    use_examine_as_baseline <- TRUE
    message("Control group not defined. Using the examine group instead (self-control).")
    penalize <- 0
  } else {
    baseline_group <- controls
    use_examine_as_baseline <- FALSE
  }

  all_groups <- unique(c(controls, examine))

  # check input
  analysis_data <- data
  if (!all(all_groups %in% unique(analysis_data[[cohort_col]]))) {
    stop("One or more specified cohort groups not found in the data.")
  }
  analysis_data <- analysis_data %>%
    dplyr::filter(!!sym(cohort_col) %in% all_groups)
  message("Running targeted analysis on cohort(s): ", paste(examine, collapse=", "))

  id_cols <- c(sample_id_col, cohort_col)
  numerical_analytes <- setdiff(names(analysis_data)[sapply(analysis_data, is.numeric)], id_cols)
  if (length(numerical_analytes) == 0) { stop("No numerical columns found to analyze.") }
  message("Analyzing ", length(numerical_analytes), " analytes across ", nrow(analysis_data), " samples...")

  # do detailed output
  if (detailed_results) {
    individual_results <- lapply(numerical_analytes, function(analyte_name) {
      result_shell <- analysis_data %>% dplyr::select(all_of(id_cols))
      x_raw <- analysis_data[[analyte_name]]
      n_complete <- sum(!is.na(x_raw))
      if (sd(x_raw, na.rm = TRUE) == 0 || n_complete < 5) { return(NULL) }

      iso <- isotree::isolation.forest(data.frame(value = x_raw), sample_size = n_complete, ntrees = ntrees)
      scores <- predict(iso, data.frame(value = x_raw), type = "score")

      baseline_mask <- analysis_data[[cohort_col]] %in% baseline_group & !is.na(x_raw)
      baseline_raws_med <- median(x_raw[baseline_mask], na.rm = TRUE)

      med <- median(x_raw, na.rm = TRUE)
      mad <- mad(x_raw, na.rm = TRUE)
      inflec <- med + (mad_multiplier * mad)
      raw_sigmoid <- 1 / (1 + exp(-sigmoid_k * (x_raw - inflec)))
      score_sigmoid_comp <- scores * raw_sigmoid
      FCvsBaseline <- x_raw / (baseline_raws_med + 1e-9)

      anomaly <- score_sigmoid_comp >= anomaly_scale & !is.na(score_sigmoid_comp)
      absolute_raw_threshold <- min(x_raw[anomaly], na.rm = TRUE)
      hit_flag <- (x_raw >= absolute_raw_threshold) & (FCvsBaseline >= fold_threshold)
      hit_flag[is.na(hit_flag)] <- FALSE

      hit_value <- rep(NA_real_, length(x_raw))
      is_hit <- hit_flag == TRUE
      if (any(is_hit)) {
        hit_value[is_hit] <- (log10(x_raw[is_hit])) * (log2(FCvsBaseline[is_hit])) * (score_sigmoid_comp[is_hit])
        hit_value[is.infinite(hit_value) | is.nan(hit_value)] <- 0

        if (!use_examine_as_baseline) {
          is_control_group <- analysis_data[[cohort_col]] %in% controls
          control_hits_mask <- is_hit & is_control_group
          if(any(control_hits_mask)) {
            hit_value[control_hits_mask] <- hit_value[control_hits_mask] * penalize
          }
        }
      }
      final_analyte_tbl <- result_shell %>%
        dplyr::mutate(
          raw_value = x_raw, if_score = scores, sigmoid = raw_sigmoid,
          score_sigmoided = score_sigmoid_comp, fold_change = FCvsBaseline,
          anomaly = anomaly, hit_flag = hit_flag, hit_value = hit_value
        )
      return(final_analyte_tbl)
    })

    names(individual_results) <- numerical_analytes
    individual_results <- individual_results[!sapply(individual_results, is.null)]

    summary_list <- lapply(names(individual_results), function(analyte_name) {
      result_df <- individual_results[[analyte_name]]
      target_hits <- result_df[result_df$hit_flag & !result_df[[cohort_col]] %in% controls, ]
      control_hits <- result_df[result_df$hit_flag & result_df[[cohort_col]] %in% controls, ]
      target_score <- sum(target_hits$hit_value, na.rm = TRUE)
      control_penalty <- sum(control_hits$hit_value, na.rm = TRUE)
      total_hit_val <- target_score + control_penalty
      num_hits_val <- nrow(target_hits)

      if (num_hits_val == 0) {
        data.frame(
          analyte = analyte_name, num_hits = 0, control_hit_penalty = control_penalty,
          total_hit_value = total_hit_val, rank_hit_value = 0,
          min_hit_raw_value = NA_real_, max_hit_raw_value = NA_real_,
          implied_cutoff = NA_real_,
          min_hit_score = NA_real_, max_hit_score = NA_real_,
          min_hit_fc = NA_real_, max_hit_fc = NA_real_
        )
      } else {
        all_hits <- result_df[result_df$hit_flag, ]
        non_hits <- result_df[!result_df$hit_flag, ]

        min_hit_val <- min(all_hits$raw_value, na.rm = TRUE)
        highest_non_hit_val <- max(non_hits$raw_value, na.rm = TRUE)

        implied_cutoff_val <- if(is.infinite(highest_non_hit_val)) {
          NA_real_
        } else {
          mean(c(min_hit_val, highest_non_hit_val))
        }

        rank_hit_val <- total_hit_val * num_hits_val
        data.frame(
          analyte = analyte_name, num_hits = num_hits_val,
          control_hit_penalty = control_penalty, total_hit_value = total_hit_val,
          rank_hit_value = rank_hit_val,
          min_hit_raw_value = min_hit_val,
          max_hit_raw_value = max(all_hits$raw_value, na.rm = TRUE),
          implied_cutoff = implied_cutoff_val,
          min_hit_score = min(target_hits$score_sigmoided, na.rm = TRUE),
          max_hit_score = max(target_hits$score_sigmoided, na.rm = TRUE),
          min_hit_fc = min(target_hits$fold_change, na.rm = TRUE),
          max_hit_fc = max(target_hits$fold_change, na.rm = TRUE)
        )
      }
    })
    summary_all_df <- do.call(rbind, summary_list)
    summary_hits_only_df <- summary_all_df %>% dplyr::filter(num_hits > 0)

    if(nrow(summary_hits_only_df) > 0) {
      summary_hits_only_df <- summary_hits_only_df[order(summary_hits_only_df$rank_hit_value, decreasing = TRUE), ]
      summary_hits_only_df$rank <- 1:nrow(summary_hits_only_df)
    }

  } else { # do summary (fast) output
    summary_list <- lapply(numerical_analytes, function(analyte_name) {
      x_raw <- analysis_data[[analyte_name]]
      n_complete <- sum(!is.na(x_raw))
      if (sd(x_raw, na.rm = TRUE) == 0 || n_complete < 5) { return(NULL) }

      iso <- isotree::isolation.forest(data.frame(value = x_raw), sample_size = n_complete, ntrees = ntrees)
      scores <- predict(iso, data.frame(value = x_raw), type = "score")

      baseline_mask <- analysis_data[[cohort_col]] %in% baseline_group & !is.na(x_raw)
      baseline_raws_med <- median(x_raw[baseline_mask], na.rm = TRUE)

      med <- median(x_raw, na.rm = TRUE)
      mad <- mad(x_raw, na.rm = TRUE)
      inflec <- med + (mad_multiplier * mad)
      raw_sigmoid <- 1 / (1 + exp(-sigmoid_k * (x_raw - inflec)))
      score_sigmoid_comp <- scores * raw_sigmoid
      FCvsBaseline <- x_raw / (baseline_raws_med + 1e-9)

      anomaly <- score_sigmoid_comp >= anomaly_scale & !is.na(score_sigmoid_comp)
      absolute_raw_threshold <- min(x_raw[anomaly], na.rm = TRUE)
      hit_flag <- (x_raw >= absolute_raw_threshold) & (FCvsBaseline >= fold_threshold)
      hit_flag[is.na(hit_flag)] <- FALSE

      temp_df <- data.frame(
        cohort = analysis_data[[cohort_col]],
        raw_value = x_raw,
        score_sigmoided = score_sigmoid_comp,
        fold_change = FCvsBaseline,
        hit_flag = hit_flag
      )

      target_hits <- temp_df[temp_df$hit_flag & !temp_df$cohort %in% controls, ]
      control_hits <- temp_df[temp_df$hit_flag & temp_df$cohort %in% controls, ]

      num_hits_val <- nrow(target_hits)

      target_hit_values <- (log10(target_hits$raw_value)) * (log2(target_hits$fold_change)) * (target_hits$score_sigmoided)
      target_hit_values[is.infinite(target_hit_values) | is.nan(target_hit_values)] <- 0
      target_score <- sum(target_hit_values, na.rm = TRUE)

      control_hit_values <- (log10(control_hits$raw_value)) * (log2(control_hits$fold_change)) * (control_hits$score_sigmoided)
      control_hit_values[is.infinite(control_hit_values) | is.nan(control_hit_values)] <- 0
      control_penalty <- sum(control_hit_values * penalize, na.rm = TRUE)

      total_hit_val <- target_score + control_penalty

      if (num_hits_val == 0) {
        data.frame(
          analyte = analyte_name, num_hits = 0, control_hit_penalty = control_penalty,
          total_hit_value = total_hit_val, rank_hit_value = 0,
          min_hit_raw_value = NA_real_, max_hit_raw_value = NA_real_,
          implied_cutoff = NA_real_,
          min_hit_score = NA_real_, max_hit_score = NA_real_,
          min_hit_fc = NA_real_, max_hit_fc = NA_real_
        )
      } else {
        all_hits <- temp_df[temp_df$hit_flag, ]
        non_hits <- temp_df[!temp_df$hit_flag, ]

        min_hit_val <- min(all_hits$raw_value, na.rm = TRUE)
        highest_non_hit_val <- max(non_hits$raw_value, na.rm = TRUE)

        implied_cutoff_val <- if(is.infinite(highest_non_hit_val)) {
          NA_real_
        } else {
          mean(c(min_hit_val, highest_non_hit_val))
        }

        rank_hit_val <- total_hit_val * num_hits_val
        data.frame(
          analyte = analyte_name, num_hits = num_hits_val,
          control_hit_penalty = control_penalty, total_hit_value = total_hit_val,
          rank_hit_value = rank_hit_val,
          min_hit_raw_value = min_hit_val,
          max_hit_raw_value = max(all_hits$raw_value, na.rm = TRUE),
          implied_cutoff = implied_cutoff_val,
          min_hit_score = min(target_hits$score_sigmoided, na.rm = TRUE),
          max_hit_score = max(target_hits$score_sigmoided, na.rm = TRUE),
          min_hit_fc = min(target_hits$fold_change, na.rm = TRUE),
          max_hit_fc = max(target_hits$fold_change, na.rm = TRUE)
        )
      }
    })

    summary_all_df <- do.call(rbind, summary_list)
    summary_hits_only_df <- summary_all_df %>% dplyr::filter(num_hits > 0)

    if(nrow(summary_hits_only_df) > 0) {
      summary_hits_only_df <- summary_hits_only_df[order(summary_hits_only_df$rank_hit_value, decreasing = TRUE), ]
      summary_hits_only_df$rank <- 1:nrow(summary_hits_only_df)
    }
  }

  # wrangle data for output
  cohort_order <- unique(c(examine, controls))

  long_raw_data <- analysis_data %>%
    dplyr::select(all_of(id_cols), all_of(numerical_analytes)) %>%
    dplyr::mutate(!!sym(cohort_col) := factor(!!sym(cohort_col), levels = cohort_order)) %>%
    dplyr::arrange(!!sym(cohort_col)) %>%
    tidyr::pivot_longer(
      cols = all_of(numerical_analytes),
      names_to = "analyte",
      values_to = "raw_value"
    )

  wide_raw_data <- long_raw_data %>%
    tidyr::pivot_wider(
      id_cols = "analyte",
      names_from = all_of(sample_id_col),
      values_from = "raw_value"
    )

  # summary data
  summary_all_df <- dplyr::left_join(summary_all_df, wide_raw_data, by = "analyte")
  summary_hits_only_df <- dplyr::left_join(summary_hits_only_df, wide_raw_data, by = "analyte")

  if (detailed_results) {
    final_output <- c(
      list(summary_hits_only = summary_hits_only_df),
      list(summary_all_analytes = summary_all_df),
      individual_results
    )
    return(final_output)
  } else {
    return(summary_hits_only_df)
  }
}
