#' Detect Significant Protein Hits
#'
#' This function identifies and ranks significant protein hits by comparing the expression levels between control and non-control groups in microarray data. It applies several filters based on the standard deviation to mean ratio, minimum mean expression, fold change, and absolute thresholds.
#'
#' @param data A data frame containing the protein expression data, where rows represent proteins and columns represent different samples. The first column should contain protein names.
#' @param group_vector An optional character vector specifying the group labels for the samples. If not provided, the function will assume the first row of the data frame contains group labels.
#' @param controls A character string specifying the label of the control group.
#' @param examine An optional character vector specifying the non-control group labels to be analyzed. If NULL, all groups except the control group are analyzed.
#' @param sdmean A numeric value specifying the threshold for the standard deviation to mean ratio used for filtering. Default is 0.7.
#' @param min_mean A numeric value specifying the minimum mean expression level required for a protein to be considered. Default is 250.
#' @param fold_threshold A numeric value specifying the fold change threshold for detecting significant hits. Default is 10.
#' @param absolute_threshold A numeric value specifying the absolute expression threshold for detecting significant hits. Default is 7000.
#' @param fold_shared An integer specifying the minimum number of non-control samples that must exceed the fold change threshold. Default is 1.
#' @param abs_shared An integer specifying the minimum number of non-control samples that must exceed the absolute threshold. Default is 1.
#'
#' @return A data frame of ranked significant protein hits, including the highest fold change, the number of potential hits, control mean expression, and a ranking score.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Example usage of detect.hits function
#' ranked_proteins <- detect.hits(data = my_data,
#'                                controls = "Control",
#'                                examine = c("Treatment1", "Treatment2"))
#' }
detect.hits <- function(data, group_vector = NULL, controls, examine = NULL,
                        sdmean = 0.7, min_mean = 250,
                        fold_threshold = 10, absolute_threshold = 7000,
                        fold_shared = 1, abs_shared = 1) {
  # Convert to a data frame if it's not already (handling tibbles)
  if ("tbl_df" %in% class(data)) {
    data <- as.data.frame(data)
  }

  # Attempt to convert the first row (excluding the first column) to numeric
  first_row_numeric <- suppressWarnings(as.numeric(as.character(data[1, -1])))
  is_first_row_numeric <- !any(is.na(first_row_numeric))

  # If group_vector is provided, use it; otherwise, decide based on the first row content
  if (!is.null(group_vector)) {
    if (length(group_vector) != ncol(data) - 1) {
      stop("The length of 'group_vector' must match the number of samples in the data (columns except the protein name column).")
    }
    cat("A group_vector has been provided.\nAnalysis will assume the first row contains array data.\n\n")
    groups <- as.character(group_vector)
    proteins <- as.character(data[, 1])
    data <- data[, -1]
  } else if (is_first_row_numeric) {
    cat("The first row is numeric (array data?) and no group_vector has been provided.\nINFO: This function can accept groups as a character vector (group_vector = '')\nor as character entries in the first row of the dataset.\nPlease provide group info as such.\n")
    stop()
  } else {
    message("The first row appears to contain group data, analyzing based on these groups")
    groups <- as.character(data[1, -1])
    proteins <- as.character(data[-c(1), 1])
    data <- data[-c(1), -1]
  }

  # Convert column data to numeric and define proteins as row names
  data[] <- lapply(data, as.numeric)
  rownames(data) <- proteins
  print(table(groups))

  # Identify the control group indices
  control_indices <- which(groups == controls)
  if (length(control_indices) == 0) {
    stop("The specified control group was not found in the group annotations.")
  }

  # If no specific non-control groups are provided, use all except the control group
  if (is.null(examine)) {
    non_control_indices <- which(groups != controls)
  } else {
    non_control_indices <- which(groups %in% examine)
  }
  if (length(non_control_indices) == 0) {
    stop("No valid non-control groups were found based on your criteria.")
  }

  # Check and print out groups being analyzed
  message(paste("Selected control group:", controls))
  message(paste("Selected analysis group:", paste(unique(groups[non_control_indices]), collapse=", ")))

  # Check if the standard deviation-to-mean ratio is greater than sdmean value
  variation <- function(row, ratio) {
    sd_to_mean_ratio <- sd(row) / mean(row)
    return(sd_to_mean_ratio > ratio)
  }

  # Filter rows based on variance in the data
  filter_variance <- apply(as.matrix(data[, non_control_indices]), 1, variation, sdmean)
  edf2 <- data[filter_variance, ]
  proteins <- proteins[filter_variance]  # adjust the 'proteins' list with the same filter
  message(paste("Filtered data based on variance:", nrow(edf2)))

  # Filter rows based on the mean value of each row in the non-control data
  filter_minmean <- apply(edf2[, non_control_indices], 1, function(row) mean(row) > min_mean)
  edf3 <- edf2[filter_minmean, ]
  proteins <- proteins[filter_minmean]
  message(paste("Filtered data based on mean value:", nrow(edf3)))

  # Calculate the mean value for controls while excluding the highest value
  edf3$ControlMean <- apply(edf3[, control_indices], 1, function(row) mean(row[row != max(row)]))

  # Exclude proteins where controls have higher values compared to the max value in non-control samples
  filter_greatercontrols <- edf3$ControlMean < apply(edf3[, non_control_indices], 1, max)
  edf4 <- edf3[filter_greatercontrols, ]
  proteins <- proteins[filter_greatercontrols]
  message(paste("Filtered data based on greater-than-controls check:", nrow(edf4)))

  # Filter based on fold change
  filter_folddif <- rowSums(edf4[, non_control_indices] > fold_threshold * edf4$ControlMean) > (fold_shared - 1)
  edf5 <- edf4[filter_folddif, ]
  proteins <- proteins[filter_folddif]
  message(paste("Filtered data based on fold difference:", nrow(edf5)))

  # Filter based on absolute difference
  filter_absdif <- rowSums(edf5[, non_control_indices] > absolute_threshold) > (abs_shared - 1)
  edf6 <- edf5[filter_absdif, ]
  proteins <- proteins[filter_absdif]
  message(paste("Filtered data based on absolute difference:", nrow(edf6)))

  # Extract absolute values of non-control samples
  abs_values_non_controls <- edf6[, non_control_indices]

  # Calculate fold change for each protein
  fold_changes <- abs_values_non_controls / edf6$ControlMean
  # Number of hits
  number_of_hits <- rowSums((fold_changes > fold_threshold) & (abs_values_non_controls > absolute_threshold))
  # Base ranking score for each non-control
  fold_changes_with_abslog <- (fold_changes * log10(abs_values_non_controls))

  # Order for raw foldchange and basescore data -- separately
  fold_change_values <- t(apply(fold_changes, 1, function(row) row[order(-row)]))
  score_foldchange_and_log10 <- t(apply(fold_changes_with_abslog, 1, function(row) row[order(-row)]))

  # Sum the scores for each non-control to create final ranking score
  rankingscore <- (rowSums(score_foldchange_and_log10) * number_of_hits * log10(edf6$ControlMean))
  message("Rank scores calculated")

  # Add ranking score, highest fold change, number of hits, and BDmean to the dataframe
  abs_values_non_controls$highest_fold_change <- fold_change_values[, 1]
  abs_values_non_controls$potential_hits <- number_of_hits
  abs_values_non_controls$control_mean <- edf6$ControlMean
  abs_values_non_controls$rankingscore <- rankingscore
  message("Supportive info added.")

  # Add a rank column to the dataframe (convenience)
  abs_values_non_controls$rank <- rank(-abs_values_non_controls$rankingscore)

  # Match the length of 'Proteins' to the filtered data
  if (length(proteins) != nrow(abs_values_non_controls)) {
    stop("Mismatch in the number of rows after final filtering. Perhaps check the data for format or structure?")
  }

  # Combine results with protein and group information
  results <- data.frame(proteins, abs_values_non_controls, check.names = FALSE)
  ranked_prots <- results[order(results$rank), ]
  return(ranked_prots)
}
