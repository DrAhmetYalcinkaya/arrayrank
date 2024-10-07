#' Identify Discordant Duplicates in Microarray Data (GPR Read by Limma)
#'
#' This function identifies discordant duplicate spots in microarray data based on user-defined fold-change and absolute difference thresholds. It operates on data typically read from GenePix Result (GPR) files and processed into an `RGList` object from the `limma` package. The function flags duplicates where the intensity measurements differ significantly, which can be useful for quality control in microarray analysis.
#'
#' @param data An `RGList` object containing microarray data, specifically with components `genes` and `R`.
#' @param fold A numeric value specifying the fold-change threshold for identifying discordance between duplicates. Defaults to `1.5`.
#' @param abs A numeric value specifying the absolute difference threshold for identifying discordance between duplicates. Defaults to `500`.
#' @param array_type A string defining the type of array being examined. Possible entries: "huprot", "chambered", "segmented".
#'
#' @return A data frame containing the discordant duplicates, including columns for array name, block, gene name, duplicate intensities, mean intensity, fold-change ratio, absolute difference, and a flag indicating discordance.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' # Assuming 'raw_data' is an RGList object obtained from read.gpr()
#' discordant_duplicates <- discordant(data = raw_data, array_type = "chambered")
#' }
discordant <- function(data, array_type, fold = 1.5, abs = 500) {
  df <- data$genes
  merged_df <- cbind(df, data$R)
  colnames(merged_df) <- basename(colnames(merged_df))

  long_df <- tidyr::pivot_longer(merged_df, cols = 6:ncol(merged_df))
  colnames(long_df)[6] <- "array_name"

  if (array_type == "chambered" || array_type == "segmented") {
    result <- long_df %>%
      dplyr::group_by(array_name, Block, Name) %>%
      dplyr::filter(dplyr::n() == 2) %>%
      dplyr::summarize(
        dup1 = dplyr::first(value),
        dup2 = dplyr::last(value),
        dup_mean = mean(value),
        ratio = ifelse(min(value) == 0, 999, round(abs(max(value) / min(value)), 2)),
        abs_diff = abs(dplyr::first(value) - dplyr::last(value)),
        .groups = 'drop'
      )
  } else if (array_type == "huprot") {
    result <- long_df %>%
      dplyr::group_by(array_name, Name) %>%
      dplyr::filter(dplyr::n() == 2) %>%
      dplyr::summarize(
        dup1 = dplyr::first(value),
        dup2 = dplyr::last(value),
        dup_mean = mean(value),
        ratio = ifelse(min(value) == 0, 999, round(abs(max(value) / min(value)), 2)),
        abs_diff = abs(dplyr::first(value) - dplyr::last(value)),
        .groups = 'drop'
      )
  } else {
    stop("Please define array_type. Accepts: huprot, chambered or segmented.")
  }

  flagged <- result %>%
    dplyr::mutate(flag = ifelse(
      (ratio > fold) &
        ((abs_diff > dup_mean * 0.25 & abs_diff > abs) &
           (dup_mean > 250)), TRUE, FALSE))

  disc <- subset(flagged, flag == TRUE)
  return(disc)
}
