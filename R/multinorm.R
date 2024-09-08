#' Normalize Data
#'
#' This function normalizes a data frame of protein expression data using different methods: quantile normalization, between-array normalization, or normalization based on a specific protein.
#'
#' @param data A data frame containing the protein expression data, where rows represent proteins and columns represent different samples.
#' @param method A character string specifying the normalization method to use. Options are `"quantile"`, `"arrays"`, or `"protein"`. Default is `"quantile"`.
#' @param protein An optional integer specifying the row index of the protein to use for normalization when `method = "protein"`.
#'
#' @return A data frame with normalized protein expression values.
#' @import limma
#' @import data.table
#' @import tidyr
#' @import writexl
#' @import readxl
#' @export
#'
#' @examples
#' \dontrun{
#' # Quantile normalization
#' norm_data <- norm(data = your_data, method = "quantile")
#'
#' # Between-array normalization
#' norm_data <- norm(data = your_data, method = "arrays")
#'
#' # Protein-based normalization
#' norm_data <- norm(data = your_data, method = "protein", protein = 1)
#' }

multinorm <- function(data, method = "quantile", protein = NULL){
  data <- as.data.frame(data)
  data[] <- lapply(data, as.numeric)

  original_col_order <- colnames(data)

  if (ncol(data) == 0 || nrow(data) == 0) {
    stop("Data contains no valid numeric columns or rows.")
  }

  subject_medians <- apply(data, 2, median)
  overall_median <- median(unlist(data))
  lows <- apply(data, 2, function(x) x < overall_median)
  n_lows <- colSums(lows)
  cut_lows <- which(n_lows > (0.9 * nrow(data)))

  excluded_columns <- NULL
  if (length(cut_lows) > 0) {
    excluded_columns <- data[, cut_lows, drop = FALSE]

    message("It appears that ", paste(colnames(data)[cut_lows], collapse = ", "),
            " are negative controls (non-sample).\nHighly suggested to exclude these indices from normalization.")
    response <- readline(prompt = "Would you like to exclude these indices from normalization? Enter 'yes' or 'no': ")

    if (tolower(response) == "yes") {
      data <- data[, -cut_lows, drop = FALSE]
      message("Excluded the identified negative control columns from normalization.")
    } else {
      excluded_columns <- NULL
      message("Proceeding without excluding any columns.")
    }
  } else {
    message("No columns identified as negative controls. Proceeding with normalization.")
  }

  if(method == "quantile"){
    output <- as.data.frame(limma::normalizeQuantiles(data))
    message("Normalized each column (observations) with limma-Quantile.")

  } else if(method == "arrays"){
    output <- as.data.frame(limma::normalizeBetweenArrays(data))
    message("Normalized between arrays (Blocks or Huprot) with limma-BetweenArrays.")

  } else if(method == "protein"){
    if(is.null(protein) || protein > nrow(data) || protein <= 0){
      stop("Please provide a valid row index for the protein parameter. Should range from 1 to ", nrow(data), " (the latter is protein count)")
    }
    protein_mean <- mean(as.numeric(data[protein,]))
    mean_sample_ratio <- (protein_mean / data[protein,])
    print(mean_sample_ratio)
    output <- sweep(data, 2, as.numeric(mean_sample_ratio), `*`)
    message(paste("Normalized according to", rownames(data)[protein], "values."))
    message(paste("NOTES: Max value in input data was:", max(data), "\nMax value in output data is:", round(max(output),2),
          "\nMax values greater than 2-3 fold of 65000 may cause overestimation in hit detection.",
          "\nIn such a case, please consider quantile or between-arrays normalization instead of protein-based normalization."))

  } else {
    stop("Please define normalization parameters correctly.")
  }

  if (!is.null(excluded_columns)) {
    output <- cbind(output, excluded_columns)
    output <- output[, original_col_order]
  }

  return(output)
}
