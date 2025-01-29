#' Normalize Data
#'
#' This function normalizes a data frame of protein expression data using different methods: quantile normalization, between-array normalization, or normalization based on a specific protein.
#'
#' @param data A data frame containing the protein expression data, where rows represent proteins and columns represent different samples.
#' @param method A character string specifying the normalization method to use. Options are "quantile", "arrays", or "ig". Default is "quantile".
#'
#' @return A data frame with normalized protein expression values.
#' @export
#'
#' @import limma
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Quantile normalization
#' norm_data <- multinorm(data = your_data, method = "quantile")
#'
#' # Between-array normalization
#' norm_data <- multinorm(data = your_data, method = "arrays")
#'
#' # Protein-based normalization
#' norm_data <- multinorm(data = your_data, method = "ig")
#' }
multinorm <- function(data, method = "quantile") {
  data <- as.data.frame(data)
  data[] <- lapply(data, as.numeric)

  if (ncol(data) == 0 || nrow(data) == 0) {
    stop("Data contains no valid numeric columns or rows.")
  }

  selected <- tolower(rownames(data)) %>%
    grep("goatanti-huig", .) %>%
    data[.,]
  if (nrow(selected) == 0) {
    select <- tolower(rownames(data)) %>%
      {data[grepl("anti-human igg", .) & grepl("ng/ul", .), ]}
    selected["anti-human igg 100ng/ul",] <- select["Anti-human IgG 100ng/ul",]
    selected["anti-human igg 25ng/ul-undiluted",] <- select["Anti-human IgG 25ng/ul",] * 4
  }

  ortalamalar <- apply(selected, 2, mean)
  beta <- median(unlist(selected))

  logistic_trans <- Vectorize(function(x, b = beta, alt = 0.8, ust = 1.5) {
    if (x < b/4) {
      a = 0.009
      return(1 / (1 + exp(-a * (x - b))))
    } else {
      a = 0.001
      return(alt + (ust - alt) / (1 + exp(-a * (b - x))))
    }
  })
  log_trans <- logistic_trans(ortalamalar)
  binary <- ifelse(log_trans > 0.01, TRUE, FALSE)

  if (sum(!binary) > 0) {
    post_data <- data[, binary, drop = FALSE]
    backup <- data[, !binary, drop = FALSE]
    message("Extremely low Ig responses detected in: ", paste(colnames(data)[!binary], collapse = ", "),
            "\nHighly suggested to exclude these indices from the data. Normalization will continue without altering these columns.")
  } else {
    post_data <- data
    backup <- NULL
    message("No columns identified as negative controls. Proceeding with normalization.")
  }

  if (method == "quantile") {
    normed <- as.data.frame(round(limma::normalizeQuantiles(post_data), 1))
    data[, binary] <- normed
    output <- round(data, 1)
    message("Normalized each column (observations) with limma-Quantile.")

  } else if (method == "arrays") {
    normed <- as.data.frame(round(limma::normalizeBetweenArrays(post_data), 1))
    data[, binary] <- normed
    output <- round(data, 1)
    message("Normalized between arrays (Blocks or Huprot) with limma-BetweenArrays.")

  } else if (method == "ig") {
    df <- data
    for (i in 1:ncol(df)) {
      k <- ifelse(log_trans < 0.1, 1, log_trans)
      df[, i] <- data[, i] * k[i]
      output <- round(df, 1)
    }
    message("Normalized for Ig values.")
  } else {
    message("Please define normalization variables correctly.")
  }

  return(output)
}
