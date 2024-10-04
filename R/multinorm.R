#' Normalize Data
#'
#' This function normalizes a data frame of protein expression data using different methods: quantile normalization, between-array normalization, or normalization based on a specific protein.
#'
#' @param data A data frame containing the protein expression data, where rows represent proteins and columns represent different samples.
#' @param method A character string specifying the normalization method to use. Options are `"quantile"`, `"arrays"`, or `"ig"`. Default is `"quantile"`.
#'
#' @return A data frame with normalized protein expression values.
#' @import limma
#' @import data.table
#' @import tidyr
#' @import writexl
#' @import readxl
#' @import dplyr
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
#' norm_data <- norm(data = your_data, method = "ig")
#' }

multinorm <- function(data, method = "quantile"){
  data <- as.data.frame(data)
  data[] <- lapply(data, as.numeric)

  if (ncol(data) == 0 || nrow(data) == 0) {
    stop("Data contains no valid numeric columns or rows.")
  }

  selected <- tolower(rownames(data)) %>%
    grep("goatanti-huig", .) %>%
    data[.,]
  if(nrow(selected) == 0){
    selected <- tolower(rownames(data)) %>%
      {data[grepl("igg", .) & grepl("100ng/ul", .), ]}
  }

  ortalamalar <- apply(selected, 2, mean)
  beta <- mean(unlist(selected))

  logistic_transform <- function(x, b = beta) {
    a <- ifelse(x < b/5, 0.005, ifelse(x < b, 0.0005, 0.00025))
    1 / (1 + exp(-a * (x - b)))
  }
  log_trans <- logistic_transform(ortalamalar)
  binary <- ifelse(log_trans > 0.075, T, F)

  if (sum(!binary) > 0) {
    post_data <- data[, binary, drop = F]
    backup <- data[, !binary, drop = F]
    message("Extremely low Ig responses detected in: ", paste(colnames(data)[!binary], collapse = ", "),
            "\nHighly suggested to exclude these indices from the data. Normalization will continue without altering these columns.")
  } else {
    post_data <- data
    backup <- NULL
    message("No columns identified as negative controls. Proceeding with normalization.")
  }

  if(method == "quantile"){
    normed <- as.data.frame(round(limma::normalizeQuantiles(post_data),1))
    data[, binary] <- normed
    output <- round(data,1)
    message("Normalized each column (observations) with limma-Quantile.")

  } else if(method == "arrays"){
    normed <- as.data.frame(round(limma::normalizeBetweenArrays(post_data),1))
    data[, binary] <- normed
    output <- round(data,1)
    message("Normalized between arrays (Blocks or Huprot) with limma-BetweenArrays.")

  } else if(method == "ig"){
    oranlar <- max(ortalamalar) / ortalamalar
    k <- round(log_trans * oranlar,3)
    k <- ifelse(k < 0.05, 1, k)
    ####placeholder for now:
    k <- ifelse(k < 0.7, 0.7, k)
    ####until better solution comes up.
    df2 <- df <- data
    for(i in 1:ncol(df)){
      df2[, i] <- df[, i] * k[i]
      output <- round(df2,1)
    }
    message("Normalized for Ig values.")
  } else {
    message("Please define normalization variables correctly.")
  }

  return(output)
}
