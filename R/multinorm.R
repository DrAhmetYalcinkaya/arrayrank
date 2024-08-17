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
  if(method == "quantile"){
    output <- as.data.frame(limma::normalizeQuantiles(data))
    print("Normalized each column (observations) with limma-Quantile.")
  } else if(method == "arrays"){
    output <- as.data.frame(limma::normalizeBetweenArrays(data))
    print("Normalized between arrays (Blocks or Huprot) with limma-Betweenarrays.")
  } else if(method == "protein"){
    if(is.null(protein) || protein > nrow(data) || protein <= 0){
      stop("Please provide a valid row index for the protein parameter.")
    }
    min_max_ratio <- (1 / (data[protein,] / max(data[protein,])))
    output <- sweep(data, 2, as.numeric(min_max_ratio), `*`)
    print(paste("Normalized according to", rownames(data[protein,]), "values."))
  } else {
    stop("Please define normalization parameters correctly.")
  }
  return(output)
}
