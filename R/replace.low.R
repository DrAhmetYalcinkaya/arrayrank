#' Replace Low Values in Data
#'
#' This function replaces values in a data frame that are below a specified offset with a corrected value. Users can choose different strategies for replacing low values.
#'
#' @param data A data frame containing numeric values, such as protein expression data.
#' @param offset A numeric value specifying the threshold below which values will be replaced. Default is 10.
#' @param strategy A character string specifying the replacement strategy: "computed" for a computed value (default), "fixed" for a fixed replacement value, or "random" for adding small random noise.
#' @param fix A numeric value used only when `strategy` is "fixed". Default is defined as the `offset` value given by the user.
#'
#' @return A data frame with low values replaced based on the specified offset and strategy.
#' @export
#'
#' @examples
#' \dontrun{
#' # Replace low values in the data frame with a computed value
#' corrected_data <- replace.low(data = your_data, offset = 10)
#'
#' # Replace low values with a fixed value
#' corrected_data <- replace.low(data = your_data, offset = 10, strategy = "fixed", fix = 5)
#'
#' # Replace low values with small random noise
#' corrected_data <- replace.low(data = your_data, offset = 10, strategy = "random")
#' }
replace.low <- function(data, offset = 10, strategy = "computed", fix = offset) {
  data <- as.data.frame(data)
  numeric_columns <- sapply(data, is.numeric)

  replaced_n <- sum(data[numeric_columns] < offset, na.rm = TRUE)

  if (strategy == "computed") {
    data[numeric_columns][data[numeric_columns] < offset] <- abs(data[numeric_columns][data[numeric_columns] < offset] / offset + offset)
  } else if (strategy == "fixed") {
    data[numeric_columns][data[numeric_columns] < offset] <- fix
  } else if (strategy == "random") {
    data[numeric_columns][data[numeric_columns] < offset] <- runif(replaced_n, min = offset * 0.5, max = offset * 1.5)
  } else {
    stop("Invalid strategy. Choose 'computed', 'fixed', or 'random'.")
  }
  message(paste("Replaced", replaced_n, "values that are lower than the selected offset value of", offset, "with", strategy, "strategy"))

  return(data)
}
