#' Perform Background Subtraction on Microarray Data
#'
#' This function applies background correction to microarray data using the "subtract" method from the `limma` package. Background correction is a crucial step in preprocessing microarray data to remove systematic noise.
#'
#' @param data An `RGList` object containing the raw microarray data, typically produced by the `limma` package's `read.maimages` function.
#'
#' @return An `RGList` object with background-corrected data.
#' @import limma
#' @import data.table
#' @import tidyr
#' @import writexl
#' @import readxl
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' corrected_data <- subtract.bg(data = your_RGList_object)
#' }
subtract.bg <- function(data){
  output <- limma::backgroundCorrect(data, method = "subtract")
  return(output)
}
