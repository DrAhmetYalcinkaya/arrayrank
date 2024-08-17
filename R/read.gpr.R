#' Read GenePix Result (GPR) Files
#'
#' This function allows users to select a directory and read all `.gpr` files within it. The function uses the `limma` package to load the raw data from GenePix Result files, which are commonly used in microarray analysis.
#'
#' @return A `RGList` object containing raw intensity data from the `.gpr` files. If no files are found or if the user cancels the directory selection, `NULL` is returned.
#' @export
#'
#' @import limma
#' @import tcltk
#' @import data.table
#' @import tidyr
#' @import writexl
#' @import readxl
#'
#' @examples
#' \dontrun{
#' # Run the function to select a directory and read .gpr files
#' raw_data <- read.gpr()
#' }
read.gpr <- function() {
  if (.Platform$OS.type == "windows") {
    directory <- choose.dir()
  } else {
    directory <- tk_choose.dir()
  }
  
  if (is.na(directory) || directory == "") {
    return(NULL)
  }
  
  files <- list.files(path = directory, pattern = "*.gpr", full.names = TRUE)
  
  if (length(files) == 0) {
    message("No .gpr files found in the selected directory.")
    return(NULL)
  }
  
  raw_data <- limma::read.maimages(files, source = "genepix")
  
  return(raw_data)
}