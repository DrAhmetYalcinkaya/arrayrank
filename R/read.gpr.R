#' Read GenePix Result (GPR) Files and Simple Background Correction if Applicable
#'
#' This function allows users to select a directory and read all `.gpr` files within it. The function uses the `limma` package to load the raw data from GenePix Result files, which are commonly used in microarray analysis.
#' Additionally, a check is performed to confirm the presence of background measurements. If present, limma::backgroundCorrect is applied with the `subtract` method.
#'
#' @param bgcorrect a boolean input specifying whether background correction should be done. Default is `TRUE`.
#'
#' @return A `RGList` object containing raw intensity or background-corrected data from the `.gpr` files. If no files are found or if the user cancels the directory selection, `NULL` is returned.
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
read.gpr <- function(bgcorrect = T) {
  if (.Platform$OS.type == "windows") {
    directory <- choose.dir()
  } else {
    directory <- tk_choose.dir()
  }

  if (is.na(directory) || directory == "") {
    return(NULL)
  }

  files <- list.files(path = directory, pattern = "*.gpr", full.names = T, include.dirs = F)

  if (length(files) == 0) {
    message("No .gpr files found in the selected directory.")
    return(NULL)
  }
  raw_data <- limma::read.maimages(files, source = "genepix")
  raw_data$targets$FileName <- basename(raw_data$targets$FileName)


  elements <- c("R", "G", "Rb", "Gb")
    if (all(elements %in% names(raw_data)) & bgcorrect == T) {
    message("Red, Green channels and background data found to exist, user requested bg correction; subtracting background.")
      post_bg <- limma::backgroundCorrect(raw_data, method = "subtract")
    } else {
    message("No background subtration performed.")
      post_bg <- raw_data
  }

  return(post_bg)
}
