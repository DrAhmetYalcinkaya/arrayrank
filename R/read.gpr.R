#' Read GenePix Result (GPR) Files and Perform Background Correction
#'
#' Reads all `.gpr` files within a selected directory. Utilizes the `limma` package to load raw data from GenePix Result files (.gpr). Defaults to performing background correction if specified.
#'
#' @param bgcorrect Logical; whether background correction should be performed. Default is `TRUE`.
#' @param fdata Character; whether to select `median` or `mean` foreground intensity values. Defaults to `"mean"`.
#'
#' @return An `RGList` object containing raw or background-corrected data from the `.gpr` files. If no files are found or if the user cancels the directory selection, `NULL` is returned.
#' @export
#'
#' @import limma
#' @import tcltk
#'
#' @examples
#' \dontrun{
#' # Run the function to select a directory and read .gpr files
#' raw_data <- read.gpr()
#' }
read.gpr <- function(fdata = "mean", bgcorrect = TRUE) {
  if (.Platform$OS.type == "windows") {
    directory <- choose.dir()
  } else {
    directory <- tk_choose.dir()
  }

  if (is.na(directory) || directory == "") {
    return(NULL)
  }

  files <- list.files(path = directory, pattern = "*.gpr", full.names = TRUE, include.dirs = FALSE)

  if (length(files) == 0) {
    message("No .gpr files found in the selected directory.")
    return(NULL)
  }

  source_type <- if (fdata == "mean") "genepix" else "genepix.median"
  raw_data <- limma::read.maimages(files, source = source_type)
  raw_data$targets$FileName <- basename(raw_data$targets$FileName)

  elements <- c("R", "G", "Rb", "Gb")
  if (all(elements %in% names(raw_data)) && bgcorrect) {
    message("Red, Green channels and background data found. Performing background correction.")
    post_bg <- limma::backgroundCorrect(raw_data, method = "subtract")
  } else {
    message("No background correction performed.")
    post_bg <- raw_data
  }

  return(post_bg)
}
