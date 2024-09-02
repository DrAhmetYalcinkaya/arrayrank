#' Store Data Frame for Export
#'
#' This function takes a data frame, adds protein names as a separate column, and prepares it for export, typically into an Excel file.
#'
#' @param data A data frame with protein expression data where rows represent proteins and columns represent different samples or calculated values. The row names should contain the protein names.
#'
#' @return A data frame with an additional column for protein names, ready for export to Excel or other formats.
#' @import data.table
#' @import tidyr
#' @import writexl
#' @import readxl
#' @import tcltk
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of store.df function
#' data_to_store <- store.df(data = your_data)
#' }
store.df <- function(data){
  protein_names <- rownames(data)
  output <- cbind(protein_names, data)
  message("Storable dataframe created for export into Excel.")
  return(output)
}
