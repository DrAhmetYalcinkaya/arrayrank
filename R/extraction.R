#' Extract Data from RGList Object
#'
#' This function extracts protein data, array information, and intensity values from an `RGList` object (typically produced by the `limma` package) and returns a combined data frame.
#'
#' @param data An `RGList` object containing gene expression data, typically produced by the `limma` package's `read.maimages` function. The object should include `genes`, `targets`, and `R` components.
#'
#' @return A data frame containing protein information, array identifiers, and corresponding intensity values, with one row per protein-array combination.
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
#' extracted_data <- extraction(data = your_RGList_object)
#' }
extraction <- function(data){
  proteins <- data$genes
  files <- data$targets
  values <- data$R
  array_count <- ncol(values)

  combined_list <- vector("list", array_count)

  for (i in 1:array_count) {
    df <- as.data.frame(proteins)
    df$array <- files[[1]][i]
    df$value <- values[, i]
    combined_list[[i]] <- df
  }
  df_out <- do.call(rbind, combined_list)
  return(df_out)
}
