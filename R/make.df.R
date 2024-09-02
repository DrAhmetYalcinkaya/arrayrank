#' Create a Data Frame for Different Array Types
#'
#' This function processes data from different types of protein arrays (chambered, segmented, huprot) and returns a wide-format data frame with processed and rounded numeric values.
#'
#' @param data A data frame or data.table containing the protein expression data. It should include columns for array, Block (for chambered/segmented arrays), Name, and value.
#' @param array_type A character string specifying the type of array to process. Options are `"chambered"`, `"segmented"`, or `"huprot"`.
#'
#' @return A wide-format data frame with processed and rounded numeric values. The column names are adjusted based on the array type.
#' @import data.table
#' @import tidyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage for chambered/segmented arrays
#' processed_data <- make.df(data = your_data, array_type = "chambered")
#'
#' # Example usage for huprot arrays
#' processed_data <- make.df(data = your_data, array_type = "huprot")
#' }
make.df <- function(data, array_type) {
  # Convert data to data.table format
  df <- data.table::setDT(data)

  # Validate array_type parameter
  if (!(array_type %in% c("chambered", "segmented", "huprot"))) {
    stop("Invalid array_type. Please specify 'chambered', 'segmented', or 'huprot'.")
  }

  # Define required columns based on array type
  required_columns <- if (array_type %in% c("chambered", "segmented")) {
    c("array", "Block", "Name", "value")
  } else {
    c("array", "Name", "value")
  }

  # Check if the data contains the required columns
  if (!all(required_columns %in% colnames(df))) {
    stop(paste("Data must contain the following columns:", paste(required_columns, collapse = ", ")))
  }

  if (array_type == "chambered" || array_type == "segmented") {
    long_result <- df[, .(Mvalue = mean(.SD$value, na.rm = TRUE)), by = .(array, Block, Name)]
    wide_result <- tidyr::pivot_wider(long_result, names_from = Name, values_from = Mvalue)
    message("Data calculated/arranged for chambered/segmented arrays")

    transposed_result <- as.data.frame(t(wide_result))
    file_names <- as.character(transposed_result[1, ])
    array <- paste0("Arr", sapply(basename(file_names), function(x) gsub("[^0-9]", "", x)))
    new_colname <- paste0(array, "-Block", transposed_result[2,])
    corrected_transposed_result <- transposed_result[-c(1, 2), ]
    colnames(corrected_transposed_result) <- new_colname
    numeric_df <- as.data.frame(lapply(corrected_transposed_result, as.numeric))
    numeric_df <- as.data.frame(apply(numeric_df, 2, function(x) round(x, 2)))
    rownames(numeric_df) <- rownames(corrected_transposed_result)

  } else if (array_type == "huprot") {
    long_result <- df[, .(Mvalue = mean(.SD$value, na.rm = TRUE)), by = .(array, Name)]
    wide_result <- tidyr::pivot_wider(long_result, names_from = Name, values_from = Mvalue)
    message("Data calculated/arranged for huprot arrays")

    transposed_result <- as.data.frame(t(wide_result))
    file_names <- as.character(transposed_result[1, ])
    array <- paste0("Arr", sapply(basename(file_names), function(x) gsub("[^0-9]", "", x)))
    new_colname <- array
    corrected_transposed_result <- transposed_result[-1, ]
    colnames(corrected_transposed_result) <- new_colname
    numeric_df <- as.data.frame(lapply(corrected_transposed_result, as.numeric))
    numeric_df <- as.data.frame(apply(numeric_df, 2, function(x) round(x, 2)))
    rownames(numeric_df) <- rownames(corrected_transposed_result)
  } else {
    stop("Please define array type. Supports 'chambered', 'segmented', or 'huprot'.")
  }

  if (nrow(numeric_df) == 0) {
    warning("The resulting data frame is empty after processing.")
  }

  return(numeric_df)
}

