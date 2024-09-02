#' Extract Data from RGList Object
#'
#' This function extracts protein data, array information, and intensity values from an `RGList` object (typically produced by the `limma` package) and returns a combined data frame. Based on user input, the dataframe can be in wide or long format.
#'
#' @param data An `RGList` object containing gene expression data, typically produced by the `limma` package's `read.maimages` function. The object should include `genes`, `targets`, and `R` components.
#' @param array_type A character string specifying the type of array to process. Options are `"chambered"`, `"segmented"`, or `"huprot"`.
#' @param format A character string specifying the type of dataframe to be used for the output. Options are `"wide"` or `"long"`.
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
extraction <- function(data, array_type, format = "wide") {
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

  # Validate array_type parameter
  if (!(array_type %in% c("chambered", "segmented", "huprot"))) {
    stop("Invalid array_type. Please specify 'chambered', 'segmented', or 'huprot'.")
  }
  # Convert data to data.table format
  df <- data.table::setDT(df_out)
  if (format == "long" & array_type != "huprot") {
    long_result <- df[, .(Mvalue = mean(.SD$value, na.rm = TRUE)), by = .(array, Block, Name)]
    output <- long_result
    message("Data calculated/arranged for chambered/segmented arrays")
    message("Long data created based on format input (user defined: long)")
  } else if (format == "long" & array_type == "huprot") {
    long_result <- df[, .(Mvalue = mean(.SD$value, na.rm = TRUE)), by = .(array, Name)]
    output <- long_result
    message("Data calculated/arranged for huprot arrays")
    message("Long data created based on format input (user defined: long)")
  } else if (format == "wide") {
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
      new_colname <- paste0(array, "Block", transposed_result[2,])
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

    colnames(numeric_df) <- gsub("\\.", "0", colnames(numeric_df))
    output <- numeric_df
    message("Wide data created based on format input (user defined: wide)")
  }

  return(output)
}
