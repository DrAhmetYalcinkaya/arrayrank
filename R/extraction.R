#' Extract Protein Data from Limma-read Array Data (RGList Object)
#'
#' Extracts protein data from different types of arrays and outputs it in either long or wide format.
#'
#' @param data An \code{RGList} object containing microarray data, specifically with `genes` (vector), `targets` (list), and `R` (matrix) representing protein data, target files, and measurement values respectively.
#' @param array_type A character string specifying the type of array used in the analysis. Must be one of "chambered", "segmented", or "huprot".
#' @param format A character string specifying the format of the output. Either "wide" (default) or "long".
#'
#' @return A data frame containing the extracted and arranged data in the specified format. The data frame contains calculated mean values for each group, with columns that vary depending on the specified `array_type` and `format`.
#' @export
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' # Assuming 'raw_data' is an RGList object obtained from read.gpr()
#' Meaning run read.gpr first!
#'
#' # Extract data in long format for chambered arrays
#' long_data <- extraction(data = raw_data, array_type = "chambered", format = "long")
#'
#' # Extract data in wide format for huprot arrays
#' wide_data <- extraction(data = raw_data, array_type = "huprot", format = "wide")
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

  if (format == "long" & array_type != "huprot") {
    mid <- df_out %>%
      dplyr::group_by(array, Block, Name, ID) %>%
      dplyr::summarise(Mvalue = mean(value, na.rm = TRUE)) %>%
      ungroup()
    long_result <- mid %>%
      group_by(array, Block, Name) %>%
      summarise(Mvalue = if(n() > 5) Mvalue = mean(Mvalue, na.rm = T) else Mvalue) %>%
      mutate(Name = make.unique(as.character(Name), sep = ".")) %>%
      ungroup()
    output <- long_result
    message("Data calculated/arranged for chambered/segmented arrays")
    message("Long data created based on format input (user defined: long)")
  } else if (format == "long" & array_type == "huprot") {
    long_result <- df_out %>%
      dplyr::group_by(array, Name, Block, Column) %>%
      dplyr::summarise(Mvalue = mean(value, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Block, -Column) %>%
      dplyr::group_by(array, Name) %>%
      dplyr::summarise(Mvalue = if (n() >= 8) mean(Mvalue, na.rm = TRUE) else Mvalue) %>%
      dplyr::mutate(Name = make.unique(as.character(Name), sep = "."))
    message("Data calculated/arranged for huprot arrays")
    message("Long data created based on format input (user defined: long)")
    output <- long_result
  } else if (format == "wide") {
    if (array_type == "chambered" || array_type == "segmented") {
      mid <- df_out %>%
        dplyr::group_by(array, Block, Name, ID) %>%
        dplyr::summarise(Mvalue = mean(value, na.rm = TRUE)) %>%
        ungroup()
      long_result <- mid %>%
        group_by(array, Block, Name) %>%
        summarise(Mvalue = if(n() > 5) Mvalue = mean(Mvalue, na.rm = T) else Mvalue) %>%
        mutate(Name = make.unique(as.character(Name), sep = ".")) %>%
        ungroup()
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
      long_result <- df_out %>%
        dplyr::group_by(array, Name, Block, Column) %>%
        dplyr::summarise(Mvalue = mean(value, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Block, -Column) %>%
        dplyr::group_by(array, Name) %>%
        dplyr::summarise(Mvalue = if (n() >= 8) mean(Mvalue, na.rm = TRUE) else Mvalue) %>%
        dplyr::mutate(Name = make.unique(as.character(Name), sep = "."))
      wide_result <- tidyr::pivot_wider(long_result, names_from = Name, values_from = Mvalue)
      message("Data calculated/arranged for huprot arrays. Duplicate protein names appended with .n --creating dataframe...")

      transposed_result <- as.data.frame(t(wide_result))
      file_names <- as.character(transposed_result[1, ])
      array <- paste0("Arr", sapply(basename(file_names), function(x) gsub("[^0-9]", "", x)))
      new_colname <- array
      corrected_transposed_result <- as.data.frame(transposed_result[-1, ])
      colnames(corrected_transposed_result) <- new_colname
      row_names <- rownames(transposed_result)
      row_names <- row_names[-1]
      numeric_df <- as.data.frame(lapply(corrected_transposed_result, as.numeric))
      numeric_df <- as.data.frame(apply(numeric_df, 2, function(x) round(x, 2)))
      rownames(numeric_df) <- row_names
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
