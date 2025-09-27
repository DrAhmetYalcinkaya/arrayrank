#' Assemble and Format Protein Data for Multiple Analysis Types
#'
#' This function prepares protein data in multiple formats. It creates wide-format
#' data for Excel, and generates both wide and tidy analysis sets for each of the
#' three input data types (raw, offset-corrected, normalized).
#'
#' @param key A data frame with at least an "Array" column and associated metadata,
#'   such as "Group" and "Sample".
#' @param wide_df Raw protein expression data (proteins as rows, samples as columns).
#' @param corr_df Offset-corrected data (same structure as wide_df).
#' @param norm_df Normalized data (same structure as wide_df).
#' @param file_name Name of the output Excel file.
#'
#' @return A named list of data frames. It includes 'Raw', 'Offset', 'Norm' (for Excel),
#'   as well as wide-format 'analysis_set' and tidy-format 'tidied_analysis_set'
#'   for each of the three data types.
#' @export
#'
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom writexl write_xlsx
#'
assemble.output <- function(key, wide_df, corr_df, norm_df,
                            file_name = "Data_complete.xlsx") {

  # wide formating for excel
  .match_key_wide <- function(key, data) {
    trans_data <- as.data.frame(t(data))
    trans_data <- tibble::rownames_to_column(trans_data, "Array")
    merged <- dplyr::left_join(trans_data, key, by = "Array")

    meta_cols <- setdiff(names(key), "Array")
    meta_rows <- t(merged[, meta_cols, drop = FALSE])
    colnames(meta_rows) <- merged$Array

    expr_rows <- t(merged[, setdiff(names(merged), c("Array", meta_cols)), drop = FALSE])
    colnames(expr_rows) <- merged$Array

    out <- rbind(meta_rows, expr_rows)
    return(as.data.frame(out))
  }

  # standard data set
  .create_tidy_df <- function(key, data) {
    trans_data <- as.data.frame(t(data))
    trans_data <- tibble::rownames_to_column(trans_data, "Array")
    tidy_df <- dplyr::left_join(key, trans_data, by = "Array")
    return(tidy_df)
  }

  # protein names addition
  .add_protein_names <- function(data) {
    protein_names <- rownames(data)
    output <- cbind(Protein = protein_names, data)
    return(output)
  }

  # file outputs (writable stuff)
  stored_raw <- .add_protein_names(.match_key_wide(key, wide_df))
  stored_offset_corr <- .add_protein_names(.match_key_wide(key, corr_df))
  stored_norm <- .add_protein_names(.match_key_wide(key, norm_df))

  # detect.hits-compatible data sets
  analysis_set_raw <- stored_raw[-2, ]
  analysis_set_offset <- stored_offset_corr[-2, ]
  analysis_set_norm <- stored_norm[-2, ]

  # isolate.hits-compatible data sets
  tidied_analysis_set_raw <- .create_tidy_df(key, wide_df)
  tidied_analysis_set_offset <- .create_tidy_df(key, corr_df)
  tidied_analysis_set_norm <- .create_tidy_df(key, norm_df)

  # Combine for in-R output
  stored_data <- list(
    # Ogoes to excel
    Raw_Excel = stored_raw,
    Offset_Excel = stored_offset_corr,
    Norm_Excel = stored_norm,

    # goes to detect.hits
    analysis_set_raw = analysis_set_raw,
    analysis_set_offset = analysis_set_offset,
    analysis_set_norm = analysis_set_norm,

    # goes to isolation forest
    tidied_analysis_set_raw = tidied_analysis_set_raw,
    tidied_analysis_set_offset = tidied_analysis_set_offset,
    tidied_analysis_set_norm = tidied_analysis_set_norm
  )

  writexl::write_xlsx(
    list(Raw = stored_raw, Offset = stored_offset_corr, Norm = stored_norm),
    file_name
  )

  message("Data written to Excel file: ", file_name)
  return(stored_data)
}
