#' Assemble Protein Data with Key Matching and Excel Export
#'
#' Matches sample keys to multiple datasets (raw, offset-corrected, normalized),
#' adds protein names, and writes each dataset as a separate sheet in a single Excel file.
#'
#' @param key A data frame with at least an "Array" column and associated metadata. Preferably, the key file will contain "Array", "Group", and "Sample" columns in that order).
#' @param wide_df Raw protein expression data (proteins as rows, samples as columns, based on the read output).
#' @param corr_df Offset-corrected data (same structure as wide_df).
#' @param norm_df Normalized data (same structure as wide_df).
#' @param file_name Name of the output Excel file (defaults to "Data_complete.xlsx").
#'
#' @return A named list of data frames (Raw, Offset, Norm) that are also written to Excel in your working directory.
#' @export
#'
#' @importFrom dplyr left_join relocate
#' @importFrom tibble rownames_to_column
#' @importFrom writexl write_xlsx
#'
#' @examples
#' \dontrun{
#' assembled_data <- assemble.output(
#'   key = sample_key,
#'   wide_df = raw_data,
#'   corr_df = corrected_data,
#'   norm_df = normalized_data,
#'   file_name = "MyData.xlsx"
#' )
#' }
assemble.output <- function(key, wide_df, corr_df, norm_df,
                            file_name = "Data_complete.xlsx") {

  .match_key <- function(key, data) {
    # Transpose so samples are rows
    trans_data <- as.data.frame(t(data))
    trans_data <- tibble::rownames_to_column(trans_data, "Array")

    # Find unmatched samples before join
    unmatched_in_data <- setdiff(trans_data$Array, key$Array)
    unmatched_in_key <- setdiff(key$Array, trans_data$Array)

    if (length(unmatched_in_data) > 0) {
      warning("!!! IMPORTANT WARNING: Arrays in data not found in key, meaning that group and sample names were not found for the following arrays:\n", paste(unmatched_in_data, collapse = ", "))
    }
    if (length(unmatched_in_key) > 0) {
      warning("Secondary warning: ", length(unmatched_in_key)," arrays in key not found in data, which may be OK if the key file contains keys for other arrays not being examined.")
    }

    # Join metadata
    merged <- dplyr::left_join(trans_data, key, by = "Array")

    # Separate metadata and protein expression
    meta_cols <- setdiff(names(key), "Array")
    meta_rows <- t(merged[, meta_cols, drop = FALSE])
    colnames(meta_rows) <- merged$Array

    expr_rows <- t(merged[, setdiff(names(merged), c("Array", meta_cols)), drop = FALSE])
    colnames(expr_rows) <- merged$Array

    out <- rbind(meta_rows, expr_rows)
    return(as.data.frame(out))
  }


  # --- Internal: add protein names ---
  .add_protein_names <- function(data) {
    protein_names <- rownames(data)
    output <- cbind(Protein = protein_names, data)
    return(output)
  }

  # --- Process all datasets ---
  stored_raw <- .add_protein_names(.match_key(key, wide_df))
  stored_offset_corr <- .add_protein_names(.match_key(key, corr_df))
  stored_norm <- .add_protein_names(.match_key(key, norm_df))

  stored_data <- list(
    Raw = stored_raw,
    Offset = stored_offset_corr,
    Norm = stored_norm
  )

  # --- Export to Excel ---
  writexl::write_xlsx(stored_data, file_name)

  message("Data written to Excel file in your working directory: ", file_name)
  return(stored_data)
}
