# Load necessary library
library(data.table)

# Specify the directory containing the .gpr files
input_directory <- "path_to_your_gpr_files"

# Specify the directory to save the modified files
output_directory <- "path_to_save_modified_gpr_files"

# Create the output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# Get a list of all .gpr files in the input directory
gpr_files <- list.files(input_directory, pattern = "\\.gpr$", full.names = TRUE)

# Function to remove content after '=' for specific keys
remove_content_after_equal <- function(file) {
  # Read all lines from the file, trying to fix any encoding issues
  lines <- readLines(file, warn = FALSE, encoding = "latin1")

  # Convert the lines to UTF-8 encoding, skipping invalid characters
  lines <- iconv(lines, from = "latin1", to = "UTF-8", sub = "")

  # Replace content after "=" for the specified fields
  lines <- gsub("ImageFiles=.*", "ImageFiles=", lines)
  lines <- gsub("JpegImage=.*", "JpegImage=", lines)

  # Construct the new file path in the output directory
  output_file <- file.path(output_directory, basename(file))

  # Write the modified lines to the new file
  writeLines(lines, output_file)
}

# Apply the function to all .gpr files
lapply(gpr_files, remove_content_after_equal)

print("Content after '=' removed in specified cells and saved to new files.")
