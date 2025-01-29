if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

library(devtools)
remove.packages("arrayrank")
detach("package:arrayrank", unload = TRUE)

devtools::install_github("DrAhmetYalcinkaya/arrayrank")
library(arrayrank)

#Example process

# 1. collect gpr data and make dataset
raw_data <- read.gpr(bgcorrect = T, fdata = "median")
wide_df <- extraction(raw_data, array_type = "huprot", format = "wide")

# 1b. identify discordant duplicates in data (to cross-check with ultimate output)
discordants <- discordant(raw_data, array_type = "huprot", fold = 4, abs = 2000)

# 2. data correction
corr_df <- replace.low(wide_df, offset = 10)
qnorm_df <- multinorm(corr_df, method= "ig")

# 3. storage and backup
qstore_df <- store.df(qnorm_df)
writexl::write_xlsx(qstore_df, "your_storage_file.xlsx")

# 4. read data from backup file
q_data <- readxl::read_xlsx(choose.files())

# 5. read group info from excel (if you want to provide a vector defining the sample groups)
groups <- readxl::read_xlsx(choose.files())
groupsv <- t(groups[,2]) ### note: read/create/edit the character vector based on your own requirements!!!

#IMPORTANT: Either add group info into the first row of the backup data (your_file_name.xlsx for this example) manually,
#or create a separate vector/excel file with the group info for each sample. Code accounts for both methods.

# 6. detect proteins with hits and rank them (in this example: the groups have been provided as a vector)
q_result <- detect.hits(q_data, controls = "your_control_group", group_vector = groupsv,
                        examine = "your_study_group", sdmean = 0.5, fold_threshold = 5, absolute_threshold = 5000)

# If group data were added into the backup data file (so there is no 'groups' vector), remove "group_vector = groups" before running.
