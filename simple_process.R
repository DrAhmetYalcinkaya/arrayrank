library(devtools)
remove.packages("arrayrank")
detach("package:arrayrank", unload = TRUE)

devtools::install_github("DrAhmetYalcinkaya/arrayrank")
library(arrayrank)

#Example process

# 1. collect gpr data and make dataset
raw_data <- read.gpr(bgcorrect = T, fdata = "median")
wide_df <- extraction(raw_data, array_type = "chambered", format = "wide")

# 1b. identify discordant duplicates in data (to cross-check with ultimate output)
discordants <- discordant(raw_data, array_type = "chambered", fold = 2, abs = 1000)

# 2. data correction
corr_df <- replace.low(wide_df, offset = 10)
qnorm_df <- multinorm(corr_df, method= "ig")

# 3. storage and backup
qstore_df <- store.df(qnorm_df)
writexl::write_xlsx(qstore_df, "your_file_name.xlsx")

# 4. read data from backup file
q_data <- readxl::read_xlsx(choose.files())

# 5. read group info from excel (if you want to provide a vector defining the sample groups)
groups <- readxl::read_xlsx(choose.files())
groupsv <- t(groups[,2]) ### note: read/create/edit the character vector based on your own requirements!!!
# adjust groups definition as required.

#IMPORTANT: Either add group info into the first column of the backup data (normed.xlsx for this example) manually,
#or create a separate vector/excel file with the group info for each sample. Code accounts for both methods.

# 6. detect proteins with hits and rank them (in this example: the groups, BD and SSc, have been provided as a vector)
q_result1 <- detect.hits(q_data, controls = "BD", group_vector = groupsv, examine = "Pan", sdmean = 0.5,
                           fold_threshold = 5, absolute_threshold = 5000)

# If group data were added into the backup data file (so there is no 'groups' vector), remove "group_vector = groups" before running.
