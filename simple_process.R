library(devtools)
devtools::install_github("DrAhmetYalcinkaya/arrayrank")
library(arrayrank)

#Example process

# 1. collect gpr data and make dataset
raw_data <- read.gpr(bgcorrect = T)
wide_df <- extraction(raw_data, array_type = "chambered", format = "wide")

# 2. data correction
corr_df <- replace.low(wide_df, offset = 10)
qnorm_df <- multinorm(corr_df, method= "quantile")

# 3. storage and backup
qstore_df <- store.df(qnorm_df)
writexl::write_xlsx(qstore_df, "normed.xlsx")

# 4. read data from backup file
q_data <- readxl::read_xlsx(choose.files())

# 5. read group info from excel (if you want to provide a vector defining the sample groups)
groups <- readxl::read_xlsx(choose.files())
groups <- groups[1,]
# adjust groups definition as required, the above assumes that groups are provided side-by-side
# in the first row of an Excel file used to list the groups of each subject (not the backup data file).

#IMPORTANT: Either add group info into the first column of the backup data (normed.xlsx for this example) manually,
#or create a separate vector/excel file with the group info for each sample. Code accounts for both methods.

# 6. detect proteins with hits and rank them (in this example: the groups, BD and SSc, have been provided as a vector)
q_result <- detect.hits(q_data, controls = "BD", group_vector = groups, examine = "SSc", sdmean = 0.7,
                           fold_threshold = 10, absolute_threshold = 8000)

# If group data were added into the backup data file, then remove "group_vector = groups" before running.
