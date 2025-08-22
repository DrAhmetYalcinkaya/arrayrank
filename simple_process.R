if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

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

# 3. data assembly and storage
keyfile <- readxl::read_xlsx(choose.files())
stored_data <- assemble.output(keyfile, wide_df,corr_df,qnorm_df)

# 4. selection of dataset to be analyzed and removal of sample ID row
q_data <- stored_data[["Norm"]][-2,]

# 5. detect proteins with hits and rank them
q_result <- detect.hits(q_data, controls = "BD", examine = "APS1", sdmean = 0.5, fold_threshold = 5, absolute_threshold = 5000)
