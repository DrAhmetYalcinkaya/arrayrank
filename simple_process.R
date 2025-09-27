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
qnorm_df <- multinorm(corr_df, method= "quantile")

# 3. data assembly and storage
keyfile <- readxl::read_xlsx(choose.files())
stored_data <- assemble.output(keyfile, wide_df,corr_df,qnorm_df)

# 4a. detect hits with numerical approach and rank
numerical_results <- detect.hits(
  stored_data$analysis_set_offset,
  controls = "BD",
  examine = "APS1",
  sdmean = 0.5,
  fold_threshold = 10,
  absolute_threshold = 5000
  )
# 4b. detect hits with isolation forest and rank
isolation_results <- isolate.hits(
  stored_data$tidied_analysis_set_offset,
  sample_id_col = "Sample",
  cohort_col = "Group",
  controls = "BD",
  examine = "APS1",
  fold_threshold = 10,
  anomaly_scale = 0.5,
  penalize = -4
)
