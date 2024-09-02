library(arrayrank)
library(readxl)

#Example process

# 1. collect gpr data and make dataset
files <- read.gpr()
bg_subbed <- subtract.bg(files)
long_format <- extraction(bg_subbed)
orig_df <- make.df(long_format, "huprot")

# 2. data correction
df <- replace.low(orig_df, offset = 10)
df_qnorm <- multinorm(df, method= "quantile")

# 3. storage and backup
df_qstore <- store.df(df_qnorm)
writexl::write_xlsx(df_qstore, "normed.xlsx")

# 4. read data from backup file
q_data <- readxl::read_xlsx(choose.files())

# 5. read group info from excel
groups <- read_xlsx(choose.files())
groups <- groups[1,-1]

# 6. detect proteins with hits and rank them
q_result <- detect.hits(q_data, controls = "BD", group_vector = groups, examine = "SSc", sdmean = 0.7,
                           fold_threshold = 10, absolute_threshold = 8000)
