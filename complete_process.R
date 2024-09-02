library(roxygen2)
library(devtools)
roxygenize()
build()
check()

remove.packages("arrayrank")
devtools::install_github("DrAhmetYalcinkaya/arrayrank")
library(arrayrank)


files <- read.gpr()
bg_subbed <- subtract.bg(files)
long_format <- extraction(bg_subbed)
orig_df <- make.df(long_format, "huprot")
df <- replace.low(orig_df, offset = 10)

dfs <- store.df(df)
writexl::write_xlsx(dfs, "non.xlsx")

df_qnorm <- multinorm(df, method= "quantile")
df_qstore <- store.df(df_qnorm)
writexl::write_xlsx(df_qstore, "normed.xlsx")

which(rownames(df) == "Anti-human IgG 1.5625ng/ul")

df_vnorm <- multinorm(df, method = "protein", protein = 1198)
df_vstore <- store.df(df_vnorm)
writexl::write_xlsx(df_vstore, "var_normed.xlsx")

library(readxl)
base_data <- readxl::read_xlsx(choose.files())
q_data <- readxl::read_xlsx(choose.files())
v_data <- readxl::read_xlsx(choose.files())

groups <- read_xlsx(choose.files())
groups <- groups[1,-1]

base_result <- detect.hits(base_data, controls = "BD", group_vector = groups, examine = "SSc", sdmean = 0.7,
                        fold_threshold = 10, absolute_threshold = 8000)

q_result <- detect.hits(q_data, controls = "BD", group_vector = groups, examine = "SSc", sdmean = 0.7,
                           fold_threshold = 10, absolute_threshold = 8000)

v_result <- detect.hits(v_data, controls = "BD", group_vector = groups, examine = "SSc", sdmean = 0.7,
                           fold_threshold = 10, absolute_threshold = 8000)
