# arrayrank

Extracts .gpr array data and detects possible elevated responses ('hits') in selected subject groups relative to controls. Can also be used to detect hits in other forms of multiplex data as long as extraction is done properly.

To load into R environment:

1.  devtools::install_github("DrAhmetYalcinkaya/arrayrank")

2.  library(arrayrank)

You will also need Limma from Bioconductor:

1.  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

2.  BiocManager::install("limma")

FAQ: 1) When I run read.gpr, I get errors about file names, directories, and extensions. Why?

--When Genepix software outputs gpr files, sometimes the 'linked' file directories are very long or the links are nonsensical. This appears to be a problem for the limma package. Since the initial read of the arrayrank package relies entirely upon limma, this error is replicated here. Please use the "correctgprstructure.R" helper function in the root directory of arrayrank. This will handle the most common problems (at 2 sites) in the file(s).
