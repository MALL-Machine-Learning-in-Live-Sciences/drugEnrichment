require(dplyr)
require(data.table)

# CCLE counts
# downloaded from: https://depmap.org/portal/download/all/
# DepMap Public 23Q2
metadata <- read.csv(
  "extdata/DepMap_23Q2/Model.csv", 
  header = TRUE
  )

counts <- data.table::fread(
  "extdata/DepMap_23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
  header = TRUE
) %>%
  tibble::column_to_rownames("V1")

colnames(counts) <- sapply(strsplit(colnames(counts), "..", fixed = T), "[", 1)
counts <- counts[match(metadata$ModelID, rownames(counts)), ]

saveRDS(counts, file = "data/ccle_brca_counts.rds")
saveRDS(metadata, file = "data/ccle_brca_clinical.rds")



