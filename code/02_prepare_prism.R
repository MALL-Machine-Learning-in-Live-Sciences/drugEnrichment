require(dplyr)

# CCLE metadata
ccle_meta <- readRDS("data/ccle_clinical.rds")

# CCLE counts
ccle_counts <- readRDS("data/ccle_counts.rds")

# PRISM repurposing dataset
# downloaded from https://depmap.org/portal/download/all/
# PRISM Repurposing 19Q4
prism <- data.table::fread(
  "extdata/PRISM_19Q4/secondary-screen-dose-response-curve-parameters.csv",
  header = TRUE) %>%
  filter(
    depmap_id %in% ccle_meta$ModelID
  ) %>% 
  dplyr::select(c(broad_id, depmap_id, auc)) %>% 
  tidyr::pivot_wider(.,
    names_from = "depmap_id",
    values_from = "auc",
    values_fn = {mean}
  ) %>% 
  tibble::column_to_rownames(var = "broad_id") %>%
  t() %>% as.data.frame()

# Match datasets
cell_lines <- intersect(rownames(ccle_counts), rownames(prism))
ccle_counts <- ccle_counts[cell_lines, ]
prism <- prism[cell_lines, ]

# Define drug PKN
drug_pkn <- read.csv(
  "extdata/PRISM_19Q4/secondary-screen-replicate-treatment-info.csv",
  header = T
) %>%
  filter(!is.na(moa)) %>%
  tidyr::separate_rows(moa, sep = ",\\s*") %>%
  select(moa, broad_id) %>%
  distinct(moa, broad_id)

write.csv(ccle_counts, file = "data/ccle_counts_matched.csv")
write.csv(prism, file = "data/prism_matched.csv")
write.csv(drug_pkn, file = "data/drug_pkn.csv")
