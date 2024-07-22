require(dplyr)

# CCLE from BRCA 
ccle_brca <- readRDS("data/ccle_brca_clinical.rds")

# PRISM repurposing dataset
# downloaded from https://depmap.org/portal/download/all/
# PRISM Repurposing 19Q4
dose_response <- data.table::fread(
  "extdata/PRISM_19Q4/secondary-screen-dose-response-curve-parameters.csv",
  header = TRUE)

dose_response_f <-
  dose_response %>%
  filter(
    depmap_id %in% ccle_brca$ModelID
  ) %>% 
  dplyr::select(c(broad_id, depmap_id, name, auc)) %>% 
  tidyr::pivot_wider(.,
    names_from = "depmap_id",
    values_from = "auc",
    values_fn = {mean}
  ) %>% 
  tibble::column_to_rownames(var = "broad_id")

treatment_info <- read.csv(
  "extdata/PRISM_19Q4/secondary-screen-replicate-treatment-info.csv",
  header = T
)

saveRDS(treatment_info, file = "data/treatment_info.rds")
saveRDS(dose_response_f, file = "data/drug_response.rds")
