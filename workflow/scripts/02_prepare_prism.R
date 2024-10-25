require(dplyr)

# Inputs 
prism_inputpath <- "extdata/PRISM_19Q4/secondary-screen-dose-response-curve-parameters.csv"

# CCLE metadata
ccle_meta <- readRDS("data/ccle_clinical.rds")

# CCLE counts
ccle_counts <- readRDS("data/ccle_counts.rds")

# PRISM repurposing dataset
# downloaded from https://depmap.org/portal/download/all/
# PRISM Repurposing 19Q4

prism2 <- read.table(
    prism_inputpath,
    header = TRUE,
    sep = ",",
    quote = '"',
    fill = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
    )

prism2 <- prism2[, c("broad_id", "depmap_id", "auc", "moa", "name", "screen_id")]

prism2 <- reshape2::dcast(
    dat = prism2,
    formula = screen_id + name + broad_id + moa ~ depmap_id,
    fun.aggregate = sum,
    value.var = "auc"
    )

prism2[prism2 == 0] <- NA

# remove duplicates
prism2 <-
    prism2 %>%
    group_by(name) %>%
    arrange(factor(screen_id, levels = c("MTS010", "MTS006", "MTS005", "HTS002"))) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(-c(screen_id, name))

prism2_m <- as.matrix(prism2[, -c(1, 2)])
rownames(prism2_m) <- prism2$broad_id
prism2_m <- t(prism2_m)

# Match datasets
cell_lines <- intersect(rownames(ccle_counts), rownames(prism2_m))
ccle_counts <- ccle_counts[cell_lines, ]
prism <- prism2_m[cell_lines, ]

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
