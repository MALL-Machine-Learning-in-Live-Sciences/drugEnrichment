library(ComplexHeatmap)
library(grid)
library(circlize)
library(dplyr)
library(reshape2)


# Inputs 
dataPath <- "data/diffexp/DEA_cohort.csv"
countsPath <- "data/diffexp/ddseq_norm_counts.csv"
metadataPath <- "data/GSE162104_processed/GSE162104_metadata.csv"
dseaResultPAth <- "data/DE_TF_PRISM/DE_TF_PRISM.csv"
dseaResultPAthsanger <- "data/DE_TF_SANGER/DE_TF_SANGER.csv"
netpath <- "data/collectri.csv"
# Otutputs
heatmap_path <- "figures/heatmap/heatmap_tfs.png"
dir.create(dirname(heatmap_path), recursive = TRUE, showWarnings = FALSE)

# Plotting 
data <- read.csv(dataPath)
counts <- read.csv(countsPath, row.names = 1)
metadata <- read.csv(metadataPath)
net <- read.csv(netpath, row.names = 1, check.names=FALSE)

annotation <- data.frame(Group = metadata$Sample_characteristics_ch1.1 )
row.names(annotation) <- metadata$Sample_title

ann_colors <- list(
  Group = c(
    "Resistente" = "#8e0000", 
    "Sensible" = "#008118"
  )  
)

tumor_type <- "Colorectal_Adenocarcinoma"
moa1 <- "EGFR inhibitor"
moa2 <- "EGFR signaling"
thr <- 0.65

DE1 <- read.csv(dseaResultPAth, header = TRUE)
DE1positive <- DE1 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value > thr) %>%
  dplyr::pull(TF)
DE1negative <- DE1 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value < (-thr)) %>%
  dplyr::pull(TF)

DE2 <- read.csv(dseaResultPAth, header = TRUE)
DE2positive <- DE2 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value > thr) %>%
  dplyr::pull(TF)
DE2negative <- DE2 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value < (-thr)) %>%
  dplyr::pull(TF)

positive_intersection <- intersect(DE1positive, DE2positive)
negative_intersection <- intersect(DE1negative, DE2negative)

intersection_union <- union(positive_intersection, negative_intersection)
net_filtered <- net[net$source %in% intersection_union, ]

net_filtered$Anotacion <- ifelse(net_filtered$source %in% positive_intersection, "Positivo", 
                                 ifelse(net_filtered$source %in% negative_intersection, "Negativo", NA))

data_filtered <- data %>%
  dplyr::filter(abs(padj) < 0.05)
net_filtered_selected <- net_filtered[net_filtered$target %in% data_filtered$Row.names, ]

net_filtered_selected$GSEA_value <- NA

net_filtered_selected$GSEA_value[net_filtered_selected$source %in% positive_intersection] <- 
  DE1$GSEA.value[match(net_filtered_selected$source[net_filtered_selected$source %in% positive_intersection], DE1$TF)]

net_filtered_selected$GSEA_value[net_filtered_selected$source %in% negative_intersection] <- 
  DE1$GSEA.value[match(net_filtered_selected$source[net_filtered_selected$source %in% negative_intersection], DE1$TF)]

net_filtered_selected_max <- net_filtered_selected %>%
  mutate(abs_GSEA_value = abs(GSEA_value)) %>%  
  group_by(target) %>%  #
  slice(which.max(abs_GSEA_value)) %>%  
  ungroup() %>%  
  select(-abs_GSEA_value)

row_annotation <- rowAnnotation(
  annotation = net_filtered_selected_max$Anotacion, 
  col = list(annotation = c("Positivo" = "#792020dd", "Negativo" = "#47a4c3d0")) 
)


counts_filtered <- counts[rownames(counts) %in% net_filtered_selected_max$target, ]
counts_filtered <- counts_filtered[, !(names(counts_filtered) %in% c("Row.names"))]

annotation <- annotation[match(colnames(counts_filtered), rownames(annotation)), , drop = FALSE]
annotation <- annotation %>%
  filter(Group != "cetuximab resistance: nd")

common_samples <- intersect(colnames(counts_filtered), rownames(annotation))
counts_filtered <- counts_filtered[, common_samples, drop = TRUE]

annotation$Group[annotation$Group == "cetuximab resistance: sensitive"] <- "Sensible"
annotation$Group[annotation$Group == "cetuximab resistance: resistant"] <- "Resistente"

col_fun <- colorRamp2(c(-2, 0, 2), c("#0081abd0", "#dddddd", "#9d1a1add"))
col_annotation <- HeatmapAnnotation(df = annotation, col = ann_colors)

library(ggplot2)

counts_scaled <- t(scale(t(as.matrix(counts_filtered)))) 


heatmap_plot <- grid.grabExpr(
  draw(Heatmap(
    as.matrix(counts_scaled),
    name = "Expression",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = "Gene Expression Heatmap",
    top_annotation = col_annotation,
    left_annotation = row_annotation
  ))
)

# # # Guardar con ggsave
ggsave(heatmap_path, plot = heatmap_plot, width = 8, height = 8, dpi = 300)

