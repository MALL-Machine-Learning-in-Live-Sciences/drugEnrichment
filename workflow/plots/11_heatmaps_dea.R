library(ComplexHeatmap)
library(grid)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Inputs 
dataPath <- "data/diffexp/DE_fmntop50_filter.csv"
countsPath <- "data/GSE162104_processed/GSE162104_counts.csv"
metadataPath <- "data/GSE162104_processed/GSE162104_metadata.csv"

# Otutputs
heatmap_path <- "figures/adicional/padj2.png"
dir.create(dirname(heatmap_path), recursive = TRUE, showWarnings = FALSE)

# Plotting 
# data <- read.csv(dataPath)
# counts <- read.csv(countsPath, row.names = 1)
# metadata <- read.csv(metadataPath)

data <- data %>%
  dplyr::filter(padj<0.05)

annotation <- data.frame(Group = metadata$Sample_characteristics_ch1.1 )
row.names(annotation) <- metadata$Sample_title
tcounts <- t(counts)
counts_filtered <- tcounts[rownames(tcounts) %in% data$Row.names, ]
#counts_filtered <- counts_filtered[, !(names(counts_filtered) %in% c("Row.names"))]

annotation <- annotation[match(colnames(counts_filtered), rownames(annotation)), , drop = FALSE]
annotation <- annotation %>%
  filter(Group != "cetuximab resistance: nd")

common_samples <- intersect(colnames(counts_filtered), rownames(annotation))
counts_filtered <- counts_filtered[, common_samples, drop = TRUE]

annotation$Group[annotation$Group == "cetuximab resistance: sensitive"] <- "Sensible"
annotation$Group[annotation$Group == "cetuximab resistance: resistant"] <- "Resistente"

ann_colors <- list(
  Group = c(
    "Resistente" = "#8e0000", 
    "Sensible" = "#008118"
  )  
)

col_fun <- colorRamp2(c(-5, 0, 5), c("#0081abd0", "#dddddd", "#9d1a1add"))
col_annotation <- HeatmapAnnotation(df = annotation, col = ann_colors)

library(ggplot2)

counts_scaled <- t(scale(as.matrix(t(counts_filtered)))) 

heatmap_plot <- grid.grabExpr(
  draw(Heatmap(
    as.matrix(counts_scaled),
    name = "Expression",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_title = "Gene Expression Heatmap",
    column_title_gp = gpar(fontfamily = "Times"),
    top_annotation = col_annotation,
    column_names_gp = gpar(fontfamily = "Times"),
    row_names_gp = gpar(fontfamily = "Times"),
    heatmap_legend_param = list(title_gp = gpar(fontfamily = "Times"),
                                labels_gp = gpar(fontfamily = "Times"))
  ))
)

# Guardar con ggsave
ggsave(heatmap_path, plot = heatmap_plot, width = 12, height = 12, dpi = 300)

