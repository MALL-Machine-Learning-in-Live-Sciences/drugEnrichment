library(ComplexHeatmap)
library(grid)

library(ComplexHeatmap)
library(circlize)

# Inputs 
dataPath <- "data/diffexp/DE_fmntop50_filter.csv"
countsPath <- "data/diffexp/ddseq_norm_counts.csv"
metadataPath <- "data/GSE162104_processed/GSE162104_metadata.csv"

# Otutputs
heatmap_path <- "figures/heatmap/heatmap_top 50.png"
dir.create(dirname(heatmap_path), recursive = TRUE, showWarnings = FALSE)

# Plotting 
data <- read.csv(dataPath)
counts <- read.csv(countsPath, row.names = 1)
metadata <- read.csv(metadataPath)

annotation <- data.frame(Group = metadata$Sample_characteristics_ch1.1 )
row.names(annotation) <- metadata$Sample_title

counts_filtered <- counts[rownames(counts) %in% data$Row.names, ]
counts_filtered <- counts_filtered[, !(names(counts_filtered) %in% c("Row.names"))]

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

col_fun <- colorRamp2(c(-2, 0, 2), c("#0081abd0", "#dddddd", "#9d1a1add"))
col_annotation <- HeatmapAnnotation(df = annotation, col = ann_colors)

library(ggplot2)

counts_scaled <- t(scale(t(counts_filtered))) 

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
ggsave(heatmap_path, plot = heatmap_plot, width = 8, height = 8, dpi = 300)

