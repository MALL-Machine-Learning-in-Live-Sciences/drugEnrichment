library(ComplexHeatmap)
library(grid)
library(circlize)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)


# Inputs 
DEnmfPath <- "data/DE_NMF_SANGER/DE_NMF_SANGER.csv"
HmatrixPath <- "data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv"

# Otutputs
heatmap_path <- "figures/heatmapJaccard/top100heatmap_fmn.png"
dir.create(dirname(heatmap_path), recursive = TRUE, showWarnings = FALSE)

# Plotting 
data <- read.csv(DEnmfPath)

tumor_type <- "Colorectal_Adenocarcinoma"
moa <- "EGFR signaling"

filterData <- data %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa, NMF_components==40, abs(GSEA.value)>0.7, p.value<0.05)


Hmatrix <- read.csv(HmatrixPath, row.names = 1)

factors_to_keep <- unique(filterData$Factors)
Hmatrix_filtered <- Hmatrix[rownames(Hmatrix) %in% factors_to_keep, ]


Hmatrix_df <- as.data.frame(Hmatrix_filtered)
Hmatrix_df$Factor <- rownames(Hmatrix_df)
Hmatrix_long <- melt(Hmatrix_df, id.vars = "Factor", variable.name = "Gene", value.name = "weight")


Hmatrix_top10 <- Hmatrix_long %>%
  group_by(Factor) %>%  
  mutate(percentile_90 = quantile(weight, 0.9, na.rm = TRUE)) %>%  
  filter(weight >= percentile_90)  

Hmatrix_top50 <- Hmatrix_long %>%
  group_by(Factor) %>%  
  slice_max(order_by = weight, n = 50, with_ties = FALSE)  


# Crear la matriz gene_matrix
gene_matrix <- Hmatrix_top50 %>%
  dplyr::select(Factor, Gene) %>%
  distinct() %>%
  mutate(Present = 1) %>%  
  pivot_wider(names_from = Gene, values_from = Present, values_fill = list(Present = 0))

factor_sets <- list()
for (i in 1:nrow(gene_matrix)) {
  factor_name <- gene_matrix$Factor[i]
  factor_sets[[factor_name]] <- colnames(gene_matrix)[-1][gene_matrix[i, -1] == 1]  
}

jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))  
  union <- length(union(set1, set2)) 
  return(intersection / union)  
}

factor_names <- names(factor_sets)
jaccard_matrix <- matrix(0, nrow = length(factor_names), ncol = length(factor_names), 
                         dimnames = list(factor_names, factor_names))

for (i in 1:length(factor_names)) {
  for (j in i:length(factor_names)) {
    jaccard_value <- jaccard_index(factor_sets[[factor_names[i]]], factor_sets[[factor_names[j]]])
    jaccard_matrix[i, j] <- jaccard_value
    jaccard_matrix[j, i] <- jaccard_value  
  }
}

filterData <- filterData[match(rownames(jaccard_matrix), filterData$Factors), ]

es_annotation <- ifelse(filterData$GSEA.value > 0, "ES > 0", "ES < 0")
names(es_annotation) <- filterData$Factors

row_ha <- rowAnnotation(ES = es_annotation,
                        col = list(ES = c("ES > 0" = "#792020dd", "ES < 0" = "#47a4c3d0")))

color_scale <- colorRamp2(c(0, 0.5, 1), c("white", "#523e72", "#302543"))  

heatmap_plot <- grid.grabExpr(
  draw(Heatmap(
    jaccard_matrix,
    name = "Jaccard",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    right_annotation = row_ha,
    col = color_scale,
    row_names_gp = gpar(fontfamily = "Times New Roman"),
    column_names_gp = gpar(fontfamily = "Times New Roman"),
    heatmap_legend_param = list(
      title_gp = gpar(fontfamily = "Times New Roman"),
      labels_gp = gpar(fontfamily = "Times New Roman"),
      legend_direction = "horizontal",  # Orientación horizontal para la leyenda
      legend_position = "lower"    
    )
  ))
)

ggsave(heatmap_path, plot = heatmap_plot, width = 5, height = 4, dpi = 1500)
