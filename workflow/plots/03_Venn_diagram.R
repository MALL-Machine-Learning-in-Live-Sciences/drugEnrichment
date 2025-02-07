library(VennDiagram)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Inputs 
DE_results_1 <- "data/DE_TF_PRISM/DE_TF_PRISM.csv"
DE_results_2 <- "data/DE_TF_SANGER/DE_TF_SANGER.csv"

# Outputs 
figure_path <- "figures/venndiagram/1prism_sanger_venn.png"

# Select a threshold
thr <- 0.65
tumor_type <- "Colorectal_Adenocarcinoma"
moa1 <- "EGFR inhibitor"
moa2 <- "EGFR signaling"

# Filter the results
DE1 <- read.csv(DE_results_1, header = TRUE)
DE1positive <- DE1 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value > thr) %>%
  dplyr::pull(TF)
DE1negative <- DE1 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa1, GSEA.value < (-thr)) %>%
  dplyr::pull(TF)

DE2 <- read.csv(DE_results_2, header = TRUE)
DE2positive <- DE2 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa2, GSEA.value > thr) %>%
  dplyr::pull(TF)
DE2negative <- DE2 %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa2, GSEA.value < (-thr)) %>%
  dplyr::pull(TF)


# Prepare data
venn_data_positive <- list(
  `PRISM` = DE1positive,
  `GDSC` = DE2positive
)

venn_data_negative <- list(
  `PRISM` = DE1negative,
  `GDSC` = DE2negative
)

library(VennDiagram)
library(grid)
library(gridExtra)

# Definir los datos para los diagramas de Venn
venn_data_positive <- list(
  `PRISM` = DE1positive,
  `GDSC` = DE2positive
)

venn_data_negative <- list(
  `PRISM` = DE1negative,
  `GDSC` = DE2negative
)

# Diagrama de Venn positivo
venn_positive <- draw.pairwise.venn(
  area1 = length(venn_data_positive$PRISM),
  area2 = length(venn_data_positive$GDSC),
  cross.area = length(intersect(venn_data_positive$PRISM, venn_data_positive$GDSC)),
  category = c("PRISM", "GDSC"),
  fill = c("#792020dd", "#ce34349e"),
  cat.fontfamily = "Times",   # Cambiar la fuente de las etiquetas de los conjuntos
  cat.cex = 1.5,              # Tamaño de la fuente de las etiquetas
  cex = 1.2,                  # Tamaño del número de intersecciones
  lwd = 0.5
)

# Diagrama de Venn negativo
venn_negative <- draw.pairwise.venn(
  area1 = length(venn_data_negative$PRISM),
  area2 = length(venn_data_negative$GDSC),
  cross.area = length(intersect(venn_data_negative$PRISM, venn_data_negative$GDSC)),
  category = c("PRISM", "GDSC"),
  fill = c("#47a4c3d0", "#86b3c2d0"),
  cat.fontfamily = "Times",   # Cambiar la fuente de las etiquetas de los conjuntos
  cat.cex = 1.5,              # Tamaño de la fuente de las etiquetas
  cex = 1.2,                  # Tamaño del número de intersecciones
  lwd = 0.5
)

# Títulos con fuente Times New Roman
positive_title <- textGrob("ES positivo", gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Times"))
negative_title <- textGrob("ES negativo", gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Times"))

# Combinar gráficos
combined_grob <- arrangeGrob(
  grobs = list(
    positive_title, venn_positive,
    negative_title, venn_negative
  ),
  ncol = 1, 
  heights = c(0.1, 1, 0.1, 1)
)

# Guardar la figura
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = figure_path,
  plot = combined_grob,
  width = 12,
  height = 8,
  dpi = 150
)