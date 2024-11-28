library(ggvenn)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Inputs 
DE_results_1 <- "data/tf_de_PRISM/gsea_results_tf.csv"
DE_results_2 <- "data/tf_de_SANGER/gsea_results_tf_SANGER.csv"

# Outputs 
figure_path <- "figures/venndiagram/prism_sanger_venn.png"

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
  `TF PRISM` = DE1positive,
  `TF SANGER` = DE2positive
)

venn_data_negative <- list(
  `TF PRISM` = DE1negative,
  `TF SANGER` = DE2negative
)

# Build venn diagrams
venn_positive <- ggvenn(
  venn_data_positive,
  fill_color = c("blue", "green"),
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5
)

venn_negative <- ggvenn(
  venn_data_negative,
  fill_color = c("red", "orange"),
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 5
)

positive_grob <- ggplotGrob(venn_positive)
negative_grob <- ggplotGrob(venn_negative)

positive_title <- textGrob("Positive ES", gp = gpar(fontsize = 16, fontface = "bold"))
negative_title <- textGrob("Negative ES", gp = gpar(fontsize = 16, fontface = "bold"))

combined_grob <- arrangeGrob(
  grobs = list(
    positive_title, positive_grob,
    negative_title, negative_grob
  ),
  ncol = 1, 
  heights = c(0.1, 1, 0.1, 1) 
)

dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = figure_path,
  plot = combined_grob,
  width = 12,
  height = 8,
  dpi = 150
)
