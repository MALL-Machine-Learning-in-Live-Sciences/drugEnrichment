require(ggpubr)
require(ggrepel)
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)

# Inputs
stats_inputpath <- "data/tf_gsea/gsea_results_tf.csv" #For TF_results
#stats_inputpath <-"data/nmf_gsea/gsea_results_nmf.csv" #For NMF results

# Outputs
outputpath <- "figures/volcanoplots/tf_volcano.png" #For TF_results
#outputpath <- "figures/volcanoplots/nmf_volcano.png" #For NMF results

# Arguments
alpha <- 0.05
thrFC <- 0.5
dpi <- 400
width <- 20
height <- 10
tumor_type <- "Colorectal_Adenocarcinoma"
moa <- "EGFR inhibitor"
selec_column_name <- "TF" #For TF_results
#selec_column_name <- "Factor" #For NMF results

# Load and filter
stats <- read.csv(stats_inputpath, header = TRUE)
stats <- stats %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa)


# Plot VolcanoPlot
toplot <-
    stats %>%
    mutate(
        fill = "ns", 
        log10adjpval = -log10(`p.value`), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` > thrFC, "up_cluster1", fill), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` < -thrFC, "up_cluster2", fill), 
        toshow = NA, # Inicializa 'toshow' con NA
        toshow = ifelse(`p.value` < alpha & abs(GSEA.value) > thrFC, selec_column_name, toshow)
    )

combined_string <- paste(tumor_type, moa)


volcano <- ggscatter(
    toplot,
    x = "GSEA.value",
    y = "log10adjpval",
    color = "fill",
    size = 5,
    label = NULL,
    repel = TRUE,
    title = combined_string
    ) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    geom_vline(xintercept = c(-thrFC, thrFC), linetype = "dashed") +
    scale_color_brewer(palette = "Set1")

dir.create(dirname(outputpath), recursive = TRUE, showWarnings = FALSE)

ggsave(
    volcano,
    device = "png",
    path = dirname(outputpath),
    filename = basename(outputpath),
    dpi = dpi,
    width = width,
    height = height
)