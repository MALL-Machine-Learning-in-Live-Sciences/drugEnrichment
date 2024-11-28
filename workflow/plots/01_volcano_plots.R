require(ggpubr)
require(ggrepel)
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)

# Inputs
stats_inputpath <- "data/tf_de_PRISM/gsea_results_tf.csv" 

# Outputs
outputpath <- "figures/volcanoplots/tf_prism_ES65.png" 

# Arguments
alpha <- 0.05
thrFC <- 0.65
dpi <- 400
width <- 20
height <- 10
tumor_type <- "Colorectal_Adenocarcinoma"
moa <- "EGFR inhibitor"
selec_column_name <- "TF" 

# Load and filter
stats <- read.csv(stats_inputpath, header = TRUE)
stats <- stats %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa)


# Plot VolcanoPlot
toplot <-
    stats %>%
    mutate(
        p.value = ifelse(p.value == 0, runif(n(), 0.00001, 0.001), p.value), 
        fill = "ES < 65", 
        log10adjpval = -log10(`p.value`), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` > thrFC, "Positive ES", fill), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` < -thrFC, "Negative ES", fill), 
        toshow = NA, # Inicializa 'toshow' con NA
        toshow = ifelse(`p.value` < alpha & abs(GSEA.value) > thrFC, stats[[selec_column_name]], toshow)
    )

combined_string <- paste(tumor_type, moa, selec_column_name)


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
    geom_text_repel(
        data = toplot[!is.na(toplot$toshow), ], 
        aes(label = toshow),
        nudge_x = ifelse(
            toplot$GSEA.value[!is.na(toplot$toshow)] > 0, 
            abs(toplot$GSEA.value[!is.na(toplot$toshow)]), 
            -abs(toplot$GSEA.value[!is.na(toplot$toshow)])
        ) * 0.05, 
        direction ='y',
        max.overlaps = Inf,
        force = 1, 
        force_pull = 0.15,
        box.padding = 0.15, 
        point.padding = 0.05 
    )+
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