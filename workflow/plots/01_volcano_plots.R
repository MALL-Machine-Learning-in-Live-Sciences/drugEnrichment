require(ggpubr)
require(ggrepel)
library(org.Hs.eg.db)
library(dplyr)
library(AnnotationDbi)

# Inputs
stats_inputpath <- "data/DE_TF_SANGER_2/DE_TF_SANGER.csv" 

# Outputs
outputpath <- "figures/Z/path_volcano.png" 

# Arguments
alpha <- 0.05
thrFC <- 0.01
dpi <- 1500
width <- 6
height <- 4
tumor_type <- "Colorectal_Adenocarcinoma"
moa <- "EGFR signaling"
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
        fill = "|ES| < 0.7", 
        log10adjpval = -log10(`p.value`), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` > thrFC, "ES positivo", fill), 
        fill = ifelse(`p.value` < alpha & `GSEA.value` < -thrFC, "ES negativo", fill), 
        toshow = NA, # Inicializa 'toshow' con NA
        toshow = ifelse(`p.value` < alpha & abs(GSEA.value) > thrFC, stats[[selec_column_name]], toshow)
    )

combined_string <- paste(tumor_type, moa, selec_column_name)


volcano <- ggscatter(
    toplot,
    x = "GSEA.value",
    y = "log10adjpval",
    color = "fill",
    size = 3,
    label = NULL,
    repel = TRUE,
    title = combined_string,
    ) +
    geom_text_repel(
        data = toplot[!is.na(toplot$toshow), ], 
        aes(label = toshow),
        color = "black", 
        family = "Times New Roman", # Oculta las etiquetas
        nudge_x = ifelse(
            toplot$GSEA.value[!is.na(toplot$toshow)] > 0, 
            abs(toplot$GSEA.value[!is.na(toplot$toshow)]), 
            -abs(toplot$GSEA.value[!is.na(toplot$toshow)])
        ) * 0.26, 
        direction ='y',
        max.overlaps = Inf,
        force = 1, 
        force_pull = 0.35,
        box.padding = 0.25, 
        point.padding = 0.15 
    )+
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth=0.5) +
    geom_vline(xintercept = c(-thrFC, thrFC), linetype = "dashed", linewidth=0.5) +
    scale_color_manual(
        values = c(
            "ES < 65" = "gray",        
            "ES positivo" = "#792020dd", 
            "ES negativo" = "#47a4c3d0"     
        ),
    name=""
    ) + 
    labs(x = "ES", y = "Log(pvalue)")+
    theme(
        axis.line = element_line(linewidth = 0.5),  # Grosor de las líneas de los ejes
        axis.ticks = element_line(linewidth = 0.5),  # Grosor de las marcas de los ejes
        axis.title = element_text(size = 10, family = "Times New Roman"),  # Fuente Times New Roman para los títulos de los ejes
        axis.text = element_text(size = 8, family = "Times New Roman"),  # Fuente Times New Roman para el texto de los ejes
        plot.title = element_text(size = 12, family = "Times New Roman"),  # Fuente Times New Roman para el título del gráfico
        legend.text = element_text(family = "Times New Roman"),  # Fuente Times New Roman para los textos de la leyenda
        legend.title = element_text(family = "Times New Roman")  # Fuente Times New Roman para el título de la leyenda
    )
     
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