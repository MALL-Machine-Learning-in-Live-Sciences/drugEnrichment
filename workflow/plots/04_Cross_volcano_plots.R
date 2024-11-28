require(ggpubr)
require(ggrepel)
library(dplyr)
library(gridExtra)

# Inputs
stats_inputpath_1 <- "data/tf_de_PRISM/gsea_results_tf.csv"
stats_inputpath_2 <- "data/tf_de_SANGER/gsea_results_tf_SANGER.csv"

# Outputs
outputpath <- "figures/volcanoplots/cross_volcanos.png"

# Arguments
alpha <- 0.05
thr <- 0.65
dpi <- 400
width <- 20
height <- 10
tumor_type <- "Colorectal_Adenocarcinoma"
moa_1 <- "EGFR inhibitor"
moa_2 <- "EGFR signaling"
selec_column_name <- "TF"

# Load and filter datasets
data_1 <- read.csv(stats_inputpath_1, header = TRUE) %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa_1)

data_2 <- read.csv(stats_inputpath_2, header = TRUE) %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa_2)

# Enriched TFs
enriched_set_1 <- data_1 %>%
  dplyr::filter(abs(GSEA.value) > thr) %>%
  dplyr::pull(TF)

enriched_set_2 <- data_2 %>%
  dplyr::filter(abs(GSEA.value) > thr) %>%
  dplyr::pull(TF)

# Function to prepare data for plotting
prepare_volcano_data <- function(data, highlight_set, title) {
  data %>%
    mutate(
      p.value = ifelse(p.value == 0, runif(n(), 0.00001, 0.001), p.value),
      log10adjpval = -log10(p.value),
      fill = "ES < 65",
      fill = ifelse(p.value < alpha & GSEA.value > thr, "Positive ES", fill),
      fill = ifelse(p.value < alpha & GSEA.value < -thr, "Negative ES", fill),
      toshow = ifelse(
        data[[selec_column_name]] %in% highlight_set & 
          (GSEA.value > thr | GSEA.value < -thr),
        data[[selec_column_name]],
        NA
      ),
      title = title
    )
}

filtered_data_1 <- data_1 %>%
  dplyr::filter(TF %in% enriched_set_2)

filtered_data_2 <- data_2 %>%
  dplyr::filter(TF %in% enriched_set_1)

toplot_data_1 <- prepare_volcano_data(
  filtered_data_1, enriched_set_2, "Sanger TFs in PRISM"
)
toplot_data_2 <- prepare_volcano_data(
  filtered_data_2, enriched_set_1, "PRISM TFs in Sanger"
)

# Function to create a single volcanoplot
create_volcanoplot <- function(toplot) {
  ggscatter(
    toplot,
    x = "GSEA.value",
    y = "log10adjpval",
    color = "fill",
    size = 5,
    label = NULL,
    repel = TRUE,
    title = unique(toplot$title)
  ) +
    geom_text_repel(
      data = toplot[!is.na(toplot$toshow), ],
      aes(label = toshow),
      nudge_x = ifelse(
        toplot$GSEA.value[!is.na(toplot$toshow)] > 0,
        abs(toplot$GSEA.value[!is.na(toplot$toshow)]),
        -abs(toplot$GSEA.value[!is.na(toplot$toshow)])
      ) * 0.05,
      direction = "y",
      max.overlaps = Inf,
      force = 1,
      force_pull = 0.15,
      box.padding = 0.15,
      point.padding = 0.05
    ) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-thr, thr), linetype = "dashed", color = "red") +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 10),
      limits = c(-1, 1)
    ) +
    scale_color_brewer(palette = "Set1")
}

# Función para crear y guardar los gráficos
save_volcano_plot <- function() {
  volcano_data_1 <- create_volcanoplot(toplot_data_1)
  volcano_data_2 <- create_volcanoplot(toplot_data_2)
  
  combined_plot <- grid.arrange(
    ggplotGrob(volcano_data_1),
    ggplotGrob(volcano_data_2),
    ncol = 1
  )
  
  dir.create(dirname(outputpath), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = outputpath,
    plot = combined_plot,
    dpi = dpi,
    width = width,
    height = height
  )
}

# Ejecutar la función para guardar los gráficos
invisible(save_volcano_plot())