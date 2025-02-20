library("msigdbr")
require(DOSE)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(reshape2)
set.seed(1234)

# nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
# rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
# offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
# offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
# child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
# parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")


# Inputs 
DEprismPath <- "data/DE_TF_PRISM/DE_TF_PRISM.csv"
DEsangerPath <- "data/DE_TF_SANGER/DE_TF_SANGER.csv"
netPath <- "data/collectri.csv"
universePath <- "data/prism_processed/counts_matched_prism.csv"

# Outputs
enrichment_plots_path <- "figures/ZenrichmentPlot/"
select_name <- "TF"
dir.create(enrichment_plots_path, showWarnings = FALSE)

# Plotting 
tfPrims <- read.csv(DEprismPath)
tfSanger <- read.csv(DEsangerPath)
net <- read.csv(netPath)
counts <- read.csv(universePath, check.names = FALSE)

stats2 <- read.csv(in2, header = TRUE)
stats2 <- stats2 %>%
  dplyr::filter(pvalue < 0.05)

# Filter the results
DE1positive <- tfPrims %>%
  dplyr::filter(Tumor.type == 'Colorectal_Adenocarcinoma', MoA == 'EGFR inhibitor', (p.value) < 0.05, (GSEA.value)>0.65) %>%
  dplyr::pull(TF)

DE2positive <- tfSanger %>%
  dplyr::filter(Tumor.type == 'Colorectal_Adenocarcinoma', MoA == 'EGFR signaling',(p.value) < 0.05, (GSEA.value)>0.65) %>%
  dplyr::pull(TF)

TF_intersection <- intersect(DE1positive, DE2positive)
GoI1 <- net %>%
  dplyr::filter(source %in% TF_intersection) %>%
  dplyr::pull(target)


  

Universe1 <- colnames(counts)
# #h_gene_sets = msigdbr(species = "human", category = "H")
h_gene_sets= msigdbr(species = "human", category = "C5",subcategory = "GO:MF")

msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()


# ORA1 <- enricher(gene = GoI1, pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize=5,
# universe = Universe1,TERM2GENE = msigdbr_t2g)
# par(mfrow=c(2,1)) 
# if(dim(ORA1)[1] < 1){
#   print("No over-represented terms were found")
#   } else {
# ORA1@result$Description <- gsub("GOMF", "", ORA1@result$Description)
# p_bar <- barplot(ORA1,font.size = 7,showCategory =10)
# p_bar <- p_bar + 
#   coord_flip() +  
#   theme(
#     text = element_text(family = "Times New Roman"),  
#     axis.text = element_text(size = 10, angle = 60, hjust = 1),  # Ajustar el tamaño de las etiquetas y orientación
#     plot.title = element_text(size = 14, face = "bold"),  # Ajustar el tamaño del título
#     plot.margin = margin(10, 10, 10, 10)  # Añadir márgenes para evitar superposición
#   ) +
#   scale_fill_gradient(
#     low = "#36274c",  
#     high = "#5b4e6f"  
#   )
# path <- paste0(enrichment_plots_path, select_name, "_tfcc_plot.png")
# ggsave(path, plot = p_bar, width = 5.25, height = 3.5, dpi = 1500)}





