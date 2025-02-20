library("msigdbr")
require(DOSE)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(reshape2)
library(dplyr)
set.seed(1234)
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")


# Inputs 
DEnmfPath <- "data/DE_NMF_SANGER/DE_NMF_SANGER.csv"
HmatrixPath <- "data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv"
in2 <- "data/diffexp/DEA_cohort.csv"
universe2 <- "data/diffexp/ddseq_norm_counts.csv"

# Outputs
enrichment_plots_path <- "figures/adicional/"
select_name <- "ORADEAneg"
dir.create(enrichment_plots_path, showWarnings = FALSE)

# Plotting 
# data <- read.csv(DEnmfPath)
# Hmatrix <- read.csv(HmatrixPath, row.names = 1, check.names = FALSE)
# stats2 <- read.csv(in2, header = TRUE)
# uni2 <- read.csv(universe2,row.names = 1, header = TRUE)


tumor_type <- "Colorectal_Adenocarcinoma"
moa <- "EGFR signaling"

filterData <- data %>%
  dplyr::filter(Tumor.type == tumor_type, MoA == moa, NMF_components==40, p.value<0.05)


colnames(Hmatrix) <- sapply(strsplit(colnames(Hmatrix), " "), `[`, 1)


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
  

stats2f <- stats2 %>%
  dplyr::filter(padj < 0.05, (log2FoldChange)>(0)) 

  #  %>%
  # dplyr::filter(Row.names %in% Hmatrix_top50$Gene)


GoI1 <- stats2f$Row.names
Universe1 <- rownames(uni2)


# # Obtain hallmark gene sets for Homo sapiens
#h_gene_sets = msigdbr(species = "human", category = "H")
# h_gene_sets_KEGG = msigdbr(species = "human", category = "C2",subcategory = "CP:KEGG")
# h_gene_sets_BP = msigdbr(species = "human", category = "C5",subcategory = "GO:BP")
h_gene_sets_MF = msigdbr(species = "human", category = "C5",subcategory = "GO:MF")
#h_gene_sets_CC = msigdbr(species = "human", category = "C5",subcategory = "GO:CC")

#msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# msigdbr_t2g_KEGG = h_gene_sets_KEGG %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# msigdbr_t2g_BP = h_gene_sets_BP %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g = h_gene_sets_MF %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
#msigdbr_t2g = h_gene_sets_CC %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

# # # Make ORA anaylisis
ORA1 <- enricher(gene = GoI1, pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize=5,
universe = Universe1,TERM2GENE = msigdbr_t2g)
ORA1@result$Description <- gsub("GOMF", "", ORA1@result$Description)

par(mfrow=c(2,1))
if(dim(ORA1)[1] < 1){
  print("No over-represented terms were found")
  } else {
p <- enrichplot::dotplot(ORA1,font.size = 7,showCategory =50)
p <- p + 
  coord_flip() +  
  theme(
    text = element_text(family = "Times New Roman"),  
    axis.text = element_text(size = 10, angle = 45, hjust = 1),  # Ajustar el tamaño de las etiquetas y orientación
    plot.title = element_text(size = 14, face = "bold"),  # Ajustar el tamaño del título
    plot.margin = margin(10, 10, 10, 10)  # Añadir márgenes para evitar superposición
  ) +
  scale_fill_gradient(
    low = "#36274c",  
    high = "#5b4e6f"  
  )
path <- paste0(enrichment_plots_path, select_name, "_ORA.png") 
ggsave(path, plot = p, width = 12, height = 5, dpi = 1000)
}

# ORA1 <- enricher(gene = GoI1, pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize=5,
# universe = Universe1,TERM2GENE = msigdbr_t2g)
# par(mfrow=c(2,1)) 
# if(dim(ORA1)[1] < 1){
#   print("No over-represented terms were found")
#   } else {
# ORA1@result$Description <- gsub("GOMF", "", ORA1@result$Description)
# p_bar <- barplot(ORA1,font.size = 7,showCategory =10 )
# p_bar <- p_bar + 
#   coord_flip() +  
#   theme(
#     text = element_text(family = "Times New Roman"),  
#     axis.text = element_text(size = 10, angle = 90, hjust = 1),  # Ajustar el tamaño de las etiquetas y orientación
#     plot.title = element_text(size = 14, face = "bold"),  # Ajustar el tamaño del título
#     plot.margin = margin(10, 10, 10, 10)  # Añadir márgenes para evitar superposición
#   ) +
#   scale_fill_gradient(
#     low = "#36274c",  
#     high = "#5b4e6f"  
#   )
# path <- paste0(enrichment_plots_path, select_name, "GOMF.png")
# ggsave(path, plot = p_bar, width = 5, height = 4, dpi = 1000)
# }

