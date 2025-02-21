library("msigdbr")
require(DOSE)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(reshape2)
library(dplyr)
set.seed(1234)


# Inputs 
in2 <- "data/diffexpB/WTdiffexp/DEA_cohort.csv"
universe2 <- "data/diffexpB/WTdiffexp/ddseq_norm_counts.csv"

# Outputs
enrichment_plots_path <- "figures/adicional/oraCohorteB/"
select_name <- "dea_oraWTneg"
dir.create(enrichment_plots_path, showWarnings = FALSE)

# Plotting 

stats2 <- read.csv(in2, header = TRUE)
uni2 <- read.csv(universe2,row.names = 1, header = TRUE)
stats2f <- stats2 %>%
  dplyr::filter(padj < 0.05, (log2FoldChange)<(0))


GoI1 <- stats2f$X
Universe1 <- rownames(uni2)


# # Obtain hallmark gene sets for Homo sapiens
#h_gene_sets = msigdbr(species = "human", category = "H")
# h_gene_sets_KEGG = msigdbr(species = "human", category = "C2",subcategory = "CP:KEGG")
# h_gene_sets_BP = msigdbr(species = "human", category = "C5",subcategory = "GO:BP")
h_gene_sets_MF = msigdbr(species = "human", category = "C5",subcategory = "GO:CC")
#h_gene_sets_CC = msigdbr(species = "human", category = "C5",subcategory = "GO:CC")

#msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# msigdbr_t2g_KEGG = h_gene_sets_KEGG %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
# msigdbr_t2g_BP = h_gene_sets_BP %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g = h_gene_sets_MF %>% dplyr::distinct(gs_name, ensembl_gene) %>% as.data.frame()
#msigdbr_t2g = h_gene_sets_CC %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

# # # Make ORA anaylisis
ORA1 <- enricher(gene = GoI1, pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize=5,
universe = Universe1,TERM2GENE = msigdbr_t2g)

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
path <- paste0(enrichment_plots_path, select_name, "_ORACC.png") 
ggsave(path, plot = p, width = 7, height = 5, dpi = 1500)
}
