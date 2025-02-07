# Differential expression 
require(dplyr)
require(DESeq2)
require(data.table)
require(dplyr)
library(org.Hs.eg.db)


# 1. Inputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------

counts_path <- "extdata/GSE162104/GSE162104_count.txt"
DSE_path = "data/DE_NMF_SANGER/DE_NMF_SANGER.csv"
h_matrix_path = "data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv"
metadata_path <- "data/GSE162104_processed/GSE162104_metadata.csv"


# 2. Outputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------

output_folder <- "data/diffexp/"
de_genes_path <- "data/diffexp/DEA_cohort.csv"
nomr_counts_path <- "data/diffexp/ddseq_norm_counts.csv"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


# 3. DEA
# ---------------------------------------------------------------------------------------------------------------------------------------------------

counts <- read.table(counts_path, header=TRUE, sep="\t", row.names=1)
metadata <- read.csv(metadata_path)
sampleTable <- metadata %>% 
            select(Sample_title, Sample_characteristics_ch1.1) %>%
            rename(sampleName = Sample_title, Group = Sample_characteristics_ch1.1)
sampleTable <- sampleTable[match(colnames(counts), sampleTable$sampleName), ]
sampleTable$Group <- as.factor(sampleTable$Group)
sampleTable$Group <- relevel(sampleTable$Group, "cetuximab resistance: resistant")

ddsHTSeq <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = sampleTable,
                                   design = ~ Group)

dds <- estimateSizeFactors(ddsHTSeq)
idx <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 10
dds <- dds[idx,]

ddsMF <- DESeq(dds)
res <- results(ddsMF,contrast=list("Group_cetuximab.resistance..sensitive_vs_cetuximab.resistance..resistant"))
resOrdered1 <- res[order(res$padj),]
my_genes1 <- subset(resOrdered1, padj<0.05 & abs(log2FoldChange)>1.0)


# 4. Annotate results
# ---------------------------------------------------------------------------------------------------------------------------------------------------
anno_all_1 <- AnnotationDbi::select(org.Hs.eg.db,keys=rownames(res),
              columns=c("SYMBOL","GENENAME"),
              keytype="ENSEMBL")
anno_all_1 <- anno_all_1[!duplicated(anno_all_1$ENSEMBL), ]  
anno_all1_2 <- anno_all_1[,-1]
rownames(anno_all1_2) <- anno_all_1[,1]
d_1 <- cbind(rownames(anno_all1_2), data.frame(anno_all1_2, row.names=NULL))
d_1 <- d_1 %>% mutate(SYMBOL = coalesce(SYMBOL,rownames(anno_all1_2)))
rownames(d_1) <- d_1[,1]
d_1 <- d_1[,-1]
d1_2 <- merge(as.data.frame(res), d_1,by = 'row.names', all = TRUE)
d1_2[,1] <- d1_2[,8]
res_A <- d1_2[,-c(8:9)]


write.csv(res_A, file = de_genes_path)

normalized_counts <- counts(ddsMF, normalized=TRUE)
df <- log(normalized_counts+1,2)
d2_2 <- merge(as.data.frame(normalized_counts), d_1,by = 'row.names', all = TRUE)
d2_2[[ncol(d2_2) - 1]] <- make.unique(as.character(d2_2$SYMBOL))
rownames(d2_2) <- d2_2$SYMBOL
d2_2 <- d2_2[, !(names(d2_2) %in% c("SYMBOL", "GENENAME"))]


write.csv(d2_2, file = nomr_counts_path)
