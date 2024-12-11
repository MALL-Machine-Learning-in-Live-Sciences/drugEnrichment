# 1. Preprocessing Sanger
# ---------------------------------------------------------------------------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(ggtext)
library(dplyr)
library(tidyr)
library(tibble)
library(PupillometryR)
library(ggeasy)
library(biomaRt)



# 1.2. Inputs
# ---------------------------------------------------------------------------------------------------------------------------------------------------

metadata_path <- 'extdata/GSE162104/GSE162104_series_matrix.txt'
counts_path <-  'extdata/GSE162104/GSE162104_count.txt'

# 1.2. Outputs
# ---------------------------------------------------------------------------------------------------------------------------------------------------

cleaned_metadata_path <- 'data/GSE162104_processed/GSE162104_metadata.csv'
cleaned_counts_path <-  'data/GSE162104_processed/GSE162104_counts.csv'
post_normalization_plot <- 'figures/boxplots/postnorm.png'

metadata_dir <- dirname(cleaned_metadata_path)
counts_dir <- dirname(cleaned_counts_path)
pre_plot_dir <- dirname(pre_normalization_plot)
post_plot_dir <- dirname(post_normalization_plot)
dirs <- c(metadata_dir, counts_dir, pre_plot_dir, post_plot_dir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

# 1.3. Preprocessing metadata
# ---------------------------------------------------------------------------------------------------------------------------------------------------
lines <- readLines(metadata_path)

headers <- c()
data <- list()

for (line in lines) {
  if (startsWith(line, "!")) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) > 2) {
      headers <- c(headers, gsub("!", "", parts[1]))  
      data <- append(data, list(parts[-1]))          
    }
  }
}

metadata <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)
colnames(metadata) <- headers

metadata[] <- lapply(metadata, function(x) gsub('"', '', x))

write.csv(metadata, cleaned_metadata_path, row.names = FALSE)

# 1.4. Preprocessing counts  
# ---------------------------------------------------------------------------------------------------------------------------------------------------

counts_matrix <- as.matrix(read.table(counts_path, header = TRUE, row.names = 1))

dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = DataFrame(sample=colnames(counts_matrix)), design = ~ 1)
dds <- estimateSizeFactors(dds) 

idx <- rowSums(counts(dds, normalized = TRUE) >= 12) >= 12
dds <- dds[idx, ]

normalized_counts <- counts(dds, normalized = TRUE)
df <- log(normalized_counts+1,2)

write.csv(t(df), cleaned_counts_path, row.names = TRUE)

# 1.4. Normalization plots
# ---------------------------------------------------------------------------------------------------------------------------------------------------


png(filename = post_normalization_plot, width = 2400, height = 1200, res = 300) 
boxplot(df,
        ylab = "log2 normalized expression",
        col = rep("chartreuse3", ncol(df)),  
        pars = list(xaxt = "n"))  

axis(1, at = c(1:ncol(df)), labels = FALSE)  
text(c(1:ncol(df)), par("usr")[3] - 1, cex = 0.7, labels = colnames(df), 
     srt = 90, xpd = TRUE, offset = 0, adj = 1)  

means <- colMeans(df)
points(means, col = "black", pch = 19)

dev.off()







