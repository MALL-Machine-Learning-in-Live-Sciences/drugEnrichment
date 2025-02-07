
# 1. Inputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
h_matrix_path = "data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv"

# 2. Outputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
his_path <- "figures/histogram/hmatrix_sanger.png"
dir.create(dirname(his_path), recursive = TRUE, showWarnings = FALSE)


# 3. Plot 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
hMatrix <- read.csv(h_matrix_path, header=TRUE, row.names=1)
valores <- as.numeric(as.matrix(hMatrix))

png(his_path)
hist(valores, 
     breaks = 20,            
     xlim = c(0, 1),         
     col = "lightblue",     
     border = "black",      
     xlab = "Valor",         
     ylab = "Frecuencia",   
     main = "Histograma de valores de hMatrix" 
)
dev.off() 