import decoupler as dc 
import pandas as pd
import os 

# Inputs 
fmn_path = 'data/DE_NMF_SANGER/TFs_matrix/TF_FMN_40_Sanger.csv'
with open(fmn_path) as file:
    fmn_matrix = pd.read_csv(file, index_col=0)

# Outputs 
barplot_name = 'TOP25_SANGER_FMN40_11.png'
folder_path = 'figures/barplot/SANGER/' 

if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    
# Select the FMN Factor
selected_factor = 'Factor11'
fmn_matrix = fmn_matrix[fmn_matrix.index.isin([selected_factor])]

# Plot 
dc.plot_barplot(
    acts=fmn_matrix,
    contrast=selected_factor,
    top=25,
    vertical=True,
    figsize=(6, 12),
    save=((folder_path+barplot_name))
)