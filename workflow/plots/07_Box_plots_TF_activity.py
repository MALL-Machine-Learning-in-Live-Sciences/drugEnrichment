import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import os 

# Inputs
tf_activity_matrix = pd.read_csv('data/GSE162104_processed/TFs_GSE162104.csv', index_col=0)
metadata =  pd.read_csv('data/GSE162104_processed/GSE162104_metadata.csv', index_col=0)
DE_sanger = pd.read_csv('data/DE_TF_SANGER/DE_TF_SANGER.csv', index_col=0)
DE_prism = pd.read_csv('data/DE_TF_PRISM/DE_TF_PRISM.csv', index_col=0)

#Outputs 
boxplot_folder = 'figures/boxplots_TFa'
os.makedirs(boxplot_folder, exist_ok=True)
box_path = os.path.join(boxplot_folder, 'boxplot_exluded_BRAF_intersect.png')

#Filtering
#metadata = metadata[metadata['Sample_characteristics_ch1.2']!='mutation (wt, braf, kras): KRAS']
metadata = metadata[metadata['Sample_characteristics_ch1.2']!='mutation (wt, braf, kras): BRAF']

DE_sanger = DE_sanger[DE_sanger['Tumor type']=='Colorectal_Adenocarcinoma']
DE_sanger = DE_sanger[DE_sanger['MoA']=='EGFR signaling']
DE_sanger = DE_sanger[np.abs(DE_sanger['GSEA value'])>0.65]

DE_prism = DE_prism[DE_prism['Tumor type']=='Colorectal_Adenocarcinoma']
DE_prism = DE_prism[DE_prism['MoA']=='EGFR inhibitor']
DE_prism = DE_prism[np.abs(DE_prism['GSEA value'])>0.65]

TF_prism = set(DE_prism['TF']).intersection(tf_activity_matrix.columns.tolist())
TF_sanger = set(DE_sanger['TF']).intersection(tf_activity_matrix.columns.tolist())
TF_instersect = TF_prism.intersection(TF_sanger)

resistant = metadata['Sample_title'][metadata['Sample_characteristics_ch1.1']=='cetuximab resistance: resistant']
sensitive = metadata['Sample_title'][metadata['Sample_characteristics_ch1.1']=='cetuximab resistance: sensitive']

tf_activity_matrix = tf_activity_matrix[TF_instersect]
# Plot
tf_activity_matrix["Group"] = tf_activity_matrix.index.map(
    lambda x: "Resistant" if x in resistant.values else "Sensitive" if x in sensitive.values else None
)

tf_activity_matrix_filtered = tf_activity_matrix.dropna(subset=["Group"])
tf_activity_melted = tf_activity_matrix_filtered.melt(
    id_vars="Group", var_name="TF", value_name="Activity"
)

plt.figure(figsize=(15, 8))  # Ajustar el tamaño de la figura
sns.boxplot(data=tf_activity_melted, x="TF", y="Activity", hue="Group")

plt.xlabel("TF")
plt.ylabel("Activity")
plt.legend(title="Group")
plt.xticks(rotation=45)  
plt.tight_layout()
plt.savefig(box_path)