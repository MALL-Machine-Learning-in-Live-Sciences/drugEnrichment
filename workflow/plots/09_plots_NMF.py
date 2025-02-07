import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import os 
from scipy.stats import kruskal

# Inputs
tf_activity_matrix = pd.read_csv('data/nmf_cohort/nmf_cohort_ulm.csv', index_col=0)
metadata =  pd.read_csv('data/GSE162104_processed/GSE162104_metadata.csv', index_col=0)
DE_sanger = pd.read_csv('data/DE_NMF_SANGER/DE_NMF_SANGER.csv', index_col=0)

#Outputs 
barplot_folder = 'figures/Boxplots_cohort/ulm_prims/'
os.makedirs(barplot_folder, exist_ok=True)

DE_sanger = DE_sanger[DE_sanger['MoA']=='EGFR signaling']
DE_sanger = DE_sanger[DE_sanger['NMF_components']==40]
DE_sanger = DE_sanger[np.abs(DE_sanger['GSEA value'])>0.6]

resistant = metadata['Sample_title'][metadata['Sample_characteristics_ch1.1']=='cetuximab resistance: resistant']
sensitive = metadata['Sample_title'][metadata['Sample_characteristics_ch1.1']=='cetuximab resistance: sensitive']

tf_activity_matrix["Group"] = tf_activity_matrix.index.map(
     lambda x: "Resistant" if x in resistant.values else "Sensitive" if x in sensitive.values else None)
tf_activity_matrix_filtered = tf_activity_matrix.dropna(subset=["Group"])

for col in DE_sanger['Factors'].to_list():
        grouped_data = {}
        for group_name, group_data in tf_activity_matrix_filtered.groupby('Group'):
            grouped_data[group_name] = group_data[col].dropna()

        groups = list(grouped_data.values())
        stat, p_value = kruskal(*groups)

        plt.figure(figsize=(8, 6))
        sns.boxplot(x='Group', y=col, data=tf_activity_matrix_filtered)

        plt.title(f"{col} (Kruskal-Wallis: H={stat:.2f}, p={p_value:.3f})")
        plt.xlabel('Group')
        plt.ylabel(col)
        
        file_path = os.path.join(barplot_folder, f"{col}.png")
        plt.savefig(file_path)
        plt.close()