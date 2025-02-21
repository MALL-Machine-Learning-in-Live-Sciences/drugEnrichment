from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import decoupler as dc
import os

# Inputs 
counts_path = 'data/prism_processed/counts_matched_prism.csv'
DE_results = 'data/tf_de_PRISM/gsea_results_tf.csv'

# Outputs
figure_name = 'jaccard_065_prism.png'
output_dir = 'figures/jaccardplots/'
os.makedirs(output_dir, exist_ok=True)

# Selecting variables 
thr = 0.65
tt = 'Colorectal_Adenocarcinoma'
moa = 'EGFR inhibitor'

# Open dfs
with open(counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(DE_results) as file:
    df_de = pd.read_csv(file)

# Filter DE results
df_de = df_de[df_de['Tumor type']==tt]
df_de = df_de[df_de['MoA']==moa]
gsea_results = df_de['TF'][np.abs(df_de['GSEA value'])>thr].tolist()

# Select the target Genes 
tf_net = dc.get_collectri(organism='human', split_complexes=False)
filtered_net =  tf_net[tf_net['source'].isin(gsea_results)]
filtered_net = filtered_net[filtered_net['target'].isin(df_counts.columns)]
grouped = filtered_net.groupby('source')['target'].apply(set).to_dict()

# Compute modified Jaccard index 
sources = list(grouped.keys())
n = len(sources)
jaccard_matrix = pd.DataFrame(0.0, index=sources, columns=sources)
for src1, src2 in combinations(sources, 2):
    set1 = grouped[src1]
    set2 = grouped[src2]
    jaccard_1 = len(set1 & set2) / len(set1)
    jaccard_2 = len(set1 & set2) / len(set2)
    jaccard_matrix.loc[src1, src2] = jaccard_1
    jaccard_matrix.loc[src2, src1] = jaccard_2 
for src in sources:
    jaccard_matrix.loc[src, src] = 1.0

# Ploting
figure_path = output_dir+figure_name
plt.figure(figsize=(8, 6))
sns.heatmap(jaccard_matrix, annot=False, cmap="coolwarm", fmt=".2f", cbar=True)
plt.title("Modified Jaccard, ES>0.65", fontsize=16)
plt.xlabel("TF", fontsize=12)
plt.ylabel("TF", fontsize=12)
plt.savefig(figure_path, dpi=300, bbox_inches="tight")
plt.close()