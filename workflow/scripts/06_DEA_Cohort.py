import os 
import pandas as pd
import scanpy as sc 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Inputs
# ----------------------------------------------------------------------------------------------------------------------------------------------------
counts_path = 'data/GSE162104_processed/GSE162104_counts.csv'
metadata_path =  'data/GSE162104_processed/GSE162104_metadata.csv'
DE_result_path ='data/DE_NMF_SANGER/DE_NMF_SANGER.csv'
h_matrix_path = 'data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv'

with open(counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(metadata_path) as file:
    metadata = pd.read_csv(file, index_col=0)
with open(DE_result_path) as file:
    DE_sanger = pd.read_csv(file, index_col=0)
with open(h_matrix_path) as file:
    h_matrix = pd.read_csv(file, index_col=0)

# 2. Outputs
# ----------------------------------------------------------------------------------------------------------------------------------------------------
results_path = 'DEA_GSE162104.csv'
filtered_result_path = 'DEA_GSE162104_filter.csv'
folder_path = 'data/GSE162104_processed/'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)
 
    
# 3. Select the Factors 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

DE_sanger = DE_sanger[DE_sanger['MoA']=='EGFR signaling']
DE_sanger = DE_sanger[DE_sanger['NMF_components']==40]
DE_sanger = DE_sanger[np.abs(DE_sanger['GSEA value'])>0.70]


# 3. Filter the counts and the metadata 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

metadata = metadata[['Sample_title', 'Sample_characteristics_ch1.1']]
metadata = metadata[metadata['Sample_characteristics_ch1.1'].isin(['cetuximab resistance: resistant','cetuximab resistance: sensitive'])]
metadata = metadata.set_index('Sample_title')
df_counts = df_counts[df_counts.index.isin(metadata.index)]


# 4. DEA
# ----------------------------------------------------------------------------------------------------------------------------------------------------

adata = sc.AnnData(
    X=df_counts.values,  
    obs=metadata, 
    var=pd.DataFrame(index=df_counts.columns)  
)

sc.tl.rank_genes_groups(adata, groupby='Sample_characteristics_ch1.1', method='wilcoxon', reference='rest', use_raw=False)


# 5. Prepare the net
# ----------------------------------------------------------------------------------------------------------------------------------------------------

h_matrix.columns = [col.split(" ")[0] for col in h_matrix.columns]
h_matrix = h_matrix.rename_axis('source', axis=0) 
h_matrix = h_matrix.reset_index() 
net = pd.melt(
    h_matrix,
    id_vars='source',          
    var_name='target',     
    value_name='weight'  
)

net = net[net['target'].isin(df_counts.columns)]
top = 50

filtered_net = (
    net.sort_values(by='weight', ascending=False)  
       .groupby('source')                          
       .head(top)                                  
       .reset_index(drop=True)                    
)
filtered_net = filtered_net[filtered_net['source'].isin(DE_sanger['Factors'].unique())]
top_genes = filtered_net['target'].unique()

rank_genes = adata.uns['rank_genes_groups']
groups = rank_genes['names'].dtype.names 
logfoldchanges = pd.DataFrame({group: rank_genes['logfoldchanges'][group] for group in groups})
logfoldchanges.index = rank_genes['names'][groups[0]]

pvals = pd.DataFrame({group: rank_genes['pvals'][group] for group in groups})
pval_threshold = 1
filtered_lgc = {}
for group in groups:
    group_df = pd.DataFrame({
        'gene': rank_genes['names'][group],
        'logfoldchange': rank_genes['logfoldchanges'][group],
        'pval': rank_genes['logfoldchanges'][group]
    })
    filtered_lgc[group] = group_df[group_df['pval'] > pval_threshold]
    
resistant = filtered_lgc['cetuximab resistance: resistant'][['gene', 'logfoldchange']].rename(columns={'logfoldchange': 'logfc_resistant'})
sensitive = filtered_lgc['cetuximab resistance: sensitive'][['gene', 'logfoldchange']].rename(columns={'logfoldchange': 'logfc_sensitive'})
combined = pd.merge(resistant, sensitive, on='gene', how='outer')

logfoldchanges_subset = combined[combined['gene'].isin(top_genes)].dropna()

logfoldchanges.to_csv((folder_path+results_path))
logfoldchanges_subset.to_csv((folder_path+filtered_result_path))