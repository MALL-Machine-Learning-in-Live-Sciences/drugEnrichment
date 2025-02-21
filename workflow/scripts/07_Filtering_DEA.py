import pandas as pd
import numpy as np 
import decoupler as dc
import os 


# 1. Inputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------

de_genes_path = "data/diffexp/DEA_cohort.csv"
DSE_path = "data/DE_NMF_SANGER/DE_NMF_SANGER.csv"
h_matrix_path = "data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv"
DETF_sanger_path = "data/DE_TF_SANGER/DE_TF_SANGER.csv"
DETF_prism_path = "data/DE_TF_PRISM/DE_TF_PRISM.csv"


# 2. Outputs 
# ---------------------------------------------------------------------------------------------------------------------------------------------------

folder_path = "data/diffexp/"
fmnFiltered_path = "DE_fmntop50_filter.csv"
tfFiltered_path = "DE_tf_filter.csv"
intersected_path = "DE_intersect_filter.csv"

# 3. Filter
# ---------------------------------------------------------------------------------------------------------------------------------------------------

with open(de_genes_path) as file:
    de_genes = pd.read_csv(file, index_col=0)  
with open(DSE_path) as file:
    results_fmn = pd.read_csv(file, index_col=0)  
with open(h_matrix_path) as file:
    h_matrix = pd.read_csv(file, index_col=0)   
with open(DETF_sanger_path) as file:
    results_tf = pd.read_csv(file, index_col=0)
with open(DETF_prism_path) as file:
    results_prism_tf = pd.read_csv(file, index_col=0)

if not os.path.exists(folder_path):
    os.makedirs(folder_path)
        
de_genes = de_genes[de_genes['padj']<0.05]

results_fmn = results_fmn[results_fmn['MoA']=='EGFR signaling']
results_fmn = results_fmn[results_fmn['NMF_components']==40]
results_fmn = results_fmn[np.abs(results_fmn['GSEA value'])>0.70]

results_tf = results_tf[results_tf['MoA']=='EGFR signaling']
results_tf = results_tf[np.abs(results_tf['GSEA value'])>0.70]
results_prism_tf = results_prism_tf[results_prism_tf['MoA']=='EGFR signaling']
results_prism_tf = results_prism_tf[np.abs(results_prism_tf['GSEA value'])>0.60]

net = dc.get_collectri(organism='human')

genes = de_genes['Row.names'].unique().tolist()
tf = results_tf['TF'].to_list()
tfp = results_tf['TF'].to_list()


net = net[net['source'].isin((set(tf).intersection(set(tfp))))]
tf_genes = net['target'].unique().tolist()
set1 = set(genes).intersection(set(tf_genes))

h_matrix = h_matrix.reset_index()
factors = results_fmn['Factors'].to_list()
h_matrix.columns = [col.split(" ")[0] for col in h_matrix.columns]
long_h = h_matrix.melt(id_vars='index', var_name='gene', value_name='value')
long_h = long_h[long_h['index'].isin(factors)]
long_h = long_h.groupby("index").apply(lambda x: x.nlargest(50, "value")).reset_index(drop=True)
set2 = set(genes).intersection(set(long_h['gene']))

filter1 = de_genes[de_genes['Row.names'].isin(set1)]
filter2 = de_genes[de_genes['Row.names'].isin(set2)]
filter3 = de_genes[de_genes['Row.names'].isin(set1.intersection(set2))]

filter1.to_csv((folder_path+fmnFiltered_path))
filter2.to_csv((folder_path+tfFiltered_path))
filter3.to_csv((folder_path+intersected_path))

