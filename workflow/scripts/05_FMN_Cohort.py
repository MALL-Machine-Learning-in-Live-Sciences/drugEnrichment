import pandas as pd 
import numpy as np
import os 
import decoupler as dc
import anndata as ad

# 1. Inputs
# ----------------------------------------------------------------------------------------------------------------------------------------------------

h_matrix_path = 'data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv'
counts_path = 'data/GSE162104_processed/GSE162104_counts.csv'

with open(h_matrix_path) as file:
    h_matrix = pd.read_csv(file, index_col=0)
    
with open(counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
    
    
# 2. Outputs
# ----------------------------------------------------------------------------------------------------------------------------------------------------
folder_path = 'data/nmf_cohort/'
outputs_path = 'data/nmf_cohort/nmf_cohort_wmean.csv'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    
    
# 3. Prepare the net
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


    
# 4. Weighted mean 
# ----------------------------------------------------------------------------------------------------------------------------------------------------
net['weight'] = pd.to_numeric(net['weight'], errors='coerce')
mat = ad.AnnData(X=df_counts.values, obs=pd.DataFrame(index=df_counts.index), var=pd.DataFrame(index=df_counts.columns))
dc.run_wmean(mat, net, source='source', target='target', weight='weight', min_n=5, seed=1,verbose=False, use_raw=False)
fmn_activity  = mat.obsm['wmean_estimate']
fmn_activity.to_csv(outputs_path, index=True)