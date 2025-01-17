import decoupler as dc
import anndata as ad
import pandas as pd
import os 

# Inputs 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

counts_path = 'data/GSE162104_processed/GSE162104_counts.csv'
with open(counts_path) as file:
    counts_matrix = pd.read_csv(file, index_col=0)

# Outputs 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

matrix_name = 'TFs_GSE162104.csv'
folder_path = 'data/GSE162104_processed/'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# TF activity inference 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

counts_matrix.columns = [col.split(" ")[0] for col in counts_matrix.columns]
net = dc.get_collectri(organism='human')
adata = ad.AnnData(X=counts_matrix.values, obs=pd.DataFrame(index=counts_matrix.index), var=pd.DataFrame(index=counts_matrix.columns))
dc.run_ulm(mat=adata,net=net,source='source',target='target',weight='weight',verbose=True,use_raw = False,batch_size=10000, min_n=5)
tfmatrix = adata.obsm['ulm_estimate']
tfmatrix.to_csv((folder_path+matrix_name))

