import decoupler as dc
import anndata as ad
import pandas as pd
import os 

# Inputs 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

h_matrix_path = 'data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv'
with open(h_matrix_path) as file:
    hmatrix = pd.read_csv(file, index_col=0)

# Outputs 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

matrix_name = 'TF_FMN_40_Sanger.csv'
folder_path = 'data/DE_NMF_SANGER/TFs_matrix/'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# TF activity inference 
# ----------------------------------------------------------------------------------------------------------------------------------------------------

hmatrix.columns = [col.split(" ")[0] for col in hmatrix.columns]
net = dc.get_collectri(organism='human')
adata = ad.AnnData(X=hmatrix.values, obs=pd.DataFrame(index=hmatrix.index), var=pd.DataFrame(index=hmatrix.columns))
dc.run_ulm(mat=adata,net=net,source='source',target='target',weight='weight',verbose=True,use_raw = False,batch_size=10000, min_n=5)
tfmatrix = adata.obsm['ulm_estimate']
tfmatrix.to_csv((folder_path+matrix_name))

