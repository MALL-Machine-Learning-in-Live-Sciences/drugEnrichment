import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import decoupler as dc 
import anndata as ad 
import os 

# Imputs 
#DE_path ='data/tf_de_PRISM/gsea_results_tf.csv'
counts_path = 'data/prism_processed/counts_matched_prism.csv'
metadata_path = 'data/prism_processed/metadata_ccle_prism.csv'
outputs_path = 'data/prism_processed/drug_response_prism.csv'
mutations_path = 'extdata/DepMap_23Q2/OmicsSomaticMutations.csv'
net_path = 'data/prism_processed/drug_net_prism.csv'

# Outputs 
folder_path = 'figures/regplots/'

# Set variables
thr = 0.65
tumor_type = "Colorectal_Adenocarcinoma"
moa = "EGFR inhibitor"

# Open files 
# with open(DE_path) as file:
#     df_de = pd.read_csv(file, index_col=0)
with open(counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(metadata_path) as file:
    df_metadata = pd.read_csv(file)
with open(outputs_path) as file:
    df_outputs = pd.read_csv(file, index_col=0)
with open(mutations_path) as file: 
    df_mutations = pd.read_csv(file)
with open(net_path) as file: 
    df_net = pd.read_csv(file)
    
# Filter TFs
# df_de = df_de[df_de['MoA'] == moa]
# df_de = df_de[df_de['Tumor type'] == tumor_type]
# df_de = df_de[np.abs(df_de['GSEA value'])>thr]
# tfs_list = df_de['TF'].to_list()
tfs_list = ['HDAC1', 'GLI2']

# Filter net
df_net = df_net[df_net['moa']==moa]

# Filter mutations 
df_metadata.loc[:, 'OncotreePrimaryDisease'] = df_metadata['OncotreePrimaryDisease'].str.replace(' ', '_')
df_metadata = df_metadata[df_metadata['OncotreePrimaryDisease'] == tumor_type]
df_mutations = df_mutations[df_mutations['ModelID'].isin(df_metadata['ModelID'].to_list())]
df_mutations = df_mutations[df_mutations['HugoSymbol'].isin(['KRAS'])]

# Filter outputs 
df_outputs = df_outputs[df_outputs.columns.intersection(df_net['broad_id'])]
df_outputs = df_outputs[df_outputs.index.isin(df_metadata['ModelID'].to_list())]
broads = df_outputs.columns.tolist()

# Filter counts + TFs activity inference 
df_counts = df_counts[df_counts.index.isin(df_metadata['ModelID'].to_list())]
tf_net = dc.get_collectri(organism='human', split_complexes=False)
adata = ad.AnnData(X=df_counts.values, obs=pd.DataFrame(index=df_counts.index), var=pd.DataFrame(index=df_counts.columns))
dc.run_ulm(mat=adata,net=tf_net,source='source',target='target',weight='weight',verbose=True,use_raw = False,batch_size=10000, min_n=5)
df_tf = adata.obsm['ulm_estimate']
df_tf = df_tf[df_tf.columns.intersection(tfs_list)]
plt.rcParams['font.family'] = 'Times New Roman'

# Regplots 
os.makedirs(folder_path, exist_ok=True)

merged = pd.merge(df_counts,df_outputs, left_index=True, right_index=True)
merged['KRAS Mutations'] = merged.index.isin(df_mutations['ModelID'])
for broad in broads:


    num_tfs = len(tfs_list)  
    cols = 3  
    rows = (num_tfs + cols - 1) // cols  

    fig, axes = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    axes = axes.flatten()  

    # Mapeo de la leyenda
    legend_labels = {True: "Presente", False: "Ausente"}
    drug_name = df_net['name'][df_net['broad_id'] == broad].values[0]
    
    for i, tf in enumerate(tfs_list):
        sns.scatterplot(
            x=broad,
            y=tf,
            hue='KRAS Mutations',
            data=merged,
            ax=axes[i],
            palette={True: "#47a4c3d0", False: "#792020dd"},
            s=50
        )
        sns.regplot(
            x=broad,
            y=tf,
            data=merged,
            scatter=False,
            line_kws={'color': '#302543'},
            ax=axes[i]
        )
        axes[i].set_xlabel(f"{drug_name.capitalize()} (AUC)", fontsize=8)
        axes[i].set_ylabel(tf, fontsize=8)
        
        # Cambiar las etiquetas de la leyenda
        handles, labels = axes[i].get_legend_handles_labels()
        new_labels = [legend_labels[eval(label)] for label in labels]
        axes[i].legend(handles, new_labels, title="Mutación KRAS", loc='upper left', fontsize=8)

    # Extraer el nombre del fármaco
    fig.suptitle(f"Regression Plots for {drug_name}", fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    # Guardar la figura
    output_file = os.path.join(folder_path, f"{broad}_regression_plots.png")
    plt.savefig(output_file, dpi=1000, bbox_inches='tight')
    plt.close()