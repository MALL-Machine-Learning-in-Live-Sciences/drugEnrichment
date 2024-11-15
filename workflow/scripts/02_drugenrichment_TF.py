import decoupler as dc
import anndata as ad
from scipy import stats
import pandas as pd
import os 


# 2. Gsea table
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Inputs path:
cleaned_counts_path = 'data/prism_processed/counts_matched_prism.csv'
cleaned_cell_line_information_path = 'data/prism_processed/metadata_ccle_prism.csv'
drug_net_path = 'data/prism_processed/drug_net_prism.csv'
cleaned_drug_response_path = 'data/prism_processed/drug_response_prism.csv'

# Outputs path
folder_path = 'data/tf_gsea/'
gsea_results_path = f"{folder_path}gsea_results_tf.csv"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

#Cormatrix folder path
cor_matrix_folder_path = folder_path+"cor_matrix/"
if not os.path.exists(cor_matrix_folder_path):
    os.makedirs(cor_matrix_folder_path)

# 2.1 Defining functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------
def compute_correlation(drug_col, stats, df_counts, df_drugs):
    correlations = []
    df_suboutputs = df_drugs[drug_col].dropna()
    for gene in df_counts.columns:
        if len(df_suboutputs.unique())>1:
            df_subgene = df_counts[gene].dropna()
            shared_cells = df_suboutputs.index.intersection(df_subgene.index)       
            correlation, _ = stats.pearsonr(df_suboutputs[shared_cells].sort_index().values,
                                    df_subgene[shared_cells].sort_index().values)
            correlations.append((gene, correlation))
        else:
            correlations.append((gene, None))

    return drug_col, correlations
    

def correlation_drugs(df_drugs, df_counts):  
    '''
    Inputs: drug matrix (cell-response) and expresión matrix (cell-expression)
    Outputs: correlation matrix: gene-drug
    '''
    from joblib import Parallel, delayed
    import pandas as pd
    from scipy import stats  


    correlation_matrix = pd.DataFrame(index=df_drugs.columns, columns=df_counts.columns)

    # Ejecutar las correlaciones en paralelo
    results = Parallel(n_jobs=-1)(delayed(compute_correlation)(drug_col, stats, df_counts, df_drugs) for drug_col in df_drugs.columns)

    # Almacenar los resultados en el DataFrame de correlación
    for drug_col, correlations in results:
        for gene, correlation in correlations:
            correlation_matrix.at[drug_col, gene] = correlation

    # Visualiza el DataFrame de correlación
    return correlation_matrix

# Funcion to do the gesea enrichment
def gsa_drug_enrichment(matrix, net, shared_element, net_group, n_min=5):
    '''
    Inputs: 
        matrix: correlation matrix drug-gene
        met: drugs grouped by moa (in this case)
        shared_element: Commun column between matrix-net
        net_pathway: name of the column with the groups in the net
        n_min: min elements to consider 
        
    Outputs: gsea matrix and p-values
    '''
    
    
    matrix = matrix.dropna(0)
    shared_elements = set(matrix.index).intersection(set(net[shared_element]))

    adata = ad.AnnData(matrix.loc[shared_elements, :].T)

    result  = dc.run_gsea(mat=adata,net=net[net[shared_element].isin(shared_elements)],source=net_group,target=shared_element,min_n=5,use_raw=False)
    result_df_gsea = adata.obsm['gsea_estimate'] 
    p_value_gsea = adata.obsm['gsea_pvals']
    
    
    return result_df_gsea, p_value_gsea



# 2.2 TF activity inference
# ----------------------------------------------------------------------------------------------------------------------------------------------------
with open(cleaned_counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(cleaned_cell_line_information_path) as file:
    df_metadata = pd.read_csv(file, index_col=0)
with open(drug_net_path) as file:
    drug_net = pd.read_csv(file, index_col=0)
with open(cleaned_drug_response_path) as file:
    df_outputs = pd.read_csv(file, index_col=0)


tf_net = dc.get_collectri(organism='human', split_complexes=False)
adata = ad.AnnData(X=df_counts.values, obs=pd.DataFrame(index=df_counts.index), var=pd.DataFrame(index=df_counts.columns))
dc.run_ulm(
    mat=adata,
    net=tf_net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
    use_raw = False,
    batch_size=10000, 
    min_n=5
)
df_tf = adata.obsm['ulm_estimate']



# 2.3 GSEA using all cellines:
# ----------------------------------------------------------------------------------------------------------------------------------------------------

cor_matrix = correlation_drugs(df_counts=df_tf,df_drugs=df_outputs)
gsea_result, pvalues = gsa_drug_enrichment(matrix=cor_matrix,net=drug_net, shared_element='broad_id',net_group='moa')

final_path = cor_matrix_folder_path+'cormatrix.csv'
cor_matrix.to_csv(final_path)

df_gsea_long = gsea_result.reset_index().melt(id_vars='index', var_name='MoA', value_name='GSEA value')
df_gsea_long.rename(columns={'index': 'TF'}, inplace=True)
df_pvalues_long = pvalues.reset_index().melt(id_vars='index', var_name='MoA', value_name='p-value')
df_pvalues_long.rename(columns={'index': 'TF'}, inplace=True)
df_resultado = pd.merge(df_gsea_long, df_pvalues_long, on=['TF', 'MoA'])
df_resultado['Cells_count'] = len(df_tf)
df_resultado['Split by'] = None
df_resultado['Tumor type'] = 'All cell lines'

# 2.3 GSEA splicing by organ:
# ----------------------------------------------------------------------------------------------------------------------------------------------------

for organ in df_metadata['OncotreeLineage'].unique():
    cell_lines = df_metadata.index[df_metadata['OncotreeLineage']==organ].to_list()
    df_sub_tf = df_tf[df_tf.index.isin(cell_lines)]
    df_sub_outputs = df_outputs[df_outputs.index.isin(cell_lines)]
    try:
        if len(df_sub_outputs)>10 and len(df_sub_tf)>10:
            cor_matrix = correlation_drugs(df_counts=df_sub_tf,df_drugs=df_sub_outputs)
            gsea_result, pvalues = gsa_drug_enrichment(matrix=cor_matrix,net=drug_net, shared_element='broad_id',net_group='moa')

            organ = organ.replace("/","_").replace(" ","_")
            final_path = cor_matrix_folder_path+organ+'_cormatrix.csv'
            cor_matrix.to_csv(final_path)
            
            df_gsea_long = gsea_result.reset_index().melt(id_vars='index', var_name='MoA', value_name='GSEA value')
            df_gsea_long.rename(columns={'index': 'TF'}, inplace=True)
            df_pvalues_long = pvalues.reset_index().melt(id_vars='index', var_name='MoA', value_name='p-value')
            df_pvalues_long.rename(columns={'index': 'TF'}, inplace=True)
            df_resultado_org = pd.merge(df_gsea_long, df_pvalues_long, on=['TF', 'MoA'])
            df_resultado_org['Cells_count'] = len(df_sub_tf)
            df_resultado_org['Split by'] = 'Organ'
            df_resultado_org['Tumor type'] = organ    
            df_resultado = pd.concat([df_resultado, df_resultado_org], ignore_index=True)
        else:
            continue
    except ValueError:
        continue

# 2.3 GSEA splicing by tumor type:
# ----------------------------------------------------------------------------------------------------------------------------------------------------


for organ in df_metadata['OncotreePrimaryDisease'].unique():
    cell_lines = df_metadata.index[df_metadata['OncotreePrimaryDisease']==organ].to_list()
    df_sub_tf = df_tf[df_tf.index.isin(cell_lines)]
    df_sub_outputs = df_outputs[df_outputs.index.isin(cell_lines)]
    try:
        if len(df_sub_outputs)>10 and len(df_sub_tf)>10:
            cor_matrix = correlation_drugs(df_counts=df_sub_tf,df_drugs=df_sub_outputs)
            gsea_result, pvalues = gsa_drug_enrichment(matrix=cor_matrix,net=drug_net, shared_element='broad_id',net_group='moa')

            organ = organ.replace("/","_").replace(" ","_")
            final_path = cor_matrix_folder_path+organ+'_cormatrix.csv'
            cor_matrix.to_csv(final_path)
            
            df_gsea_long = gsea_result.reset_index().melt(id_vars='index', var_name='MoA', value_name='GSEA value')
            df_gsea_long.rename(columns={'index': 'TF'}, inplace=True)
            df_pvalues_long = pvalues.reset_index().melt(id_vars='index', var_name='MoA', value_name='p-value')
            df_pvalues_long.rename(columns={'index': 'TF'}, inplace=True)
            df_resultado_org = pd.merge(df_gsea_long, df_pvalues_long, on=['TF', 'MoA'])
            df_resultado_org['Cells_count'] = len(df_sub_tf)
            df_resultado_org['Split by'] = 'Organ'
            df_resultado_org['Tumor type'] = organ    
            df_resultado = pd.concat([df_resultado, df_resultado_org], ignore_index=True)
        else:
            continue
    except ValueError:
        continue


df_resultado.to_csv(gsea_results_path)