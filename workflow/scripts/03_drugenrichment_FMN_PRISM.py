import decoupler as dc
import anndata as ad
from scipy import stats
import pandas as pd
import os 
from sklearn.decomposition import NMF


# 3. Gsea table
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Inputs path:
cleaned_counts_path = 'data/prism_processed/counts_matched_prism.csv'
cleaned_cell_line_information_path = 'data/prism_processed/metadata_ccle_prism.csv'
drug_net_path = 'data/prism_processed/drug_net_prism.csv'
cleaned_drug_response_path = 'data/prism_processed/drug_response_prism.csv'

# Outputs path
folder_path = 'data/nmf_gsea/'
gsea_results_path = f"{folder_path}gsea_results_nmf_100f.csv"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

#Cormatrix folder path
cor_matrix_folder_path = folder_path+"cor_matrix/"
if not os.path.exists(cor_matrix_folder_path):
    os.makedirs(cor_matrix_folder_path)

nmf_matrix_folder_path = folder_path+"nmf_matrix/"
if not os.path.exists(nmf_matrix_folder_path):
    os.makedirs(nmf_matrix_folder_path)
    
# 3.1 Selecting tumor types + factors
# ----------------------------------------------------------------------------------------------------------------------------------------------------

#  Avalible tumor types 

# 'Colorectal Adenocarcinoma', 'Melanoma',
# 'Bladder Urothelial Carcinoma', 'Non-Small Cell Lung Cancer',
# 'Ovarian Epithelial Tumor', 'Invasive Breast Carcinoma',
# 'Pancreatic Adenocarcinoma', 'Diffuse Glioma', 'Sarcoma, NOS',
# 'Renal Cell Carcinoma', 'Ewing Sarcoma', 'Fibrosarcoma',
# 'Osteosarcoma', 'Pleural Mesothelioma', 'Prostate Adenocarcinoma','Rhabdoid Cancer', 'Neuroblastoma',
# 'Adenosquamous Carcinoma of the Pancreas', 'Rhabdomyosarcoma',
# 'Intracholecystic Papillary Neoplasm','Head and Neck Squamous Cell Carcinoma', 'Embryonal Tumor',
# 'Anaplastic Thyroid Cancer', 'Ampullary Carcinoma',
# 'Endometrial Carcinoma','Intraductal Papillary Neoplasm of the Bile Duct',
# 'Hepatocellular Carcinoma', 'Esophagogastric Adenocarcinoma',
# 'Lung Neuroendocrine Tumor', 'Esophageal Squamous Cell Carcinoma',
# 'Pancreatic Neuroendocrine Tumor', 'Chondrosarcoma',
# 'Uterine Sarcoma/Mesenchymal','Poorly Differentiated Thyroid Cancer', 'Leiomyosarcoma',
# 'Breast Ductal Carcinoma In Situ', 'Hepatoblastoma','Well-Differentiated Thyroid Cancer',
# 'Undifferentiated Pleomorphic Sarcoma/Malignant Fibrous Histiocytoma/High-Grade Spindle Cell Sarcoma',
# 'Bladder Squamous Cell Carcinoma', 'Urethral Cancer','Medullary Thyroid Cancer'

#  Avalible organ

# 'Bowel', 'Skin', 'Bladder/Urinary Tract', 'Lung',
# 'Ovary/Fallopian Tube', 'Breast', 'Pancreas', 'CNS/Brain',
# 'Soft Tissue', 'Kidney', 'Bone', 'Pleura', 'Prostate',
# 'Peripheral Nervous System', 'Biliary Tract', 'Head and Neck',
# 'Thyroid', 'Ampulla of Vater', 'Uterus', 'Liver',
# 'Esophagus/Stomach'

tumor_type_list = ['Colorectal Adenocarcinoma']
factors = [100] #as list, its posible to try any number
iterations = 2000
    
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

def run_NMF(adataX, n_components=2, iters=200, random_state=0):
    """Run NMF on adata

    Run NMF on adata given number of components and random state

    Args:
        adataX (np.array): num of cells x num of genes
        iters (int): max number of iters for NMF algorithm
        n_components (int): number of components for NMF
        random_state (int): random seed
    
    Returns:
        Metagenes and factors

    """
    model = NMF(n_components, init='random', max_iter=iters, random_state=random_state)
    W = model.fit_transform(adataX)
    H = model.components_

    # Format objects
    W = pd.DataFrame(W, index=adataX.index, columns=[f"Factor{i}" for i in range(1, n_components + 1)])
    H = pd.DataFrame(H, index=[f"Factor{i}" for i in range(1, n_components + 1)], columns=adataX.columns) 

    return W, H

with open(cleaned_counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(cleaned_cell_line_information_path) as file:
    df_metadata = pd.read_csv(file, index_col=0)
with open(drug_net_path) as file:
    drug_net = pd.read_csv(file, index_col=0)
with open(cleaned_drug_response_path) as file:
    df_outputs = pd.read_csv(file, index_col=0)

# 2.4 GSEA splicing by organ:
# ----------------------------------------------------------------------------------------------------------------------------------------------------

df_resultado = pd.DataFrame()


for organ in df_metadata['OncotreeLineage'].unique():
    if organ in tumor_type_list:
        cell_lines = df_metadata.index[df_metadata['OncotreeLineage']==organ].to_list()
        df_sub_counts = df_counts[df_counts.index.isin(cell_lines)]
        for i in factors:
            w , h = run_NMF(df_sub_counts, n_components=i, iters=iterations)
            df_sub_outputs = df_outputs[df_outputs.index.isin(cell_lines)]
            try:
                if len(df_sub_outputs)>10 and len(df_sub_counts)>10:
                    cor_matrix = correlation_drugs(df_counts=w,df_drugs=df_sub_outputs)
                    gsea_result, pvalues = gsa_drug_enrichment(matrix=cor_matrix,net=drug_net, shared_element='broad_id',net_group='moa')

                    organ = organ.replace("/","_").replace(" ","_")
                    final_path_cormatrix = cor_matrix_folder_path+organ+f"_{i}_cormatrix.csv"
                    cor_matrix.to_csv(final_path_cormatrix)
                    final_nmf_matrix_path_w = nmf_matrix_folder_path+organ+f"_{i}_w_matrix.csv"
                    final_nmf_matrix_path_h = nmf_matrix_folder_path+organ+f"_{i}_h_matrix.csv"
                    w.to_csv(final_nmf_matrix_path_w)
                    h.to_csv(final_nmf_matrix_path_h)
                    
                    df_gsea_long = gsea_result.reset_index().melt(id_vars='index', var_name='MoA', value_name='GSEA value')
                    df_gsea_long.rename(columns={'index': 'Factors'}, inplace=True)
                    df_pvalues_long = pvalues.reset_index().melt(id_vars='index', var_name='MoA', value_name='p-value')
                    df_pvalues_long.rename(columns={'index': 'Factors'}, inplace=True)
                    df_resultado_org = pd.merge(df_gsea_long, df_pvalues_long, on=['Factors', 'MoA'])
                    df_resultado_org['NMF_components'] = i
                    df_resultado_org['Cells_count'] = len(w)
                    df_resultado_org['Split by'] = 'Organ'
                    df_resultado_org['Tumor type'] = organ    
                    df_resultado = pd.concat([df_resultado, df_resultado_org], ignore_index=True)
                else:
                    continue
            except ValueError:
                continue

# 2.5 GSEA splicing by tumor type:
# ----------------------------------------------------------------------------------------------------------------------------------------------------

for organ in df_metadata['OncotreePrimaryDisease'].unique():
    if organ in tumor_type_list:
        cell_lines = df_metadata.index[df_metadata['OncotreePrimaryDisease']==organ].to_list()
        df_sub_counts = df_counts[df_counts.index.isin(cell_lines)]
        for i in factors:
            w , h = run_NMF(df_sub_counts, n_components=i, iters=iterations)
            df_sub_outputs = df_outputs[df_outputs.index.isin(cell_lines)]
            try:
                if len(df_sub_outputs)>10 and len(df_sub_counts)>10:
                    cor_matrix = correlation_drugs(df_counts=w,df_drugs=df_sub_outputs)
                    gsea_result, pvalues = gsa_drug_enrichment(matrix=cor_matrix,net=drug_net, shared_element='broad_id',net_group='moa')

                    organ = organ.replace("/","_").replace(" ","_")
                    final_path_cormatrix = cor_matrix_folder_path+organ+f"_{i}_cormatrix.csv"
                    cor_matrix.to_csv(final_path_cormatrix)
                    final_nmf_matrix_path_w = nmf_matrix_folder_path+organ+f"_{i}_w_matrix.csv"
                    final_nmf_matrix_path_h = nmf_matrix_folder_path+organ+f"_{i}_h_matrix.csv"
                    w.to_csv(final_nmf_matrix_path_w)
                    h.to_csv(final_nmf_matrix_path_h)
                    
                    df_gsea_long = gsea_result.reset_index().melt(id_vars='index', var_name='MoA', value_name='GSEA value')
                    df_gsea_long.rename(columns={'index': 'Factors'}, inplace=True)
                    df_pvalues_long = pvalues.reset_index().melt(id_vars='index', var_name='MoA', value_name='p-value')
                    df_pvalues_long.rename(columns={'index': 'Factors'}, inplace=True)
                    df_resultado_org = pd.merge(df_gsea_long, df_pvalues_long, on=['Factors', 'MoA'])
                    df_resultado_org['NMF_components'] = i
                    df_resultado_org['Cells_count'] = len(w)
                    df_resultado_org['Split by'] = 'Organ'
                    df_resultado_org['Tumor type'] = organ    
                    df_resultado = pd.concat([df_resultado, df_resultado_org], ignore_index=True)
                else:
                    continue
            except ValueError:
                continue


df_resultado.to_csv(gsea_results_path)