import decoupler as dc
import anndata as ad
from scipy import stats
import pandas as pd
import os 
from sklearn.decomposition import NMF
import Drug_Enrichment as DE


# 3. Gsea table
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Inputs path:
cleaned_counts_path = 'extdata/DepMap_23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv'
cell_line_information_path = 'extdata/DepMap_23Q2/Model.csv'
cleaned_cell_line_information_path = 'data/sanger_processed/metadata_ccle_sanger.csv'
drug_net_path = 'data/sanger_processed/drug_net_sanger.csv'
cleaned_drug_response_path = 'data/sanger_processed/drug_response_sanger.csv'


# Outputs path
folder_path = 'data/40DE_NMF_SANGER/'
gsea_results_path = f"{folder_path}DE_NMF_SANGER.csv"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

#Cormatrix folder path
cor_matrix_folder_path = folder_path+"cor_matrix_SANGER/"
if not os.path.exists(cor_matrix_folder_path):
    os.makedirs(cor_matrix_folder_path)

nmf_matrix_folder_path = folder_path+"nmf_matrix_SANGER/" 
if not os.path.exists(nmf_matrix_folder_path):
    os.makedirs(nmf_matrix_folder_path)
    
# 3.1 Selecting tumor types + factors
# ----------------------------------------------------------------------------------------------------------------------------------------------------

tumor_type_list = ['Colorectal Adenocarcinoma'] # See avalible tumors at the end of the script 
factors = [40] #As list, its posible to try any number
iterations = 2000
    
# 2.1 Defining functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------

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
                     
                    gsea_result, pvalues, cor_matrix = DE.run_drug_enrichment(df_drugs=df_sub_outputs, df_counts=w, drug_net=drug_net, shared_elements='name', group='moa', min_elements=5, number_of_threads=-1, return_cormatrix=True)

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
                    gsea_result, pvalues, cor_matrix = DE.run_drug_enrichment(df_drugs=df_sub_outputs, df_counts=w, drug_net=drug_net, shared_elements='name', group='moa', min_elements=5, number_of_threads=-1, return_cormatrix=True)


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