import Drug_Enrichment as DE
import decoupler as dc
import pandas as pd
import os 
import anndata as ad


df_resultado = pd.DataFrame()
#Select tumor types(TT), Avalible TT are listed at the end of the code 
type_list =['Colorectal Adenocarcinoma']

# 2. Gsea table
# ----------------------------------------------------------------------------------------------------------------------------------------------------

# Inputs path:
counts_path = 'extdata/DepMap_23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv'
cell_line_information_path = 'extdata/DepMap_23Q2/Model.csv'
cleaned_cell_line_information_path = 'data/prism_processed/metadata_ccle_prism.csv'
drug_net_path = 'data/prism_processed/drug_net_prism.csv'
cleaned_drug_response_path = 'data/prism_processed/drug_response_prism.csv'

# Outputs path
folder_path = 'data/DE_TF_PRISM/'
gsea_results_path = f"{folder_path}DE_TF_PRISM.csv"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

#Cormatrix folder path
cor_matrix_folder_path = folder_path+"cor_matrix/"
if not os.path.exists(cor_matrix_folder_path):
    os.makedirs(cor_matrix_folder_path)

with open(counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
with open(cell_line_information_path) as file:
    df_metadata = pd.read_csv(file, index_col=0)
with open(cleaned_cell_line_information_path) as file:
    df_cleaned_metadata = pd.read_csv(file, index_col=0)
with open(drug_net_path) as file:
    drug_net = pd.read_csv(file, index_col=0)
with open(cleaned_drug_response_path) as file:
    df_outputs = pd.read_csv(file, index_col=0)


tf_net = dc.get_collectri(organism='human', split_complexes=False)



# 2.3 GSEA splicing by organ:
# ----------------------------------------------------------------------------------------------------------------------------------------------------

for organ in df_metadata['OncotreeLineage'].unique():
    if organ in type_list:
        # Select all cell-lines avalible at CCLE of one TT
        cell_lines = df_metadata.index[df_metadata['OncotreeLineage']==organ].to_list()
        df_counts_organ = df_counts[df_counts.index.isin(cell_lines)]
        df_counts_organ.columns = [col.split(" ")[0] for col in df_counts_organ.columns]
        
        # TF activity inference
        adata = ad.AnnData(X=df_counts_organ.values, obs=pd.DataFrame(index=df_counts_organ.index), var=pd.DataFrame(index=df_counts_organ.columns))
        dc.run_ulm(mat=adata,net=tf_net,source='source',target='target',weight='weight',verbose=True,use_raw = False,batch_size=10000, min_n=5)
        df_tf = adata.obsm['ulm_estimate']
        
        # Match both df
        sub_celllines = df_cleaned_metadata.index[df_cleaned_metadata['OncotreeLineage']==organ].to_list()      
        df_sub_outputs = df_outputs[df_outputs.index.isin(sub_celllines)]
        df_sub_tf = df_tf[df_tf.index.isin(sub_celllines)]
        try:
            if len(df_sub_outputs)>10 and len(df_sub_tf)>10:
                gsea_result, pvalues, cor_matrix = DE.run_drug_enrichment(df_drugs=df_sub_outputs, df_counts=df_sub_tf, drug_net=drug_net, shared_elements='broad_id', group='moa', min_elements=5, number_of_threads=-1, return_cormatrix=True)

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
    if organ in type_list:
        cell_lines = df_metadata.index[df_metadata['OncotreePrimaryDisease']==organ].to_list()
        df_counts_organ = df_counts[df_counts.index.isin(cell_lines)]
        df_counts_organ.columns = [col.split(" ")[0] for col in df_counts_organ.columns]

        adata = ad.AnnData(X=df_counts_organ.values, obs=pd.DataFrame(index=df_counts_organ.index), var=pd.DataFrame(index=df_counts_organ.columns))
        dc.run_ulm(mat=adata,net=tf_net,source='source',target='target',weight='weight',verbose=True,use_raw = False,batch_size=10000, min_n=5)
        df_tf = adata.obsm['ulm_estimate']
        
        sub_celllines = df_cleaned_metadata.index[df_cleaned_metadata['OncotreePrimaryDisease']==organ].to_list()      
        df_sub_outputs = df_outputs[df_outputs.index.isin(sub_celllines)]
        df_sub_tf = df_tf[df_tf.index.isin(sub_celllines)]
        
        try:
            if len(df_sub_outputs)>10 and len(df_sub_tf)>10:
                gsea_result, pvalues, cor_matrix = DE.run_drug_enrichment(df_drugs=df_sub_outputs, df_counts=df_sub_tf, drug_net=drug_net, shared_elements='broad_id', group='moa', min_elements=5, number_of_threads=-1, return_cormatrix=True)
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
