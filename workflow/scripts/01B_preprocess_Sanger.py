# 1. Preprocessing Sanger
# ---------------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd
import os 

# 1.2. Inputs
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# CCLE counts
# DepMap Public 23Q2
# downloaded from: https://depmap.org/portal/download/all/

ccle_counts_path = "extdata/DepMap_23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
ccle_cell_information = "extdata/DepMap_23Q2/Model.csv"

# GDSC dataset 
# downloaded from https://depmap.org/portal/download/all/

sanger_inputs_path = "extdata/Sanger/sanger-dose-response.csv"

# GDSC compunds
# downloaded from https://www.cancerrxgene.org/downloads/bulk_download
drugs_info_path = "extdata/Sanger/screened_compounds_rel_8.5.csv"


# 1.3 Outputs
# ---------------------------------------------------------------------------------------------------------------------------------------------------

folder_path = 'data/sanger_processed/'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Matrix with cell-lines*Gene_expression
cleaned_counts_path = 'data/sanger_processed/counts_matched_sanger.csv'
# Metadata of cells used in prism
cleaned_cell_line_information_path = 'data/sanger_processed/metadata_ccle_sanger.csv'
# Df matching broad_id+name+moa
drug_net_path = 'data/sanger_processed/drug_net_sanger.csv'
# Matrix with cell-lines*auc
cleaned_drug_response_path = 'data/sanger_processed/drug_response_sanger.csv'

# Reading dataframes:
with open(ccle_counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
# Reading dataframes:
with open(ccle_cell_information) as file:
    df_cellinfo = pd.read_csv(file, index_col=0)
# Reading dataframes:
with open(sanger_inputs_path) as file:
    df_outputs = pd.read_csv(file)
# Reading dataframes:
with open(drugs_info_path) as file:
    df_drugs = pd.read_csv(file)
    
 
# 1.4 Cleaning Sanger
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Filter drugs with known targets
df_drugs = df_drugs.dropna(subset='TARGET')
df_drugs = df_drugs[df_drugs['TARGET_PATHWAY']!='Unclassified']
df_drugs = df_drugs[df_drugs['TARGET_PATHWAY']!='Other']
df_drugs = df_drugs[df_drugs['DRUG_NAME'].str.lower().isin(df_drugs['DRUG_NAME'].str.lower())]
df_drugs = df_drugs[['DRUG_NAME','TARGET_PATHWAY']].drop_duplicates()

df_drugs['DRUG_NAME'] = df_drugs['DRUG_NAME'].str.lower()
df_outputs['DRUG_NAME'] = df_outputs['DRUG_NAME'].str.lower()
#rename drugs columns
df_drugs.rename(columns={'DRUG_NAME': 'name', 'TARGET_PATHWAY': 'moa'}, inplace=True)

#Filter sanger by drug net
shared_drugs = set(df_outputs['DRUG_NAME']).intersection(set(df_drugs['name']))
df_drugs = df_drugs[df_drugs['name'].isin(shared_drugs)]
df_outputs = df_outputs[df_outputs['DRUG_NAME'].isin(shared_drugs)]

#Filter by dataset release 
seen_drugs = set()  
result_rows = [] 

for data in ["GDSC2", "GDSC1"]:
    current_level = df_outputs[df_outputs['DATASET'] == data]
    
    new_drug = current_level[~current_level['DRUG_NAME'].isin(seen_drugs)]
    
    seen_drugs.update(new_drug['DRUG_NAME'])
    result_rows.append(new_drug)

result_df = pd.concat(result_rows, ignore_index=True)
result_df['BROAD_ID'] = result_df['BROAD_ID'].str.split(',', n=1).str[0]

#Fitlter by cells in CCLE
shared_cells = set(df_outputs['ARXSPAN_ID']).intersection(set(df_counts.index))
result_df = result_df[result_df['ARXSPAN_ID'].isin(shared_cells)]

# Transform to matrix
df_pivot = result_df.groupby(['BROAD_ID', 'DRUG_NAME'])['auc'].max().unstack()
df_pivot.index.name = 'depmap_id'

# Saving dfs
df_pivot.to_csv(cleaned_drug_response_path)
df_drugs.to_csv(drug_net_path)   

 
# 1.4 Cleaning CCLE
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Cleaning counts df
df_counts = df_counts[df_counts.index.isin(shared_cells)]
df_counts.columns = [col.split(" ")[0] for col in df_counts.columns]

# Cleaning metadata df
df_cellinfo = df_cellinfo[df_cellinfo.index.isin(shared_cells)]

# Saving dfs
df_counts.to_csv(cleaned_counts_path)
df_cellinfo.to_csv(cleaned_cell_line_information_path)



