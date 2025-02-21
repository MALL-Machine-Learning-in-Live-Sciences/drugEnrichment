# 1. Preprocessing PRISM
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

# PRISM repurposing dataset
# downloaded from https://depmap.org/portal/download/all/
# PRISM Repurposing 19Q4

prism_inptus_path = "extdata/PRISM_19Q4/secondary-screen-dose-response-curve-parameters.csv"

# 1.3 Outputs
# ---------------------------------------------------------------------------------------------------------------------------------------------------

folder_path = 'data/prism_processed/'

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Matrix with cell-lines*Gene_expression
cleaned_counts_path = 'data/prism_processed/counts_matched_prism.csv'
# Metadata of cells used in prism
cleaned_cell_line_information_path = 'data/prism_processed/metadata_ccle_prism.csv'
# Df matching broad_id+name+moa
drug_net_path = 'data/prism_processed/drug_net_prism.csv'
# Matrix with cell-lines*auc
cleaned_drug_response_path = 'data/prism_processed/drug_response_prism.csv'



# 1.4 Cleaning ccle
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Reading dataframes:
with open(ccle_counts_path) as file:
    df_counts = pd.read_csv(file, index_col=0)
# Reading dataframes:
with open(ccle_cell_information) as file:
    df_cellinfo = pd.read_csv(file, index_col=0)
    # Reading dataframes:
with open(prism_inptus_path) as file:
    df_outputs = pd.read_csv(file)
    
shared_cells = set(df_outputs['depmap_id']).intersection(df_counts.index)

# Cleaning counts df
df_counts = df_counts[df_counts.index.isin(shared_cells)]
df_counts.columns = [col.split(" ")[0] for col in df_counts.columns]

# Cleaning metadata df
df_cellinfo = df_cellinfo[df_cellinfo.index.isin(shared_cells)]

# Saving dfs
df_counts.to_csv(cleaned_counts_path)
df_cellinfo.to_csv(cleaned_cell_line_information_path)



# 1.5 Cleaning prism
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# Cleaning prism df
df_outputs = df_outputs[df_outputs['depmap_id'].isin(shared_cells)]
df_outputs['screen_id'] = pd.Categorical(df_outputs['screen_id'], categories=["MTS010", "MTS006", "MTS005", "HTS002"], ordered=True)
df_outputs_sorted_sorted = df_outputs.sort_values(by='screen_id')

seen_broad_ids = set()
result_rows = []

# Filter by screen id
for screen_id in ["MTS010", "MTS006", "MTS005", "HTS002"]:

    current_level = df_outputs_sorted_sorted[df_outputs_sorted_sorted['screen_id'] == screen_id]
    new_broad_ids = current_level[~current_level['broad_id'].isin(seen_broad_ids)]
    seen_broad_ids.update(new_broad_ids['broad_id'])
    result_rows.append(new_broad_ids)

df_outputs_filtered = pd.concat(result_rows, ignore_index=True)

# Transform to matrix
df_pivot = df_outputs_filtered.groupby(['depmap_id', 'broad_id'])['auc'].max().unstack()

# Drug_net
drug_net = df_outputs_filtered[['broad_id','name','moa']].drop_duplicates()
drug_net = df_outputs_filtered[['broad_id', 'name', 'moa']].drop_duplicates()
drug_net['moa'] = drug_net['moa'].str.split(',')  
drug_net = drug_net.explode('moa')  
drug_net['moa'] = drug_net['moa'].str.strip() 

# Saving dfs
df_pivot.to_csv(cleaned_drug_response_path)
drug_net.to_csv(drug_net_path)
