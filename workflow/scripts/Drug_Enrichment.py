import decoupler as dc
import anndata as ad

# 0.1 Defining functions
# ----------------------------------------------------------------------------------------------------------------------------------------------------
def compute_correlation(drug_col, stats, df_counts, df_drugs, min_len):
    correlations = []
    df_suboutputs = df_drugs[drug_col].dropna()
    for gene in df_counts.columns:
        if len(df_suboutputs.unique())>min_len:
            df_subgene = df_counts[gene].dropna()
            shared_cells = df_suboutputs.index.intersection(df_subgene.index)       
            correlation, _ = stats.pearsonr(df_suboutputs[shared_cells].sort_index().values,
                                    df_subgene[shared_cells].sort_index().values)
            correlations.append((gene, correlation))
        else:
            correlations.append((gene, None))

    return drug_col, correlations
    
def correlation_drugs(df_drugs, df_counts, threads, min_len):  
    '''
    Inputs: drug matrix (cell-response) and expresión matrix (cell-expression)
    Outputs: correlation matrix: gene-drug
    '''
    from joblib import Parallel, delayed
    import pandas as pd
    from scipy import stats  


    correlation_matrix = pd.DataFrame(index=df_drugs.columns, columns=df_counts.columns)

    results = Parallel(n_jobs=threads)(delayed(compute_correlation)(drug_col, stats, df_counts, df_drugs, min_len) for drug_col in df_drugs.columns)

    for drug_col, correlations in results:
        for gene, correlation in correlations:
            correlation_matrix.at[drug_col, gene] = correlation

    return correlation_matrix

def gsea_enrichment(matrix, net, shared_element, net_group, n_min):
    '''
    Inputs: 
        matrix: correlation matrix 
        net: long dataset with features and groups  
        shared_element: commun column between matrix-net
        net_pathway: name of the column with the groups in the net
        n_min: min elements to consider 
        
    Outputs: gsea matrix and p-values
    '''
    
    matrix = matrix.dropna(0)
    shared_elements = set(matrix.index).intersection(set(net[shared_element]))

    adata = ad.AnnData(matrix.loc[shared_elements, :].T)

    result  = dc.run_gsea(mat=adata,net=net[net[shared_element].isin(shared_elements)],source=net_group,target=shared_element,min_n=n_min,use_raw=False)
    result_df_gsea = adata.obsm['gsea_estimate'] 
    p_value_gsea = adata.obsm['gsea_pvals']
    
    return result_df_gsea, p_value_gsea

def run_drug_enrichment(df_drugs, df_counts, drug_net, shared_elements, group, min_elements=5, number_of_threads=-1, return_cormatrix=False, min_len=1):
    '''
    Computes drug enrichment based on gene expression and drug response data.

    Intputs
        df_drugs : pandas.DataFrame Matrix of drug response data with dimensions (n_cells, m_drug_response).
        df_counts : pandas.DataFrame Matrix of gene expression data with dimensions (n_cells, m_genes).
        drug_net : pandas.DataFrame Long-format dataset containing groups of drugs.
        shared_elements : str. Column name in `drug_net` shared with the df_drugs columns
        group : str. Column name in `drug_net` containing the groups.
        min_elements : int. Minimum number of elements required to compute drug enrichment.
        number_of_threads : int. Number of threads to use for parallel computation of the correlation matrix.
        return_cormatrix : bool, optional. If `True`, the function will also return the correlation matrix. Defaults to `False`.

    Outputs
        Drug enrichment matrix 
        P-values for enrichment matrix 
        Correlation matrix (optional if `return_cormatrix=True`
    '''
    cormatrix = correlation_drugs(df_drugs, df_counts, number_of_threads, min_len)
    df_results, df_pvalue = gsea_enrichment(matrix=cormatrix, net=drug_net, shared_element=shared_elements, net_group=group, n_min=min_elements)
    if return_cormatrix==False:
        return df_results, df_pvalue
    else:
        return df_results, df_pvalue, cormatrix


    