import decoupler as dc
import pandas as pd 
import os 
import matplotlib as mpl
import matplotlib.pyplot as plt

# Inputs
hmatrixpath = 'data/DE_NMF_SANGER/nmf_matrix_SANGER/Colorectal_Adenocarcinoma_40_h_matrix.csv'
desaresutlpath = 'data/DE_NMF_SANGER/DE_NMF_SANGER.csv'
logfctablepath = 'data/diffexp/DEA_cohort.csv'

# Outputs 
folderpath = 'figures/lfcvolcanos/'
os.makedirs(folderpath, exist_ok=True)


with open(hmatrixpath) as file:
    hmatrix = pd.read_csv(file, index_col=0)
with open(desaresutlpath) as file:
    dsea = pd.read_csv(file, index_col=0)
with open(logfctablepath) as file:
    logfcdf = pd.read_csv(file, index_col=0)
    
mpl.rcParams['font.family'] = 'Times New Roman'


tumor_type = "Colorectal_Adenocarcinoma"
# Outputs 
folderpath = 'figures/lfcvolcanos/'
os.makedirs(folderpath, exist_ok=True)


with open(desaresutlpathp) as file:
    desaprism = pd.read_csv(file, index_col=0)
with open(desaresutlpath) as file:
    dsea = pd.read_csv(file, index_col=0)
with open(logfctablepath) as file:
    logfcdf = pd.read_csv(file, index_col=0)
    
mpl.rcParams['font.family'] = 'Times New Roman'


tumor_type = "Colorectal_Adenocarcinoma"

filterData = dsea[
    (dsea['Tumor type'] == tumor_type) &
    (dsea['MoA'] == 'EGFR signaling') &
    (dsea['GSEA value'].abs() > 0.65) &
    (dsea['p-value'] < 0.05)
]

ilterDataPrism = desaprism[
    (desaprism['Tumor type'] == tumor_type) &
    (desaprism['MoA'] == 'EGFR inhibitor') &
    (desaprism['GSEA value'].abs() > 0.65) &
    (desaprism['p-value'] < 0.05)
]




factors_to_keep = filterData['Factors'].unique()

Hmatrix_filtered = hmatrix[hmatrix.index.isin(factors_to_keep)]

Hmatrix_df = Hmatrix_filtered.copy()
Hmatrix_df['Factor'] = Hmatrix_df.index

Hmatrix_long = pd.melt(Hmatrix_df, id_vars='Factor', 
                       var_name='Gene', value_name='weight')

Hmatrix_long['percentile_90'] = Hmatrix_long.groupby('Factor')['weight'] \
                                            .transform(lambda x: x.quantile(0.9))
Hmatrix_top10 = Hmatrix_long[Hmatrix_long['weight'] >= Hmatrix_long['percentile_90']]
Hmatrix_top10['Gene'] = [e.split(" ")[0] for e in Hmatrix_top10['Gene']]

for e in Hmatrix_top10['Factor'].unique():
    subgenes = Hmatrix_top10['Gene'][Hmatrix_top10['Factor']==e].unique()
    sublogfcdf = logfcdf[logfcdf['Row.names'].isin(subgenes)]
        
    logfc_series = sublogfcdf.set_index('Row.names')['log2FoldChange']
    pval_series  = sublogfcdf.set_index('Row.names')['pvalue']

    contrast_label = "Resistentes"

    logfc_single = pd.DataFrame([logfc_series.to_dict()], index=[contrast_label])
    pvals_single = pd.DataFrame([pval_series.to_dict()], index=[contrast_label])
    

    fig = dc.plot_volcano(
        logFCs=logfc_single,   
        pvals=pvals_single,    
        contrast=contrast_label,  
        top=0,
        sign_thr=0.05,
        lFCs_thr=0.5,
        return_fig=True,
    )
    ax = fig.axes[0]
    ax.set_title("")
    
    # Agregar un título en la parte inferior del gráfico
    # La posición (0.5, 0.02) coloca el texto centrado horizontalmente y cerca del borde inferior.
    fig.text(0.5, 0.95, e, ha='center', va='top', fontsize=12, fontname='Times New Roman')
    fig.set_size_inches(6,4)

    fig.savefig((folderpath+e+'.png'), dpi=1000)
    plt.close(fig)