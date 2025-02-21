import decoupler as dc
import pandas as pd 
import os 
import matplotlib as mpl
import matplotlib.pyplot as plt

# Inputs
desaresutlpathp = 'data/DE_TF_PRISM/DE_TF_PRISM.csv'
desaresutlpath = 'data/DE_TF_SANGER/DE_TF_SANGER.csv'
logfctablepath = 'data/diffexp/DEA_cohort.csv'

# Outputs 
folderpath = 'figures/ZZlfcvolcanos/'
os.makedirs(folderpath, exist_ok=True)


with open(desaresutlpathp) as file:
    desaprism = pd.read_csv(file, index_col=0)
with open(desaresutlpath) as file:
    dsea = pd.read_csv(file, index_col=0)
with open(logfctablepath) as file:
    logfcdf = pd.read_csv(file)
    
mpl.rcParams['font.family'] = 'Times New Roman'


tumor_type = "Colorectal_Adenocarcinoma"


filterData = dsea[
    (dsea['Tumor type'] == tumor_type) &
    (dsea['MoA'] == 'EGFR signaling') &
    (dsea['GSEA value'].abs() > 0.65) &
    (dsea['p-value'] < 0.05)
]

filterDataPrism = desaprism[
    (desaprism['Tumor type'] == tumor_type) &
    (desaprism['MoA'] == 'EGFR inhibitor') &
    (desaprism['GSEA value'].abs() > 0.65) &
    (desaprism['p-value'] < 0.05)
]



# Intersección positiva (GSEA value > 0)
intersect_positive_TF = set(
    filterData[filterData['GSEA value'] > 0]['TF']
).intersection(set(
    filterDataPrism[filterDataPrism['GSEA value'] > 0]['TF']
))

# Intersección negativa (GSEA value < 0)
intersect_negative_TF = set(
    filterData[filterData['GSEA value'] < 0]['TF']
).intersection(set(
    filterDataPrism[filterDataPrism['GSEA value'] < 0]['TF']
))
net = dc.get_collectri(organism='human')

positive_targets = net['target'][net['source'].isin(intersect_positive_TF)].unique()
negative_targets = net['target'][net['source'].isin(intersect_negative_TF)].unique()
 
sublogfcdf = logfcdf[logfcdf['Row.names'].isin(positive_targets)]
    
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
    return_fig=True
)

ax = fig.axes[0]
ax.set_title("")

fig.text(0.5, 0.95, 'TF positivos', ha='center', va='top', fontsize=12, fontname='Times New Roman')
fig.set_size_inches(5,3.6)

plt.savefig((folderpath+'tfsn33gsitivos.png'), dpi=1500)
