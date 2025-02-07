import decoupler as dc
import pandas as pd 
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Inputs
netPath = "data/prism_processed/drug_net_prism.csv"
cormatrixPath = "data/DE_TF_PRISM/cor_matrix/Colorectal_Adenocarcinoma_cormatrix.csv"


# Outputs 
plotPath = "figures/enrichmentPlots/STAT3PR.png"
directory = os.path.dirname(plotPath)
if not os.path.exists(directory):
    os.makedirs(directory)

df = pd.read_csv(cormatrixPath, index_col=0)
net = pd.read_csv(netPath)

cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["#47a4c3d0", "white", "#792020dd"])

plt.rcParams['font.family'] = 'Times New Roman'
fig, arr = dc.plot_running_score(df, "STAT3", net, 'EGFR inhibitor', source='moa', target='broad_id', cmap=cmap_custom, figsize=(5, 5), dpi=100, return_fig=True, save=None)
plt.title("Factor11")
plt.savefig(plotPath, dpi=1000)
