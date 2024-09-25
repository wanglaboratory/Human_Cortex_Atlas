import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bbknn import bbknn
import scrublet as scr

# Print system information
import sys
print(sys.executable)
print(sys.version)
print(sys.version_info)

sc.settings.set_figure_params(dpi=150)

# Step 1: Load data from Cell Ranger
adata = sc.read_10x_mtx('./Data/adult_human_cortex/AHCaggr24/filtered_feature_bc_matrix/', cache=False)

# Calculate mitochondrial percentage and counts
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Merge sample info
info = pd.read_csv("./Data/sample_info_v1.csv")
adata.obs['cell'] = adata.obs.index
adata.obs[['barcode', 'batch']] = adata.obs['cell'].str.split('-', expand=True)
adata.obs['batch'] = pd.to_numeric(adata.obs['batch'])
adata.obs = adata.obs.merge(info, left_on='batch', right_on='number').set_index('cell')

# Save raw data
adata.write("./Data/raw_v1.h5ad")
sc.pl.violin(adata, ['n_counts'], jitter=0.4, multi_panel=False)

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[(adata.obs['n_genes'].between(600, 10000)) & (adata.obs['percent_mito'] < 0.2)]

# Doublet detection
def run_scrublet(batch_str):
    raw_filt = adata[adata.obs.batch == batch_str, :]
    scrub = scr.Scrublet(raw_filt.X, expected_doublet_rate=0.1)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    return raw_filt.obs.index[predicted_doublets]

df_cell = adata.obs.copy()
df_cell['doublet'] = 'false'

for batch in adata.obs['batch'].unique():
    filt_cell = run_scrublet(batch)
    df_cell.loc[filt_cell, 'doublet'] = 'true'
    print(f"Batch finished: {batch}")

adata.obs['doublet'] = df_cell['doublet']
adata = adata[adata.obs.doublet == 'false']

# Step 2: Dimension reduction and clustering
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack', n_comps=100)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=100)

# BBKNN and UMAP
bbknn(adata, batch_key='sample', n_pcs=80)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Save results
adata.write("./Data/HC_HVG2000_v1.h5ad")
sc.pl.umap(adata, color=['percent_mito', 'n_genes'], legend_loc='on data')
sc.pl.umap(adata, color='leiden', legend_loc='on data')

# Step 3: remove low quality leiden cluster '17','5'
adata = adata[~adata.obs.leiden.isin(['17', '5'])]
sc.tl.umap(adata)  # Recalculate UMAP after filtering
sc.tl.leiden(adata, resolution=1)  # Recalculate leiden clustering
adata.write("./Data/HC_HVG2000_filt_lowQualityCluster.h5ad")
adata.obs.to_csv('./Data/HC_HVG2000_metadata.csv')


