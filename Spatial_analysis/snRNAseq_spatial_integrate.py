import os
import sys
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

# Load raw and processed data
adata_FL_raw = sc.read("./save_data/FL_s1_obje.h5ad")
adata_FL = sc.read('./save_data/FL_stagate_v0.h5ad')

# Create AnnData object for spatial transcriptomics
adata_FL_st = sc.AnnData(X=adata_FL_raw[adata_FL.obs_names, :].X,
                          obs=adata_FL.obs,
                          var=adata_FL_raw.var,
                          uns=adata_FL.uns,
                          obsm=adata_FL.obsm,
                          varm=adata_FL.varm,
                          obsp=adata_FL.obsp)

adata_FL_st.var = adata_FL_st.var.set_index('features')

# Load single-cell data
adata = sc.read('./data/adata_clusters2_noraw_v1.h5ad')
adata_raw = sc.read('./data/adata_clusters2_raw.h5ad')

# Create AnnData object for single-cell data
adata_FL_sc = sc.AnnData(X=adata_raw.X,
                          obs=adata.obs,
                          var=adata_raw.var,
                          uns=adata.uns,
                          obsm=adata.obsm,
                          varm=adata_raw.varm,
                          obsp=adata.obsp)

# Filter for Frontal Lobe and visualize UMAP
adata_FL_sc = adata_FL_sc[adata_FL_sc.obs['MainLoc'] == 'Frontal lobe']
sc.pl.umap(adata_FL_sc, color='clusters2_v1', frameon=False, title='')

# Integrate single-cell and spatial transcriptomics datasets
adata_sc = adata_FL_sc.copy()

# Identify common variable names
list_of_variable_names = adata_FL_st.var.index.intersection(adata_sc.var.index)
adata_subset = adata_FL_st[:, list_of_variable_names]

list_of_variable_names = adata_sc.var.index.intersection(adata_FL_st.var.index)
adata_sc_subset = adata_sc[:, list_of_variable_names]

# Normalize and scale single-cell data
sc.pp.normalize_total(adata_sc_subset)
sc.pp.log1p(adata_sc_subset)
sc.pp.calculate_qc_metrics(adata_sc_subset, percent_top=None, inplace=True)
sc.pp.regress_out(adata_sc_subset, ['total_counts'])
sc.pp.scale(adata_sc_subset)

# Normalize and scale spatial transcriptomics data
sc.pp.normalize_total(adata_subset)
sc.pp.log1p(adata_subset)
sc.pp.calculate_qc_metrics(adata_subset, percent_top=None, inplace=True)
sc.pp.regress_out(adata_subset, ['total_counts'])
sc.pp.scale(adata_subset)

# Add technology labels
adata_subset.obs['tech'] = 'st_FL'
adata_sc_subset.obs['tech'] = 'sc_FL'

# Combine datasets
combine_adata = adata_subset.concatenate(adata_sc_subset, batch_key='dataset', batch_categories=['st', 'scrna'])

# Perform PCA and Harmony integration
npc = 30
sc.tl.pca(combine_adata, n_comps=npc)
sce.pp.harmony_integrate(combine_adata, 'tech')
print('Harmony integration finished')

# Compute neighbors and UMAP
sc.pp.neighbors(combine_adata, n_neighbors=20, use_rep='X_pca_harmony')
sc.tl.umap(combine_adata)

# Visualize combined UMAP
sc.pl.umap(combine_adata, color=['dataset', 'clusters2_v1'], legend_loc='right margin', frameon=False)

# Subset datasets for prediction
combine_adata_st = combine_adata[combine_adata.obs['dataset'] == 'st'].copy()
combine_adata_sc = combine_adata[combine_adata.obs['dataset'] == 'scrna'].copy()

# Predict cell types using Nearest Neighbors
n_neighbors = 20
neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(combine_adata_sc.obsm['X_umap'])
distances, indices = neigh.kneighbors(combine_adata_st.obsm['X_umap'])
indices_all = np.concatenate(indices)

# Create a DataFrame for predicted cell types
df = pd.DataFrame(np.take(list(adata_FL_sc.obs.clusters2_v1), indices))
df.index = adata_FL_st.obs_names

# Function to get the most frequent prediction
def max_occurrence(row):
    return row.value_counts().idxmax()

max_elements = df.apply(max_occurrence, axis=1)
new_df = pd.DataFrame(max_elements, columns=['Predict_cell_type'])

# Assign predictions to AnnData object
adata_FL_st.obs['Predict_cell_type'] = new_df['Predict_cell_type'].astype('category')

# Plot predicted cell types on t-SNE
clusters = adata_FL_st.obs['Predict_cell_type'].cat.categories
fig, axes = plt.subplots(4, 6, figsize=(20, 14))
for i, cluster in enumerate(clusters):
    ax = axes[i // 6, i % 6]
    sc.pl.tsne(adata_FL_st, color=['mclust'], groups='', ax=ax, size=2, show=False, frameon=False, title=cluster, legend_loc='')
    sc.pl.tsne(adata_FL_st, color=['Predict_cell_type'], groups=[cluster], ax=ax, size=2, show=False, frameon=False, title=cluster, legend_loc='')

# Imputation
print("Preparing for imputation...")
combine_adata_st = combine_adata[combine_adata.obs['dataset'] == 'st'].copy()
combine_adata_sc = combine_adata[combine_adata.obs['dataset'] == 'scrna'].copy()
adata_sc_imputation = adata_sc.copy()
singlecell_rawexpr = adata_sc_imputation.X
shape1 = singlecell_rawexpr.shape[1]

# Nearest Neighbors for imputation
n_neighbors = 20
neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(combine_adata_sc.obsm['X_umap'])
distances, indices = neigh.kneighbors(combine_adata_st.obsm['X_umap'])
indices_all = np.concatenate(indices)

# Impute values
adata_st_imputation = []
for i in range(0, indices_all.shape[0], 100000):
    print(f'Processing indices from {i} to {i + 100000}')
    test1 = singlecell_rawexpr[indices_all[i:i + 100000]]
    test2 = test1.reshape(-1, n_neighbors, shape1)
    test3 = np.mean(test2, axis=1)
    adata_st_imputation.append(test3)

# Concatenate imputed values
adata_imput_expr = np.concatenate(adata_st_imputation, axis=0)

# Create AnnData for imputed values
adata_st_imput = sc.AnnData(X=adata_imput_expr,
                             obs=adata_FL_st.obs,
                             uns=adata_FL_st.uns,
                             obsm=adata_FL_st.obsm)

adata_st_imput.var_names = adata_sc.var_names

# Save the imputed data
adata_st_imput.write('./save_data/FL_st_imputation.h5ad')
