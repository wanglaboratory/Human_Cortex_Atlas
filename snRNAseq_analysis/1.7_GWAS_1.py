import scanpy as sc
import pandas as pd
import cellex

# Set figure parameters for better visibility
sc.settings.set_figure_params(dpi=120)

# Load the raw data
adata_raw = sc.read('./data/All_region_cortex_raw.h5ad')

# Create a combined identifier for Main Location and Cluster
adata_raw.obs['MainLoc_clusters2_v1'] = adata_raw.obs['MainLoc'] + '_' + adata_raw.obs['clusters2_v1']

# Create a DataFrame from the AnnData object, transposing the expression matrix
data = pd.DataFrame(adata_raw.X.T, index=adata_raw.var['gene_ids'], columns=adata_raw.obs_names)

# Extract the metadata for the ESObject
metadata = adata_raw.obs['MainLoc_clusters2_v1']

# Create an ESObject for enrichment analysis
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)

# Compute enrichment statistics
eso.compute(verbose=True)

# Display the top results of the enrichment analysis
print(eso.results["esmu"].head())

# Save the results to a CSV file
eso.save_as_csv(path='./data/', file_prefix='MainLoc_clusters2_v1_human_brain_cells', verbose=True)
