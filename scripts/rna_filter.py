import muon as mu
import scanpy as sc

# Load in data
adata = sc.read_h5ad(snakemake.input.rna_anndata) # type: ignore

# Filter for quality control on general values
# The lambda function sets the observation row x to be filtered the variable in 
# quotes by the value input by the snakemake params value

# Threshold below a given mitochondria percent
mu.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= snakemake.params.mito_percent_thresh)
# Threshold below a given doublet score
mu.pp.filter_obs(adata, 'doublet_score', lambda x: x < snakemake.params.doublet_thresh)   
# Threshold above a given genes per cell threshold
adata = adata[adata.obs['n_genes_by_counts'] > snakemake.params.min_genes_per_cell].copy()
# Threshold below a given ribosome threshold
mu.pp.filter_obs(adata, 'pct_counts_rb', lambda x: x <= snakemake.params.ribo_percent_thresh)

if adata.n_obs != 0:
    # Normalize data
    sc.pp.normalize_total(adata)

    # Save the CPM data
    adata.layers['cpm']=adata.X.copy() 

    # Logarithmize the data
    sc.pp.log1p(adata, layer='cpm')

    # Save the normalized-log data
    adata.layers['data']=adata.X.copy() 

# Write out filtered anndata object
adata.write(filename=snakemake.output.rna_anndata, compression='gzip') # type: ignore