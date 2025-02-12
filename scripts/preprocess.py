import os
import pandas as pd
import scipy
import scanpy as sc

"""
This script takes the output of either Cellbender or CellRanger and processes into a Anndata object
with the parameters upon which quality control filtering can be done.
"""

# Read the samples table once
samples = pd.read_csv(snakemake.input.input_table)

# Extract the metadata for the specific sample in one step
metadata = samples[samples['Sample_ID'] == snakemake.params.sample].iloc[0]

"""Preprocess the RNA data"""

# Read the single-cell data
adata = sc.read_10x_h5(snakemake.input.rna_anndata)

# Ensure unique variable names (MUST BE DONE FIRST!)
adata.var_names_make_unique()

# Make a raw counts layer
adata.layers['counts'] = adata.X.copy()
adata.raw = adata

# Add mitochondrial and ribosomal markers
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars = ['rb', 'mt'], 
    percent_top = None, 
    log1p = False, 
    inplace = True
    )

# Run scrublet to identify doublets
sc.external.pp.scrublet(
    adata, 
    expected_doublet_rate = (adata.n_obs / 1000) * 0.008)
adata.obs.drop('predicted_doublet', axis=1, inplace=True)
adata.obs['cell_barcode'] = adata.obs_names

# Add metadata to the AnnData object directly from the metadata dataframe
metadata_dict = {
    'batch'             : metadata['Use_batch'],
    'sex'               : metadata['Sex'],
    'age'               : metadata['Age'],
    'pmi'               : metadata['PMI'],
    'ethnicity'         : metadata['Ethnicity'],
    'race'              : metadata['Race'],
    'brain_bank'        : metadata['Brain_bank'],
    'homogenization'    : metadata['Homogenizing_batch'],
    'library'           : metadata['Library_batch'],
    'seq'               : metadata['Sequencing_batch'],
    'sample'            : metadata['Sample_ID'],
    'short diagnosis'   : metadata['Short Diagnosis']
}

for key, value in metadata_dict.items():
    adata.obs[key] = value

# Normalize data
sc.pp.normalize_total(adata)

# Save the CPM data
adata.layers['cpm'] = adata.X.copy() 

# Logarithmize the data
sc.pp.log1p(adata)

# Save the normalized-log data
adata.layers['data'] = adata.X.copy() 

# Calculate cell cycle()
cell_cycle_genes = [x.strip() for x in open('/data/CARD_singlecell/SN_atlas/input/lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
try:
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes, 
        g2m_genes=g2m_genes)
except:
    print("can't map genes")

# Calculate nearest neighbors
sc.pp.neighbors(adata)

# Calculate UMAP
sc.tl.umap(adata)

# Save the AnnData object
adata.write(
    filename = snakemake.output.rna_anndata, 
    compression = 'gzip')
