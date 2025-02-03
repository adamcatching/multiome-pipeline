import anndata as ad
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Open rna
adata = sc.read_h5ad(snakemake.input.rna_anndata)

# Read in parameters
cell_type = snakemake.params.cell_type
disease_name = snakemake.params.disease
control_name = snakemake.params.control

# Subset the data for the cell type
adata = adata[adata.obs['cell_type'] == cell_type].copy()

# Get the list of cells enriched for disease state
cell_type_control = cell_type_atac.obs['Primary Diagnosis'] == control_name
cell_type_disease = cell_type_atac.obs['Primary Diagnosis'] == disease_name

# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='sample',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=10
)

pdata.obs['age'] = pdata.obs.age.astype('float')

# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)

# Abreviate diagnosis
pdata.obs['diagnosis'] = pdata.obs['Primary Diagnosis']

# Normalize ages
ages = pdata.obs.age
pdata.obs['normalage'] = (ages-np.min(ages))/(np.max(ages)-np.min(ages))-.5

dc.get_metadata_associations(
    pdata,
    obs_keys = ['normalage', 'diagnosis', 'psbulk_n_cells', 'psbulk_counts'],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)

pdata_genes = dc.filter_by_expr(
    pdata, 
    group='diagnosis', 
    min_count=10, 
    min_total_count=15
    )

pdata = pdata[:, pdata_genes].copy()

inference = DefaultInference(n_cpus=8)

dds = DeseqDataSet(
    adata=astrocytes,
    design_factors=['normalage', 'diagnosis', 'batch'],
    refit_cooks=True,
    inference=inference,
)

# Compute LFCs
dds.deseq2()

# Extract contrast between normal and DLB
stat_res = DeseqStats(
    dds,
    contrast=["diagnosis", 'DLB', 'control'],
    inference=inference,
)

# Compute Wald test
stat_res.summary()

# Extract results
DGE_results_df = stat_res.results_df
DGE_results_df['-log10_padj'] = -np.log10(DGE_results_df['padj'])
DGE_results_df.to_csv(snakemake.output.output_figure)

# Plot 
dc.plot_volcano_df(
    astrocyte_DLB_results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    lFCs_thr=1,
    sign_thr=1e-2,
    figsize=(4, 4),
    dpi=600,
    save=snakemake.output.output_figure
)