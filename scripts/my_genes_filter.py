# IMPORTS
import os
import sys
import scanpy as sc

# DEFINE - file paths
input_path = snakemake.input.work_dir 

# LOAD
# - Params
params_dict = {
    'min_genes_per_cell'    : snakemake.params.min_genes_per_cell,
}

# - Data
adata_03 = sc.read_h5ad(snakemake.input.merged_rna_anndata_pre_filter) 

# RUN

# # CHECK
# # LOOP - print all params in dict
# for key, value in params_dict.items():
#     print(f'{key} = {value}')

# LOGGING
# - Get current date and time
formatted_time = snakemake.params.formatted_time

# # DEFINE - file paths
# merged_rna_anndata_pre_filter = snakemake.input.merged_rna_anndata_pre_filter
# file_name = os.path.basename(merged_rna_anndata_pre_filter)
# output_path = os.path.join(input_path, 'src', 'output', f'file_name-{formatted_time}')
# print('output_path =', output_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#* TODO - rewrite for Python script 
# TEST
# DEFINE - file pathsname
output_path = os.path.join('src', 'output', os.path.basename(snakemake.input.merged_rna_anndata_pre_filter).replace('.h5ad', '-ef-TEST.h5ad'))

# SAVE - new anndata object to a file w/the threshold filtered
adata = adata_03[adata_03.obs['n_genes_by_counts'] > 250].copy()
adata.write_h5ad(output_path, compression = 'gzip')

# CHECK - file path
print('Saved file to:', output_path)
# 03_filtered_anndata_rna-ef-TEST.h5ad