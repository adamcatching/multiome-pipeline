#!/usr/bin/env python
# coding: utf-8

# STARTED:      1/16/25
# 
# LAST UPDATED: 2/20/25
# 
# By Eugene Fong
# 
# # GOAL(S)
# 
# - EDA on the SN data
# 
# # TODO
# 
# - make a new conda env?
# ```
# conda update -n base -c conda-forge conda
# ```
# 
# # IMPORTS

# In[26]:


import numpy as np
import os
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import missingno as msno

# from scipy.sparse import csr_matrix
# print(ad.__version__)


# # LOAD DATA
# 
# Actual source:
# ```
# data_dir = '/data/CARD_singlecell/Brain_atlas/SN_Multiome/'
# 
# data_dir+'batch{batch}/Multiome/{sample}-ARC/outs/01_{sample}_anndata_object_atac.h5ad', 
# ```
# 
# SN controls table:
# ```
# /vf/users/CARD_singlecell/SN_control_atlas/input/SN_control_samples.csv
# ```

# In[27]:


# # DEFINE FILE PATHS - CHANGE THIS!!!
input_file = '/vf/users/CARD_singlecell/SN_control_atlas/multiome-pipeline/input/SN_control_samples.csv'
# input_file = os.path.join('..', 'input', 'SN_control_samples.csv')

# CHECK
print('input_file =', input_file)


# In[28]:


# LOAD - matrix file into a Pandas DF
df = pd.read_csv(input_file) 
df


# In[29]:


# DROP - the first 2 columns
df = df.drop(df.columns[:2], axis=1)
df


# ### FIX - typo

# In[30]:


print(df['Brain_bank'].unique())


# In[31]:


# CHECK - how many typos?
len(df[df['Brain_bank'] == 'Havard'])


# In[32]:


# Replace the wrong spelling
df['Brain_bank'] = df['Brain_bank'].str.replace('Havard', 'Harvard')

# CHECK
len(df[df['Brain_bank'] == 'Havard'])


# In[33]:


# CHANGE - name is just too long to fit on graphs

# Add a next line 
df['Brain_bank'] = df['Brain_bank'].str.replace('Human Brain and Spinal Fluid Resource Center', 'Human Brain and Spinal \nFluid Resource Center')

# CHECK
print(df['Brain_bank'].unique())


# ### METADATA
# 
# # #? Q) Is there a data dictionary somewhere?
# 

# In[34]:


print(f'Shape of DF: \n{df.shape[0]} rows X {df.shape[1]} columns')


# In[35]:


df.describe()


# In[36]:


df.info()


# In[37]:


df.drop_duplicates()
for column in df.columns:
    unique_values = df[column].unique()
    print(f'Unique values in {column}: \n{unique_values}\n')


# ### PLOTTING
# 
# CHECK - for missing data 1st

# In[38]:


# PLOT - missing data using missingno
msno.matrix(df)
plt.show()


# In[39]:


msno.heatmap(df)
plt.show()


# In[40]:


msno.bar(df)
plt.show()


# ### PLOTTING
# 
# Pick some col's to look at

# In[41]:


# Set Seaborn style
sns.set(style = 'whitegrid')

column_name = 'Age'

# PLOT - histogram
sns.histplot(df[column_name], bins = 5, kde = False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name} (years)')
plt.ylabel('Frequency')

# SHOW - plot
plt.show()


# #? Q) What units are the PMIs?
# 
# #? A) Hours

# In[42]:


# Set Seaborn style
sns.set(style = 'whitegrid')

column_name = 'PMI'

# PLOT - histogram
sns.histplot(df[column_name], bins = 5, kde = False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name} (yrs)')
plt.ylabel('Frequency')

# SHOW - plot
plt.show()


# In[43]:


# Set Seaborn style
sns.set(style='whitegrid')

column_name = 'Ethnicity'

# PLOT - histogram
sns.histplot(df[column_name], bins=5, kde=False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name}')
plt.ylabel('Frequency')

# ROTATE - x-axis labels
plt.xticks(rotation = 45)

# SHOW - plot
plt.show()


# In[44]:


# Set Seaborn style
sns.set(style='whitegrid')

column_name = 'Race'

# PLOT - histogram
sns.histplot(df[column_name], bins=5, kde=False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name}')
plt.ylabel('Frequency')

# ROTATE - x-axis labels
plt.xticks(rotation = 45)

# SHOW - plot
plt.show()


# In[ ]:


# Set Seaborn style
sns.set(style='whitegrid')

column_name = 'Brain_bank'

# PLOT - histogram
sns.histplot(df[column_name], bins=5, kde=False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name}')
plt.ylabel('Frequency')

# SEPARATE - x-axis labels
plt.xticks(rotation=45, ha='right')

# ROTATE - x-axis labels
plt.xticks(rotation = 45)

# SHOW - plot
plt.show()


# In[46]:


# Set Seaborn style
sns.set(style='whitegrid')

column_name = 'Sex'

# PLOT - histogram
sns.histplot(df[column_name], bins=5, kde=False)

# ADD - title and labels
plt.title(f'{column_name} Distribution')
plt.xlabel(f'{column_name}')
plt.ylabel('Frequency')

# SHOW - plot
plt.show()


# In[47]:


#* TODO - change into a pie chart?


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ANNDATA SECTION
# 
# SOURCE: [https://anndata.readthedocs.io/en/stable/tutorials/notebooks/getting-started.html](https://anndata.readthedocs.io/en/stable/tutorials/notebooks/getting-started.html)

# ### DEFINE PATHS

# In[48]:


input_file = os.path.join('/data/CARD_singlecell/Brain_atlas/SN_Multiome/control_atlas/05_annotated_anndata_rna.h5ad')


# ### INITIALIZING ANNDATA OBJECT
# 
# ADAM: try importing with `sc.read.h5ad`

# In[54]:


adata = sc.read_h5ad(input_file)
adata


# In[56]:


adata.obs


# In[60]:


# FIX - the wrong spelling
adata.obs['brain_bank'] = adata.obs['brain_bank'].str.replace('Havard', 'Harvard')

# CHECK
len(adata.obs[adata.obs['brain_bank'] == 'Havard'])


# ### PLOT - using scanpy

# In[ ]:


# PLOT - umap by cell type
sc.pl.umap(adata, color='cell_type')


# In[ ]:


# FILTER - by celltype = T-cells
TC = adata[adata.obs['cell_type'] == 'TC']
TC


# In[ ]:


# PLOT - the T-cell filtered data
sc.pl.umap(
    TC, 
    color = 'cell_type')


# ### FILTERED UMAPS
# 
# List of variables to plot:
# ```
# - young VS old (binned)
# - 'age' on a continuous spectrum
# - 'pmi'
# - filter out PC??? --> ask Adam
# ```

# In[ ]:


# Get PMI range
adata.obs['pmi'].describe()


# In[ ]:


# PLOT - UMAP by PMI
sc.pl.umap(adata, color = 'pmi')


# In[ ]:


# Get age range
adata.obs['age'].describe()


# In[ ]:


# PLOT - UMAP by age
sc.pl.umap(adata, color='age')


# In[ ]:


# CALCULATE - Age stats
print('Describe Age column = \n', df.describe()['Age'], '\n')

age_max = df.describe()['Age']['max']
print('age_max =', age_max)

age_min = df.describe()['Age']['min']
print('age_min =', age_min)

# Math
print(f'{age_max} - {age_min} = {age_max - age_min}')
print(f'{age_max - age_min} / 2 = {(age_max - age_min) / 2}')
print(f'{age_min} + {(age_max - age_min) / 2} = {age_min + (age_max - age_min) / 2}')

# NOTE - ok to use the 50% row
age_half = df.describe()['Age']['50%']
print('age_half =', age_half)


# In[ ]:


# CHECK
adata.obs['age']


# In[ ]:


# EXTRACT - Age col, split into young/old, unique vals only, convert to a list, 
age_old = df[df['Age'] >= age_half]['Age'].unique().tolist()
print('age_old =', age_old)
print('\n')
age_young = df[df['Age'] < age_half]['Age'].unique().tolist()
print('age_young =', age_young)


# In[ ]:


# CREATE - new 'age_old' col, based on whether it's higher/lower than the age_half pt
adata.obs['age_old'] = adata.obs['age'] >= age_half
adata.obs['age_old']


# In[ ]:


# PLOT - UMAP of old VS young ages
sc.pl.umap(
    adata, 
    color = 'age_old'
    )


# In[ ]:


# EDIT - adata obj creating a new col called `age_bin' and the value depends on age being >= the halfway point
adata.obs['age_bin'] = np.where(
    adata.obs['age_old'], 
    f'old (>= {age_half})', 
    f'young (< {age_half})')
adata.obs['age_bin']


# In[ ]:


# PLOT - Customize UMAP of old VS young ages
fig = sc.pl.umap(
    adata, 
    color = 'age_bin', 
    title = 'Control group ages: \nYoung VS Old',
    return_fig = True
    )
ax = fig.axes[0]
ax.legend_.set_title('Age (yrs)')


# In[ ]:


# Define a function to check for existing dir's and create them if missing
def my_makedirs(path):
    """
    Provide a full file path, this will grab just the parent directory and check if it exists, if not, it will create it.
    """
    # Get the parent directory path
    path = os.path.dirname(path)
    
    # If-else logic to check if they exist, if not, create the directory
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f"Created: {path}")
    else:
        print(f"\'/{path}\' folder is already there")

# Define output file path for my testing, base on the original file name
output_file = os.path.join('output', os.path.basename(input_file).replace('.h5ad', '-ef.h5ad'))
print('output_file =', output_file)

# Call function to check for or create that path to the output file location
my_makedirs(output_file)


# ### QC

# In[ ]:


# Mitochondrial genes "MT-" for human
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith (('RPS', 'PRL'))


# In[ ]:


sc.pp.calculate_qc_metrics(
    adata,
    qc_vars = ['mt', 'ribo'],
    inplace = True,
    log1p = True
)


# In[ ]:


# PLOT - violin plot of QC metrics
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter = 0.4,
    multi_panel = True,
)


# Useful to consider QC metrics jointly by inspecting a scatter plot colored by `pct_counts_mt`

# In[ ]:


sc.pl.scatter(
    adata,
    'total_counts',
    'n_genes_by_counts',
    color = 'pct_counts_mt'
)


# # FILTER??? - Check snakemake rules 1st

# In[ ]:


# sc.pp.filter_cells(adata, min_genes = 100)

# #? Q) How to check this?


# In[ ]:


# # CHECK - starting # of genes
# print(f'# of genes pre-filtering: {adata.n_vars}')

# # Remove genes detected in < 3 cells (along w/the 0 count genes)
# sc.pp.filter_genes(adata, min_cells = 3)

# # CHECK - remaining genes
# print(f'# of genes post-filtering: {adata.n_vars}')


# ### DIFFERENTIALLY-EXPRESSED (DE) GENES AS MARKERS
# 
# Can also calc marker genes/cluster and then look up whether we can link those marker genes to known biology (ie. cell types and/or states). Can be done using simple stat tests (ie. Wilcoxon and t-test) for each cluster VS the rest.

# In[ ]:


# CHECK
adata


# In[ ]:


# # CHECK - already done by Adam
# sc.pp.log1p(adata)


# In[ ]:


#! WARNING - It seems you use rank_genes_groups on the raw count data. Please logarithmize your data before calling rank_genes_groups.
#* NOTE - need to specify the layer as 'data' and `use_raw = False`
#! WARNING - took 32 mins - ~3 hrs, might crash in sinteractive (run as sbatch.sh script then reload the adata obj afterwards)

# Obtain cluster-specific DE genes
sc.tl.rank_genes_groups(
    adata, 
    groupby = 'age_bin',
    method  = 'wilcoxon',
    layer   = 'data' ,      # ADAM: already normalized + logarithmized
    use_raw = False,
)


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # SAVE + RELOAD

# In[ ]:


# Define output file path for my testing, base on the original file name
output_file = os.path.join('output', os.path.basename(input_file).replace('.h5ad', '-ef-rank_genes_groups.h5ad'))
print('output_file =', output_file)


# In[ ]:


# # Save the AnnData object
# adata.write_h5ad(
#     filename = output_file, 
#     compression = 'gzip'
#     )


# In[ ]:


# RELOAD - my saved AnnData object
adata = sc.read_h5ad(output_file)
adata


# In[ ]:


import pprint

pprint.pprint(adata.uns)


# In[ ]:


# EXTRACT - list of ranked genes (I think it's returning highest expression in BOTH old/young groups)
rank_genes_groups_names_array = adata.uns['rank_genes_groups']['names']
print('rank_genes_groups_names_array.shape =', rank_genes_groups_names_array.shape)
print('rank_genes_groups_names_array[:30] =', rank_genes_groups_names_array[:30])


# In[ ]:


# CHECK - see the DF of the rank genes groups 
sc.get.rank_genes_groups_df(adata, group = None)


# In[ ]:


adata.var


# In[ ]:


# KEEP!!!
# IMPORTS
import matplotlib.colors as mcolors

# Log scaling of the palette
norm = mcolors.LogNorm()


# In[ ]:


# # TEST - EXAMPLE - scaling the vmin/vmax

# # # Make mock col w/log-normally distributed values
# # adata.obs['lognormal'] = np.random.lognormal(3, 1, adata.shape[0])

# # PLOT - heatmap
# sc.pl.umap(
#     adata,
#     color = 'lognormal', 
#     s = 20,
#     norm = norm,
# )

# # #? Q) What does this part do?
# # #? A) Cleans up the 'lognormal' col it created and drops it
# # adata.obs.drop(
# #     'lognormal', 
# #     axis = 1,
# #     inplace = True
# #     )

# #* TODO - Vmin/Vmax - normalize axes based on the range in my dataset
# # - by separately plotting these to see


# In[ ]:


#* TODO - find where the cols in rank_genes_groups comes frmo


# ### #? Q) What is the score?
# 
# #? A) Wilcoxon test results

# In[ ]:


# FILTER - young ranked gene groups
dc_cluster_genes_young = sc.get.rank_genes_groups_df(
    adata, 
    group = 'young (< 57.9)'
).head(5)['names']

# PLOT - heatmap of top 5 DE genes from the young group
sc.pl.umap(
    adata,
    color       = dc_cluster_genes_young,
    legend_loc  = 'on data',
    frameon     = True,
    ncols       = 3
)

#! WARNING - nearly all 0s

#* TODO - Vmin/Vmax - normalize axes based on the range in my dataset
# - by separately plotting these to see


# In[ ]:


# FILTER - young group, just names
dc_cluster_genes_young = sc.get.rank_genes_groups_df(
    adata, 
    group = 'young (< 57.9)'
).head(5)['names']

dc_cluster_genes_young


# In[ ]:


# Normalizing legend

# PLOT - heatmap of top 5 DE genes from the young group
sc.pl.umap(
    adata,
    color       = dc_cluster_genes_young,
    s           = 20,       # Size of points
    norm        = norm,     # Log normalization of color palette
    legend_loc  = 'on data',
    frameon     = True,
    ncols       = 3
)


# In[ ]:


# FILTER - old group, just names
dc_cluster_genes_old = sc.get.rank_genes_groups_df(
    adata, 
    group = 'old (>= 57.9)'
).head(5)['names']

dc_cluster_genes_old


# In[ ]:


# Normalizing legend

# PLOT - heatmap of top 5 DE genes from the old group
sc.pl.umap(
    adata,
    color       = dc_cluster_genes_old,
    s           = 20,       # Size of points
    norm        = norm,     # Log normalization of color palette
    legend_loc  = 'on data',
    frameon     = True,
    ncols       = 3,
)


# In[ ]:


# TEST - scaling the vmin/vmax using quartiles?

# IMPORTS
from functools import partial

# PLOT - heatmap of top 5 DE genes from the young group
sc.pl.umap(
    adata,
    # color       = [dc_cluster_genes_young, dc_cluster_genes_old],     #! ERROR
    color       = dc_cluster_genes_young,
    legend_loc  = 'on data',
    frameon     = True,
    ncols       = 5,
    vmin        = 'p90',
    vmax        = 'p90',
)

#* TODO - Vmin/Vmax - normalize axes based on the range in my dataset
# - by separately plotting these to see


# In[ ]:


# # %pip install monkeybread

# # Remove existing installation
# %pip uninstall monkeybread

# # Install fresh copy
# %pip install --upgrade git+https://github.com/immunomind/monkeybread.git


# In[ ]:


import monkeybread as mb 

# PLOT - volcano plot
mb.plot.volcano_plot(
    adata,
    groups = None,
)


# #* TODO - Xylena: pseudobulk each category (young/old) then run w/the tools Adam used
# - avg's expression of each cell type in each sample
# - should have a tool in scanpy
# - also in Adam's pipeline
# 
# ADAM
# - break up the adata obj into indiv clusters by cell type for the trajectory analysis
# - redo the DGE by cell type, 
# - use the leiden_2 cluster

# #* TODO - where is the list of genes from ranking the genes coming from?

# In[ ]:





# #* TODO - pseudo-bulking: collapsing into counts the mean expression for that sample
# - sum is better for genes that are lowly expressed by highly variable betw samples
# - means is better for genes that are highly expressed genes
# 

# In[ ]:


# IMPORTS
import matplotlib.pyplot as plt

# Do DEA (alreadt done!)
# - groupby = 'cell_type'

# PLOT - volcano plot
plt.figure(figsize = (10, 8))
sc.pl.volcano_plot(
    adata,
    groupby = 'age_bin',
    min_pctl = 0.1,
    logfc_threshold = 0.25,
    pval_threshold = 0.05,
    show = True,
)


# In[ ]:


# Lit review on pseudotime and trajectory analysis tools, run by some suggestions by Xylena before starting on one

