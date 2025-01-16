import scanpy as sc
import pandas as pd
import numpy as np
import snapatac2 as snap

# Open the data containing the important diagnosis samples
important_diagnosis_df = pd.read_csv('/data/CARD_singlecell/SN_atlas/input/SN_diagnoses2024 - control50_PD_DLB_only.csv')
important_diagnosis_df = important_diagnosis_df[['SampleID', 'PrimaryDiagnosis']].iloc[:-3]
important_diagnosis_df['PrimaryDiagnosis'].drop_duplicates()

# Change column names for easier merge
important_diagnosis_df.columns = ['sample', 'primary diagnosis']

# Import sample metadata
samples = pd.read_csv('/data/CARD_singlecell/SN_atlas/input/SNsamples.csv')
samples = samples.replace({np.NaN: 'None'})

# From the metadata csv extract the ID of the sample and the batch folder it is saved under
sample_names = samples['Sample_ID'].to_list()

# Remove the dataset
important_diagnosis_df['sample'] = [str(x) for x in important_diagnosis_df['sample']]
important_diagnosis_df = important_diagnosis_df[important_diagnosis_df['sample'].isin(sample_names)]
samples_to_keep = important_diagnosis_df['sample'].to_numpy()
important_diagnosis_df['sample']

# Read list of atac data locations
atac_anndata = ['/data/CARD_singlecell/SN_atlas/data/samples/' + x + '/03_' + x + '_anndata_filtered_atac.h5ad' for x in samples_to_keep]

# Read in snapATAC2 datasets into a list of anndata objects in read only
list_of_anndata = [(samples_to_keep[i], snap.read(atac_anndata[i])) for i in range(len(atac_anndata))]

anndataset = snap.AnnDataSet(
    adatas=list_of_anndata,
    filename='/data/CARD_singlecell/SN_atlas/data/atlas/02_filtered_anndata_atac.h5ad'
)

# Select variable features
snap.pp.select_features(anndataset, n_features=250000, n_jobs=60)

# Spectral MDS analysis
snap.tl.spectral(anndataset)

# Batch correction›
snap.pp.mnc_correct(anndataset, batch="sample", key_added='X_spectral')

# Perform k-nearest neighbors
snap.pp.knn(anndataset)

# Cluster 
snap.tl.leiden(anndataset)

# Calculate umap
snap.tl.umap(anndataset)

# Write out the umap coordinates
umap_df = anndataset.obsm['X_umap']
umap_df = pd.DataFrame(umap_df)
umap_df.to_csv(snake.output.atac_umap)

# Write out the selected bins and count of each bin
var_df = pd.merge(
    left=pd.DataFrame(anndataset.var['count']),
    right=pd.DataFrame(anndataset.var['selected']),
    left_index=True,
    right_index=True
)
var_df.columns= ['count', 'selected']
var_df.to_csv(snake.output.atac_var)



# Be kind, rewind
anndataset.close()