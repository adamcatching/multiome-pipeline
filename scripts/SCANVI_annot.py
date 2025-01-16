import numpy as np
import scanpy as sc
import scvi

# Open the RNA merged and filtered
adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

# Make a copy of the AnnData atlas that only contains variable genes
filtered_adata = adata[:, adata.var['highly_variable']].copy()

model = scvi.model.SCVI.load(snakemake.input.model, adata = filtered_adata)

SCANVI_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=filtered_adata,
    labels_key="cell_type",
    unlabeled_category='Unknown'
    )


SCANVI_model.train(
    accelerator='gpu', 
    max_epochs=1000,
    early_stopping=True,
    early_stopping_patience=20,
    devices=1
    )

SCANVI_LATENT_KEY = "X_scANVI"
SCANVI_PREDICTIONS_KEY = "C_scANVI"

elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']
elbo.to_csv(snakemake.output.model_history, index=False)#snakemake.output.model, index=False)


adata.obsm[SCANVI_LATENT_KEY] = SCANVI_model.get_latent_representation(filtered_adata)
adata.obs[SCANVI_PREDICTIONS_KEY] = SCANVI_model.predict(adata)

model.save('data/CARD_singlecell/SN_atlas/wnn_pipeline/models/rna', overwrite=True) 
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.umap(adata, min_dist=0.3)

sc.tl.leiden(adata, flavor='igraph', resolution=1, key_added='leiden_SCANVI')

adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')
