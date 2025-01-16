import muon as mu
import scanpy as sc


mdata = mu.read(snakemake.input.anndata)


muon.pp.neighbors(mdata, key_added='wnn', metric='cosine', low_memory=True)

sc.tl.leiden(mdata, neighbors_key='wnn', resolution=0.3)

mdata.mod['rna'].obs['leiden'] = mdata.mod['atac'].obs['leiden'] = mdata.obs['leiden']

mu.tl.umap(mdata, neighbors_key='wnn', min_dist=0.3)

mdata.write_h5mu(filename=snakemake.output.muon_object, compression='gzip') # type: ignore