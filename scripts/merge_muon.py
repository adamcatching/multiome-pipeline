import muon as mu
import scanpy as sc

mdata = mu.MuData({"RNA": sc.read_h5ad(snakemake.input.merged_rna_anndata), "ATAC": sc.read_h5ad(snakemake.input.merged_atac_anndata)})

mdata.write_h5mu(snakefile.output.merged_multiome)