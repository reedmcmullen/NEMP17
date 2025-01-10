#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

#Define variables and settings.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP25/scanpy_NEMP17'
os.chdir(directory_path)
save_name = 'NEMP17'

#Load dataset.
print('Loading the dataset...')
anndata_initial = directory_path + f'{save_name}_initial.h5ad'
adata = sc.read_h5ad(anndata_initial)

#Define QC cutoffs.
print('Filtering the dataset...')
n_genes_cutoff = 9000
total_counts_cutoff = 30000
mito_cutoff = 10
ribo_cutoff = 5

#Plot violin plots of QC metrics with cutoffs.
fig, axes = plt.subplots(1, 4, figsize=(16, 6))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], stripplot=False, show=False)
axes[0].axhline(n_genes_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'total_counts', ax=axes[1], stripplot=False, show=False)
axes[1].axhline(total_counts_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_mito', ax=axes[2], stripplot=False, show=False)
axes[2].axhline(mito_cutoff, color='red', linestyle='--')
sc.pl.violin(adata, 'pct_counts_ribo', ax=axes[3], stripplot=False, show=False)
axes[3].axhline(ribo_cutoff, color='red', linestyle='--')
plt.savefig(os.path.join(directory_path, 'figures', f'_{save_name}_qc_metrics_cutoffs.png'), bbox_inches='tight')

#Filter outlier cells from each AnnData object based on QC metric plots.
adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
adata = adata[adata.obs.total_counts < total_counts_cutoff, :]
adata = adata[adata.obs.pct_counts_mito < mito_cutoff, :]
adata = adata[adata.obs.pct_counts_ribo < ribo_cutoff, :]

#Preprocess the dataset.
print('Preprocessing the dataset...')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="GEMwell", flavor='seurat')
sc.pl.highly_variable_genes(adata, save=f'{save_name}.png')
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata)
res = 1
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=res)

#Save the results.
print('Saving preprocessed AnnData object')
anndata_preprocessed = directory_path + f'{save_name}_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')
