#Load required packages and modules.
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch

#Define working directory and settings.
sc.set_figure_params(figsize=(6, 6), frameon=True)
torch.set_float32_matmul_precision("high")
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17'
os.chdir(directory_path)
save_name='NEMP17'

#Load dataset.
print('Loading dataset...')
concat_preprocessed = directory_path + '/ref_query_concat_preprocessed.h5ad'
adata_concat = sc.read_h5ad(concat_preprocessed)

#Set up dataset for scVI model.
scvi.model.SCVI.setup_anndata(adata_concat, layer="counts", batch_key="sample_id")

#Set up scVI model.
arches_params = dict(use_layer_norm="both", use_batch_norm="none", encode_covariates=True, dropout_rate=0.2, n_layers=2)
scvi_model = scvi.model.SCVI(adata_concat, **arches_params)

#Train scVI model.
scvi_model.train(early_stopping=True, max_epochs=250, train_size=0.75, check_val_every_n_epoch=5)

#Save the scVI model.
scvi_model.save(directory_path + '/concat_scvi_model/', overwrite=True)
