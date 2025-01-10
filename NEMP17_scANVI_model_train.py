#! /usr/bin/env python3
#Load required packages and modules.
import os
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
import random
from scvi.model.utils import mde

#Define variables and settings.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/refmap_NEMP17'
os.chdir(directory_path)
save_name = 'NEMP17'

#Load the dataset.
print('Loading dataset...')
concat_preprocessed = directory_path + '/ref_query_concat_preprocessed.h5ad'
adata_concat = sc.read_h5ad(concat_preprocessed)

#Load the scVI model.
print('Loading scVI model...')
scvi_model = scvi.model.SCVI.load(directory_path + '/concat_scvi_model/', adata_concat)

#Set up labels to transfer.
print('Setting up dataset...')
adata_concat.obs["labels_scANVI"] = 'Unknown'
ref_mask = adata_concat.obs["dataset"] == "reference"
adata_concat.obs["labels_scANVI"][ref_mask] = adata_concat.obs.age_bracket_by_region_by_cell_class[ref_mask].values

#Set up scANVI model.
print('Setting up scANVI model...')
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, adata=adata_concat, unlabeled_category="Unknown", labels_key='labels_scANVI')

#Training scANVI model.
print('Training scANVI model...')
scanvi_model.train(max_epochs=200, n_samples_per_label=100)

#Save scANVI model.
print('Saving scANVI model...')
scanvi_model.save(directory_path + '/concat_scanvi_model/', overwrite=True)

#Save dataset.
print('Saving dataset...')
concat_preprocessed = directory_path + '/ref_query_concat_preprocessed.h5ad'
adata_concat.write(concat_preprocessed, compression='gzip')



