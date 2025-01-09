#! /usr/bin/env python3
#Import required packages and modules.
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

#Define variables and settings.
directory_path = '/wynton/home/pollenlab/reedmcmullen/projects/NEMP17/scanpy_NEMP17'
os.chdir(directory_path)
save_name = "NEMP17"

#Load the dataset.
print('Loading the dataset...')
anndata_preprocessed = directory_path + '/NEMP17_preprocessed.h5ad'
adata = sc.read_h5ad(anndata_preprocessed)

#Add metadata.
print('Adding metadata...')
adata.obs['time_point'] = adata.obs['GEMwell'].apply(lambda x: 'day 2' if 1 <= x <= 8 else 'day 7' if 9 <= x <= 16 else None)
adata.obs['pool_tx'] = adata.obs['GEMwell'].apply(lambda x: 'intraspecies' if x in [1,3,5,7,9,11,13,15] else 'interspecies' if x in [2,4,6,8,10,12,14,16] else None)
adata.obs['fgf2_tx'] = adata.obs['GEMwell'].apply(lambda x: '0 ng/mL' if x in [1,2,5,6,9,10,13,14] else '20 ng/mL' if x in [3,4,7,8,11,12,15,16] else None)
adata.obs['human_indiv_set'] = adata.obs['GEMwell'].apply(lambda x: 'set 1' if x in [1,2,3,4,9,10,11,12] else 'set 2' if x in [5,6,7,8,13,14,15,16] else None)

#Save the dataset.
print('Saving the dataset...')
anndata_preprocessed = directory_path + '/NEMP17_preprocessed.h5ad'
adata.write(anndata_preprocessed, compression='gzip')

#Plotting metadata
print('Plotting metadata...')

#UMAP plot.
plt.rcParams["figure.figsize"] = (6, 6)
sc.pl.umap(adata, color=['time_point', 'pool_tx', 'fgf2_tx', 'human_indiv_set'], save=f'_{save_name}_fixed_effects.png')

#Barplots.
adata.obs['time_point'] = pd.Categorical(adata.obs['time_point'], categories=['day 2', 'day 7'], ordered=True)
data = adata.obs['time_point'].value_counts(normalize=True).reindex(['day 2', 'day 7'])
data = round(data, 3) * 100
colors = adata.uns['time_point_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Time Point')
plt.title('Time Point Composition')
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_timepoint_barchart_initial.png'), bbox_inches='tight')

adata.obs['pool_tx'] = pd.Categorical(adata.obs['pool_tx'], categories=['intraspecies', 'interspecies'], ordered=True)
data = adata.obs['pool_tx'].value_counts(normalize=True).reindex(['intraspecies', 'interspecies'])
data = round(data, 3) * 100
colors = adata.uns['pool_tx_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('Pool Treatment')
plt.title('Pool Treatment Composition')
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_pooltx_barchart_initial.png'), bbox_inches='tight')

adata.obs['fgf2_tx'] = pd.Categorical(adata.obs['fgf2_tx'], categories=['0 ng/mL', '20 ng/mL'], ordered=True)
data = adata.obs['fgf2_tx'].value_counts(normalize=True).reindex(['0 ng/mL', '20 ng/mL'])
data = round(data, 3) * 100
colors = adata.uns['fgf2_tx_colors']
ax = data.plot.bar(color=colors)
for container in ax.containers:
    ax.bar_label(container)
ax.grid(False)
plt.ylabel('Percent of Total Cells')
plt.xlabel('FGF2 Treatment')
plt.title('FGF2 Treatment Composition')
plt.savefig(os.path.join(directory_path, 'figures', f'{save_name}_fgf2tx_barchart_initial.png'), bbox_inches='tight')



