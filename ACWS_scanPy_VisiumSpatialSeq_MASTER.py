#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 14:31:34 2020

This is the main script for processing scRNA-Seq Data through ScanPy.

@author: smith
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from skimage import io
import matplotlib
import openpyxl
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama as sca
import countDEGs as cd
import ACWS_filterCells as fcx

new = True #if new analysis, set to True read adata from 10X mtx or cache). If re-analyzing data, set to false and read from results_file path.

###SET DIRECTORY TO READ/WRITE DATA.
#THIS SHOULD BE THE DIRECTORY CONTAINING THE .MTX DATA FILE AND .TSV BARCODES & FEATURE FILES:
BaseDirectory = '/d1/studies/scanPy/spatial/March2021_FinalAnalysis/'
sampleName = 'Spatial_Stressed_Cropped' #This is used for name result output files
batches = False #Set to True if you need to do batch correction (i.e. if samples were taken to core and sequenced at different times)
import os
os.chdir(BaseDirectory)
%logstart -o scanpy_log.txt

###SET SCANPY SETTINGS:
results_file = os.path.join(BaseDirectory, sampleName + '_scanpy_results.h5ad')  # the file that will store the analysis results
results_file_partial = os.path.join(BaseDirectory, sampleName + '_scanpy_adata_prefiltering.h5ad')

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.n_jobs=8 #use parallel processing when possible
sc.logging.print_header()
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')
matplotlib.rcParams.update({'text.usetex': False, 'font.family': 'stixgeneral', 'mathtext.fontset': 'stix',})

###LOAD DATA
adatac1 = sc.read_visium('/d1/studies/scanPy/spatial/slide_1_NDB_ctrl_1/',
                          count_file='/d1/studies/scanPy/spatial/slide_1_NDB_ctrl_1/slide_1_NDB_1_ctrl_ST_filtered_feature_bc_matrix.h5',
                          library_id='ctrl1')
adatac1.var_names_make_unique()
adatac1.obs['sample']='Ctrl1'
adatac2 = sc.read_visium('/d1/studies/scanPy/spatial/slide_1_NDB_ctrl_2/',
                         count_file='/d1/studies/scanPy/spatial/slide_1_NDB_ctrl_2/slide_1_NDB_2_ctrl_ST_filtered_feature_bc_matrix.h5',
                         library_id='ctrl2')
adatac2.var_names_make_unique()
adatac2.obs['sample']='Ctrl2'

adatas1 = sc.read_visium('/d1/studies/scanPy/spatial/slide_1_NDB_stress_1/',
                         count_file='/d1/studies/scanPy/spatial/slide_1_NDB_stress_1/Slide_1_NDB_stress_1_ST_filtered_feature_bc_matrix.h5',
                         library_id='stress1')
adatas1.var_names_make_unique()
adatas1.obs['sample']='Stress1'

adatas2 = sc.read_visium('/d1/studies/scanPy/spatial/slide_1_NDB_stress_2/',
                         count_file='/d1/studies/scanPy/spatial/slide_1_NDB_stress_2/slide_1_NDB_stress_2_ST_filtered_feature_bc_matrix.h5',
                         library_id='stress2')
adatas2.var_names_make_unique()
adatas2.obs['sample']='Stress2'

vmin=2000
vmax=25000

left=2000
right=3050
c1spx = adatac1.obsm['spatial'][:,1]
adatac1_sub = adatac1[(adatac1.obsm['spatial'][:,0]>left) & (adatac1.obsm['spatial'][:,0]<right)]
adatac1_sub.var["mt"] = adatac1_sub.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adatac1_sub, qc_vars=["mt"], inplace=True)
sc.pl.spatial(adatac1_sub, vmin=vmin, vmax=vmax, color=['total_counts'], color_map=color_map, title='Control1 Total Counts', save='_control1_' + str(left) + '-' + str(right) + '_norm')

left=1900
right=2950
adatac2_sub = adatac2[(adatac2.obsm['spatial'][:,0]>left) & (adatac2.obsm['spatial'][:,0]<right)]
adatac2_sub.var["mt"] = adatac2_sub.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adatac2_sub, qc_vars=["mt"], inplace=True)
sc.pl.spatial(adatac2_sub, vmin=vmin, vmax=vmax, color=['total_counts'], color_map=color_map, title='Control2 Total Counts', save='_control2_' + str(left) + '-' + str(right) + '_norm')

left = 1900
right=2950
adatas1_sub = adatas1[(adatas1.obsm['spatial'][:,0]>left) & (adatas1.obsm['spatial'][:,0]<right)]
adatas1_sub.var["mt"] = adatas1_sub.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adatas1_sub, qc_vars=["mt"], inplace=True)
sc.pl.spatial(adatas1_sub, vmin=vmin, vmax=vmax, color=['total_counts'], color_map=color_map, title='Stress1 Total Counts', save='_stress1_' + str(left) + '-' + str(right) + '_norm')

left=2000
right=3050
adatas2_sub = adatas2[(adatas2.obsm['spatial'][:,0]>left) & (adatas2.obsm['spatial'][:,0]<right)]
adatas2_sub.var["mt"] = adatas2_sub.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adatas2_sub, qc_vars=["mt"], inplace=True)
sc.pl.spatial(adatas2_sub, vmin=vmin, vmax=vmax, color=['total_counts'], color_map=color_map, title='Stress2 Total Counts', save='_Stress2_' + str(left) + '-' + str(right) + '_norm')

###MERGE & CONCATENATE ADATAS:
adatas = [adatac1_sub, adatac2_sub, adatas1_sub, adatas2_sub]

spadata = adatac1_sub.concatenate(adatac2_sub,adatas1_sub,adatas2_sub, batch_key='library_id', uns_merge='unique')
samples = spadata.obs['sample'].tolist()
conditions=[]
for s in samples:
    if s.startswith('Ctrl'):
        conditions.append(0)
    elif s.startswith('Stress'):
        conditions.append(1)
spadata.obs['condition']=conditions
spadata.obs['sample']=spadata.obs['sample'].astype('category')
spadata.obs['condition']=spadata.obs['condition'].astype('category')

g1n='Control'
g2n='Stress'

spadata.write('SpAdata_ConcatenatedOnly.h5ad')
###BASIC QC:
spadata.var["mt"] = spadata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(spadata, qc_vars=["mt"], inplace=True)

sc.pl.violin(spadata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_' + str(sampleName) + '_plot_percentMito_raw')
sc.pl.scatter(spadata, x='total_counts', y='pct_counts_mt', save='_' + str(sampleName) + '_mito_counts_raw')
sc.pl.scatter(spadata, x='total_counts', y='n_genes_by_counts', save='_' + str(sampleName) + '_genes_counts_raw')


###SAVE CONCATENATED OBJECT:
spadata.write('SpAdata_Concatenated_Prefiltering.h5ad')

###IF YOU NEED TO LOAD UNCOMMENT THIS OUT:
spadata = sc.read('SpAdata_Concatenated_Prefiltering.h5ad')
###BASIC FILTERING OF SPOTS BASED ON TOTAL COUNTS & EXPRESSED GENES:
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(spadata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(spadata.obs["total_counts"][spadata.obs["total_counts"] < 20000], kde=False, bins=40, ax=axs[1])
sns.distplot(spadata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(spadata.obs["n_genes_by_counts"][spadata.obs["n_genes_by_counts"] < 6000], kde=False, bins=60, ax=axs[3])
fig.savefig('figures/Covariants_Histogram.png')


#sc.pp.filter_cells(spadata, min_counts=2000)
#sc.pp.filter_cells(spadata, max_counts=60000)
spadata = spadata[spadata.obs["pct_counts_mt"] < 15].copy()
#spadata = spadata[spadata.obs['pct_counts_ribo'] < 7].copy()
print(f"Cells after MT filter: {spadata.n_obs}")

#sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

###PLOT % MITO & TOTAL COUNTS:
sc.pl.violin(spadata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_' + str(sampleName) + '_plot_percentMito_filtered')
sc.pl.scatter(spadata, x='n_genes_by_counts', y='pct_counts_mt', save='_' + str(sampleName) + '_mito_counts_filtered')
sc.pl.scatter(spadata, x='total_counts', y='n_genes_by_counts', save='_' + str(sampleName) + '_genes_counts_filtered')

sc.pp.regress_out(['total_counts'])
###THE BELOW CALCULATIONS MAY HELP YOU DETERMINE CRITERIA FOR FILTERING BASED ON MAX GENES EXPRESSED and % mt-DNA.
###THE DEFAULTS HERE WORK WELL FOR MY DATA, BUT YOU ARE HIGHLY ENCOURAGED TO TEST DIFFERENT NUMBERS.

#Max genes (1SD over mean):
n_genes_mean = np.mean(spadata.obs.n_genes_by_counts)
n_genes_std = np.std(spadata.obs.n_genes_by_counts)
n_genes_suggestedMax = n_genes_mean + n_genes_std
print("Suggested max genes = " + str(n_genes_suggestedMax))

#Max genes by percentile:
max_genes_percentile=95
max_genes_perc = np.percentile(spadata.obs.n_genes_by_counts, max_genes_percentile)
print(str(max_genes_percentile) + "% of cells express fewer than " + str(round(max_genes_perc,4)) + " genes")

#Max mito (1SD over mean):
mito_mean = np.mean(spadata.obs.pct_counts_mt)
mito_std = np.std(spadata.obs.pct_counts_mt)
mito_suggestedMax = mito_mean + mito_std
print("Suggested max Mito% = " + str(mito_suggestedMax))

#Max mito% by percentile:
mito_percentile = 95
mito_perc = np.percentile(spadata.obs.pct_counts_mt, mito_percentile)
print(str(mito_percentile) + " percent of cells have less than " + str(round(mito_perc,4)) + " % mt-DNA")
###FILTER - TUNE THESE PARAMETERS BASED ON INFO FROM LAST STEP
max_genes = 5815
spadata = spadata[spadata.obs['n_genes_by_counts'] < max_genes, :].copy()
sc.pl.scatter(spadata, x='total_counts', y='pct_counts_mt', save='_' + str(sampleName) + '_filtered_mito_counts', title='Filtered < ' + str(max_genes) +  ' total counts')
sc.pp.filter_genes(spadata, min_cells=3)

spadata.raw = spadata

###SAVE THE CURRENT ADATA STATE, CAN REVERT TO THIS WHILE TUNING PARAMETERS:
spadata.write(results_file_partial)
###IF YOU LATER WANT TO REVERT TO THE PRE-FILTERED RESULTS FILE USE:
spadata = sc.read(results_file_partial)


###NORMALIZE & LOG TRANSFORM
sc.pp.normalize_total(spadata, target_sum=1e4, max_fraction=.1, exclude_highly_expressed=True)
sc.pp.log1p(spadata)


###PRE-PROCESS DATA, SELECT HIGHLY VARIABLE GENES. 
###YOU WILL LIKELY WANT TO PLAY WITH THESE VARIABLES AND SEE HOW IT AFFECTS RESULTS.
sc.pp.highly_variable_genes(spadata)
                            
min_mean = .02
max_mean = 3.5
min_disp = 0.3
if batches:
    sc.pp.highly_variable_genes(spadata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, batch_key='sample')
elif not batches:
    sc.pp.highly_variable_genes(spadata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp) 
"""This adds following dimensions to the data:
    'highly_variable', boolean vector (adata.var)
    'means', float vector (adata.var)
    'dispersions', float vector (adata.var)
    'dispersions_norm', float vector (adata.var)"""
    

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(spadata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(spadata.obs["total_counts"][spadata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(spadata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(spadata.obs["n_genes_by_counts"][spadata.obs["n_genes_by_counts"] < 6000], kde=False, bins=60, ax=axs[3])
fig.savefig('figures/Covariants_Histogram_Normalized.' + sc.settings.file_format_figs)

spadata.var.to_excel(os.path.join(BaseDirectory, 'AdataVar_prefiltering.xlsx'))

###PLOT SOME BASIC QC BEFORE BATCH CORRECTION
sc.pl.highest_expr_genes(spadata, n_top=50, save='_highest_expr_genes_pre-batchCorrection')

#LOG TRANSFORM EXPRESSION:  
sc.pp.combat(spadata, key='sample')
spadata.var.to_csv('Spadata_Var_Combat.csv')

print(np.sum(np.isnan(spadata.X)))
print(np.any(spadata.X.sum(axis=0) == 0))
print(np.any(spadata.X.sum(axis=1) == 0))
###PLOT DISPERSIONS AND MEANS TO HELP REFINE:
sc.pl.scatter(spadata, x='n_cells', y='dispersions', save='_' + str(sampleName) + '_dispersions', title='Dispersions')
sc.pl.scatter(spadata, x='n_cells', y='dispersions_norm', save='_' + str(sampleName) + '_dispersions_norm', title='Dispersions Normalized')
sc.pl.scatter(spadata, x='n_cells', y='means', save='_' + str(sampleName) + '_means', title='Means')

###WRITE TABLE OF RESULTS BEFORE FILTERING:
spadata.var.to_excel(os.path.join(BaseDirectory, 'Adata_var_raw_preFiltering.xlsx'))

###UPDATE ADATA.RAW TO CONTAIN FILTERED DATA BEFORE HVG SELECTION
spadata.raw = spadata

###ACTUALLY DO THE FILTERING
sc.pl.highly_variable_genes(spadata, save='_' + str(sampleName) + '_highlyVariableGenes')
spadata = spadata[:, spadata.var['highly_variable']].copy()
#sc.pp.regress_out(spadata, ['total_counts'], n_jobs=12)

#Batch effect correction:
if batches:
    sc.pp.combat(adata, key='batch')

###WRITE EXCEL FILE OF HIGHLY VARIABLE GENES:
spadata.var.to_excel(os.path.join(BaseDirectory, str(sampleName) + 'HighlyVariableGenes_minMean' + str(min_mean) + '_maxMean' + str(max_mean) + '_min_disp' + str(min_disp) + '.xlsx'))

###SCALE EACH GENE TO UNIT OF VARIANCE, CLIP VALUES EXCEEDING MAX VARIANCE:
sc.pp.scale(spadata, max_value=10)
###EXAMPLE HOW TO PLOT EXPRESSION DISTRIBUTION INDIVIDUAL GENES:
sc.pl.violin(adata, 'Oprm1', save='_' + str(sampleName) + '_Oprm1')


###CREATE LIST OF GENES TO LABEL ON PCA/UMAP PLOTS:
labeled_genes_var = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Pdgfra', 'Flt1']
labeled_genes = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Ntsr2', 'Pdgfra', 'Tmem119', 'Flt1', 'Pcbp3']

#labeled_genes = ['Oprm1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Ntsr2', 'Pdgfra', 'Gpr17', 'Tmem119', 'C1qc', 'Cldn5', 'Flt1', 'Dbi']


#RUN PCA
n_comps = 30
sc.tl.pca(spadata, n_comps=n_comps, svd_solver='arpack')
sc.pl.pca(spadata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_comps) + 'comps_PCA_labeled')
sc.pl.pca_variance_ratio(spadata, log=True, save='_' + str(n_comps) + '_ncomps_PCA_VarianceRatio')
###LOOK AT VARIANCE_RATIO OUTPUT FILE BEFORE PROCEEDING TO NEXT STEP
n_pcs = 7
sc.pl.pca_overview(spadata, save='_' + str(sampleName) + '_' + str(n_pcs) + 'PCs_PCA_Overview')
sc.pl.pca_loadings(spadata, components=list(np.arange(1, n_pcs+1)), save='_' + str(sampleName) + '_' + str(n_pcs) + '_PCs')

#COMPUTING NEIGHBORHOOD GRAPH:
#Uses PCA representation of data matrix
###LOOK AT SAVED FIGURE W/ SUFFIX _PCA_VarianceRatio AND CHOOSE NUMBER OF PCs BEFORE APEX (HERE ~20)
n_neighbors = 25
min_dist = .1
sc.pp.neighbors(spadata, n_neighbors=n_neighbors, n_pcs=n_pcs)

###UMAP EMBEDDING:
sc.tl.umap(spadata, min_dist=min_dist)
vmin=0 #minimum to scale umap plots to
vmax_raw=None #maximum to scale umap plots to
color_map='plasma' #see options for colormaps at https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
sc.pl.umap(spadata, color=labeled_genes, color_map=color_map, vmin=vmin, vmax=vmax_raw, save='_' + str(sampleName) + '_' + str(n_neighbors) + 'neighbors_' + str(n_pcs) + 'PCs_umap_' + str(min_dist) + 'minDist_' + str(color_map) + '_labeled_raw')
###THIS REQUIRES GENES IN LABELED_GENES TO BE IN HIGHLY VARIABLE LIST:
vmax_filt=None #maximum to scale filtered umap plots to
sc.pl.umap(spadata, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, save='_' + str(sampleName) + '_nNeighbors' + str(n_neighbors) + '_nPCs' + str(n_pcs) + '_umap_filtered')


###ASSIGN CLUSTERS WITH LEIDEN ALGORITHM:
resolution = 0.5
iterations = -1 #only change this if you know what you're doing
sc.tl.leiden(spadata, resolution=resolution, n_iterations=iterations)
labeled_genes.insert(0,'leiden')
sc.pl.umap(spadata, color=labeled_genes, color_map=color_map, vmax=vmax_raw, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden')
labeled_genes_var.insert(0,'leiden')
sc.pl.umap(spadata, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered')
sc.pl.umap(spadata, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_leiden_only')

###ASSIGN CLUSTERS WITH LOUVAIN ALGORITHM:
resolution = 0.5
sc.tl.louvain(spadata, resolution=resolution)
labeled_genes.insert(0,'louvain')
sc.pl.umap(spadata, color=labeled_genes, color_map=color_map, vmax=vmax_raw, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_louvain')
labeled_genes_var.insert(0,'louvain')
sc.pl.umap(spadata, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_louvain_filtered')
sc.pl.umap(spadata, color=['louvain'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_louvain_only')

###PLOT UMAP WITH COUNTS & GENES LABELED:
sc.pl.umap(spadata, color=["total_counts", "n_genes_by_counts", "leiden"], color_map=color_map, wspace=0.5, save='_counts_genes_resolution' + str(resolution))
sc.pl.umap(spadata, color=["total_counts", "n_genes_by_counts", "leiden"], use_raw=False, color_map=color_map, wspace=0.5, save='_counts_genes_resolution' + str(resolution) + '_filtered')
sc.pl.umap(spadata, color=["total_counts", "n_genes_by_counts", "library_id"], use_raw=True, color_map=color_map, wspace=0.5, save='_counts_by_libraryID')

leiden_colors = spadata.uns['leiden_colors']
if '#8c564b' in leiden_colors:
    pos = leiden_colors.index('#8c564b')
    leiden_colors.remove('#8c564b')
    leiden_colors.insert(pos, '#fff000')
    spadata.uns['leiden_colors']=leiden_colors

###PLOT SPATIAL RECONSTRUCTION:
size=1
alpha=0.75
vmin=0
vmax=5
#sc.pl.spatial(spadata, library_id='stress1', img_key="hires", color=["total_counts", "n_genes_by_counts"], size=size, color_map=color_map, save='_spatialReconstruction')
#sc.pl.spatial(spadata, library_id='stress1', img_key="hires", color=["leiden"], size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution))
#sc.pl.spatial(spadata, library_id='stress1', img_key="hires", color=["louvain"], size=size, color_map=color_map, save='_spatialReconstruction_louvainClusters_resolution' + str(resolution))
#sc.pl.spatial(spadata, library_id='stress1', img_key="hires", vmin=vmin, vmax=vmax, color=labeled_genes_var, use_raw=False, size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution) + '_size_' + str(size) + '_alpha_' + str(alpha) + '_labeledGenes_filtered_scaled')

stress1 = spadata[spadata.obs['sample']=='Stress1'].copy()
sc.pl.spatial(stress1, library_id='stress1', img_key="hires", color=["total_counts", "n_genes_by_counts", "leiden"], size=size, color_map=color_map, save='_spatialReconstruction_stress1_wCounts_Clusters')
sc.pl.spatial(stress1, library_id='stress1', img_key="hires", alpha=alpha, color=["leiden"], title='Stress1', size=size, color_map=color_map, save='_spatialReconstruction_stress1_alpha'+str(alpha))
sc.pl.spatial(stress1, library_id='stress1', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_stress1_labeledGenes_alpha'+str(alpha))

stress2= spadata[spadata.obs['sample']=='Stress2'].copy()
sc.pl.spatial(stress2, library_id='stress2', img_key="hires", color=["total_counts", "n_genes_by_counts", "leiden"], size=size, color_map=color_map, save='_spatialReconstruction_stress2_wCounts_Clusters')
sc.pl.spatial(stress2, library_id='stress2', img_key="hires", alpha=alpha, color=["leiden"], title='Stress2', size=size, color_map=color_map, save='_spatialReconstruction_stress2_alpha'+str(alpha))
sc.pl.spatial(stress2, library_id='stress2', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_stress2_labeledGenes_alpha'+str(alpha))

ctrl1= spadata[spadata.obs['sample']=='Ctrl1'].copy()
sc.pl.spatial(ctrl1, library_id='ctrl1', img_key="hires", color=["total_counts", "n_genes_by_counts", "leiden"], size=size, color_map=color_map, save='_spatialReconstruction_Ctrl1_wCounts_Clusters')
sc.pl.spatial(ctrl1, library_id='ctrl1', img_key="hires", alpha=alpha, color=["leiden"], title='Control1', size=size, color_map=color_map, save='_spatialReconstruction_ctrl1_alpha'+str(alpha))
sc.pl.spatial(ctrl1, library_id='ctrl1', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_ctrl1_labeledGenes_alpha'+str(alpha))

ctrl2= spadata[spadata.obs['sample']=='Ctrl2'].copy()
sc.pl.spatial(ctrl2, library_id='ctrl2', img_key="hires", color=["total_counts", "n_genes_by_counts", "leiden"], size=size, color_map=color_map, save='_spatialReconstruction_Ctrl2_wCounts_Clusters')
sc.pl.spatial(ctrl2, library_id='ctrl2', img_key="hires", color=["leiden"], title='Control2', size=size, color_map=color_map, save='_spatialReconstruction_Ctrl2_wCounts_alpha'+str(alpha))
sc.pl.spatial(ctrl2, library_id='ctrl2', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_ctrl2_labeledGenes_alpha'+str(alpha))



labeled_genes_var.insert(0,'leiden')
sc.pl.spatial(stress1, library_id='stress1', img_key="hires", wspace=0.1, vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_stress1_labeledGenesClus_alpha'+str(alpha))
sc.pl.spatial(stress2, library_id='stress2', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_stress2_labeledGenesClus_alpha'+str(alpha))
sc.pl.spatial(ctrl1, library_id='ctrl1', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_ctrl1_labeledGenesClus_alpha'+str(alpha))
sc.pl.spatial(ctrl2, library_id='ctrl2', img_key="hires", vmax=vmax, alpha=alpha, color=labeled_genes_var, size=size, color_map=color_map, save='_spatialReconstruction_ctrl2_labeledGenesClus_alpha'+str(alpha))


###ISOLATE ONE CLUSTER & RE-CLUSTER
labeled_genes_var = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Pdgfra', 'Tmem119', 'Flt1', 'Dbi']
labeled_genes = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Ntsr2', 'Pdgfra', 'Tmem119', 'Flt1', 'Pcbp3', 'Dbi']


markers = ['Slc17a7', 'Tac1', 'Zic1', 'Nptxr', 'Lhx8', 'Nrgn', 'Cldn11', 'Fam131c', 'Hmgb2', 'Slc6a13']

sc.pl.spatial(stress1, library_id='stress1', img_key="hires", vmax=vmax, alpha=alpha, color=markers, size=size, color_map=color_map, save='_spatialReconstruction_stress1_MarkerGenes')
sc.pl.spatial(stress2, library_id='stress2', img_key="hires", vmax=vmax, alpha=alpha, color=markers, size=size, color_map=color_map, save='_spatialReconstruction_stress2_MarkerGenes')
sc.pl.spatial(ctrl1, library_id='ctrl1', img_key="hires", vmax=vmax, alpha=alpha, color=markers, size=size, color_map=color_map, save='_spatialReconstruction_ctrl1_labeledGenes_MarkerGenes')
sc.pl.spatial(ctrl2, library_id='ctrl2', img_key="hires", vmax=vmax, alpha=alpha, color=markers, size=size, color_map=color_map, save='_spatialReconstruction_ctrl2_labeledGenes_MarkerGenes')




cluster = 5
adata_clu = spadata[spadata.obs.leiden==str(cluster)].copy()
resolution=1.0
sc.tl.leiden(adata_clu, resolution=resolution, n_iterations=iterations)
labeled_genes.insert(0,'leiden')
sc.pl.umap(adata_clu, color=labeled_genes, color_map=color_map, vmax=vmax_raw, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_cluster' + str(cluster))
labeled_genes_var.insert(0,'leiden')
sc.pl.umap(adata_clu, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered_cluster' + str(cluster))
sc.pl.umap(adata_clu, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_leiden_only')

libraries = ['ctrl1', 'ctrl2', 'stress1', 'stress2']
for library in libraries:
    sc.pl.spatial(adata_clu, img_key="hires", library_id=library, color=["leiden"], size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution) + '_cluster' + str(cluster) + '_' + library)
    sc.pl.spatial(adata_clu, img_key="hires", library_id=library, color=labeled_genes_var, use_raw=False, size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution) + 'labeledGenesVar_cluster' + str(cluster) + '_' + library)
    sc.pl.spatial(adata_clu, img_key="hires", library_id=library, color=labeled_genes, size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution) + 'labeledGenes_cluster' + str(cluster) + '_' + library)

###RE-CLUSTER WITHOUT ISOLATING
sc.tl.leiden(adata, restrict_to=('leiden', [str(cluster)]), resolution=resolution, n_iterations=iterations)
labeled_genes_var.insert(0, 'leiden_R')
sc.pl.umap(adata, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered_cluster' + str(cluster) + '_full')
sc.pl.umap(adata, color=['leiden_R'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_leiden_only_reclustered_' + str(cluster))
sc.pl.spatial(adata, img_key="hires", alpha=0.5, color=["leiden_R"], size=size, color_map=color_map, save='_spatialReconstruction_leidenClusters_resolution' + str(resolution) + '_cluster' + str(cluster) + '_full_alpha0.5')

###SEPARATE LOUVAIN CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs = list(zip(spadata.obs['condition'], spadata.obs['louvain'].astype('int')))
spadata.obs['pairs'] = pairs
spadata.obs['louvain'] = spadata.obs['louvain'].values.remove_unused_categories()

###SEPARATE LEIDEN CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs_leid = list(zip(spadata.obs['condition'], spadata.obs['leiden'].astype('int')))
spadata.obs['pairs_leiden'] = pairs_leid
spadata.obs['leiden'] = spadata.obs['leiden'].values.remove_unused_categories()


#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts = spadata.obs['pairs'].value_counts().sort_index()
print(counts)

#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts_leid = spadata.obs['pairs_leiden'].value_counts().sort_index()
print(counts_leid)

###SPLIT DATA BY GROUP TO EXAMINE CLUSTERS FOR EACH GROUP INDIVIDUALLY:
g1n='Control'
g2n='Stress'
spadata_g1 = spadata[spadata.obs['condition']==0].copy()
spadata_g2 = spadata[spadata.obs['condition']==1].copy()
sc.pl.umap(spadata_g1, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_filtered_' + str(g1n))
sc.pl.umap(spadata_g2, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_filtered_' + str(g2n))
sc.pl.umap(spadata_g1, color=['leiden'], title=g1n, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_only_' + str(g1n))
sc.pl.umap(spadata_g2, color=['leiden'], title=g2n, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_only' + str(g2n))

###PLOT NUMBEROF CELLS IN EACH CLUSTER:
counts_g1 = counts_leid[:int(len(counts_leid)/2)]
counts_g2 = counts_leid[int(len(counts_leid)/2):]
cg1df = pd.DataFrame(counts_g1)
cg1df.reset_index(drop=True, inplace=True)
cg2df = pd.DataFrame(counts_g2)
cg2df.reset_index(drop=True, inplace=True)
cat = pd.concat([cg1df,cg2df],axis=1, ignore_index=True)
cat.columns=[g1n,g2n]
cf = cat.plot.bar(grid=False)
cf.set_ylabel('# Cells')
cf.set_xlabel('Leiden Cluster')
cf.set_title('Number of cells per cluster')
fig = cf.get_figure()
fig.savefig('figures/CellsPercluster.' + sc.settings.file_format_figs)

###RUN STATS, NOTE 't-test' can be changed to 'wilcoxon':
###COMPARE EXPRESSION BY CONDITION (EQUIVALENT TO BULK-SEQ):
spadata.obs.condition=spadata.obs.condition.astype('category')
sc.settings.verbosity = 2
method = 't-test' #t-test, t-test_overestim_var, wilcoxon, or logreg

sc.tl.rank_genes_groups(adata, 'condition', method=method)
sc.pl.rank_genes_groups(adata, groupby='condition', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_conditions')
###FIND UPREGULATED GENES IN EACH CLUSTER COMPARED TO ALL OTHER CLUSTERS - LOUVAIN:
sc.tl.rank_genes_groups(adata, 'louvain', method=method)
sc.pl.rank_genes_groups(adata, groupby='louvain', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters_louvain')
###FIND UPREGULATED GENES IN EACH CLUSTER SEPARATED SEPARATED BY TREATMENT CONDITION - LOUVAIN:
sc.tl.rank_genes_groups(adata, 'pairs', method=method)
sc.pl.rank_genes_groups(adata, groupby='pairs', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters_grouped_louvain')
###FIND UPREGULATED GENES IN EACH CLUSTER COMPARED TO ALL OTHER CLUSTERS - LOUVAIN:
sc.tl.rank_genes_groups(adata, 'leiden', method=method)
sc.pl.rank_genes_groups(adata, groupby='leiden', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters_leiden')
###FIND UPREGULATED GENES IN EACH CLUSTER SEPARATED SEPARATED BY TREATMENT CONDITION - LEIDEN:
sc.tl.rank_genes_groups(adata, 'pairs_leiden', method=method)
sc.pl.rank_genes_groups(adata, groupby='pairs_leiden', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters_grouped_leiden')


###CALCULATE CLUSTER STATISTICS USING HIGHLY VARIABLE GENES, MAKE TABLES OF MARKER GENES FOR EACH CLUSTER:
method = 'logreg' #t-test, wilcoxon, or logreg
cluster_method = 'leiden'
n_genes=4192 #set to adata.var.shape[0] to include all genes

if method=='logreg':
    sc.tl.rank_genes_groups(spadata, cluster_method, use_raw=False, pts=True, n_genes=n_genes, method=method, max_iter=1000)
else:
    sc.tl.rank_genes_groups(spadata, cluster_method, use_raw=False, pts=True, n_genes=n_genes, method=method)
sc.pl.rank_genes_groups(spadata, groupby=cluster_method, use_raw=False, n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_' +  cluster_method + '_clusters')
pd.DataFrame(spadata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(spadata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_' + cluster_method + '_cluster_table_' + str(n_genes) + 'genes.xlsx'), engine='openpyxl')
#make table with p-values included
result = spadata.uns['rank_genes_groups']
groups = result['names'].dtype.names
if method != 'logreg':
    pval_table = pd.DataFrame(
            {group + '_' + key[:2]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_filtered_corrected.xlsx'), engine='openpyxl')
elif method=='logreg':
        pval_table = pd.DataFrame(
            {group + '_' + key[:2]: result[key][group]
            for group in groups for key in ['names', 'scores']}).head(n_genes)
        pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_coefs_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_filtered_corrected.xlsx'), engine='openpyxl')

###CALCULATE CLUSTER STATISTICS USING ALL GENES, MAKE TABLES OF MARKER GENES FOR EACH CLUSTER:
method = 't-test' #t-test, wilcoxon, or logreg
cluster_method = 'leiden'
n_genes=2000 #set to adata.var.shape[0] to include all genes

if method=='logreg':
    sc.tl.rank_genes_groups(spadata, cluster_method, use_raw=True, pts=True, n_genes=n_genes, method=method, max_iter=1000)
else:
    sc.tl.rank_genes_groups(spadata, cluster_method, use_raw=True, pts=True, n_genes=n_genes, method=method)
sc.pl.rank_genes_groups(spadata, groupby=cluster_method, use_raw=True, n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_' +  cluster_method + '_clusters')
pd.DataFrame(spadata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(spadata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_' + cluster_method + '_cluster_table_' + str(n_genes) + 'genes_raw.xlsx'), engine='openpyxl')
#make table with p-values included
result = spadata.uns['rank_genes_groups']
groups = result['names'].dtype.names
if method != 'logreg':
    pval_table = pd.DataFrame(
            {group + '_' + key[:2]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_raw_corrected.xlsx'), engine='openpyxl')
elif method=='logreg':
        pval_table = pd.DataFrame(
            {group + '_' + key[:2]: result[key][group]
            for group in groups for key in ['names', 'scores']}).head(n_genes)
        pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_coefs_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_raw_corrected.xlsx'), engine='openpyxl')


###DIFFERENTIAL EXPRRESSION OF GENES WITHIN CLUSTERS:
pairs = list(zip(spadata.obs['condition'], spadata.obs['louvain'].astype('int')))
pairs_set = list(set(pairs))
s = sorted(pairs_set)
half = int((len(s)/2))
list1 = s[:half]
list2 = s[half:]
lz_louvain = list(zip(list1, list2))

pairs_leid = list(zip(spadata.obs['condition'], spadata.obs['leiden'].astype('int')))
pairs_leid_set = list(set(pairs_leid))
s_leid = sorted(pairs_leid_set)
half_leid = int((len(s_leid)/2)+1)
list1_leid = s_leid[:half_leid-1]
list2_leid = s_leid[half_leid-1:]
lz_leiden = list(zip(list1_leid, list2_leid))

#IMPORTANT: INSPECT LZ_LOUVAIN AND LZ_LEID TO MAKE SURE THEY ARE CORRECTLY ALIGNED. 
#EMPTY CLUSTERS IN ONE GROUP WILL CAUSE PROBLEMS.
lz_louvain

lz_leiden


###CALCULATE GENES UPREGULATED IN GROUP 2:
method = 't-test' #t-test, wilcoxon, or logreg
cluster_method = 'leiden'
n_genes = spadata.var.shape[0]

if cluster_method=='louvain':
    list2compare=lz_louvain
elif cluster_method=='leiden':
    list2compare=lz_leiden
else:
    raise ValueError("Invalid cluster method selected")

###CALCULATE GENES UPREGULATED IN GROUP 2 USING RAW DATA:
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(spadata, 'pairs_' + cluster_method, groups=[str(i[1])], reference=str(i[0]), use_raw=True, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(spadata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, use_raw=True, method=method)        
    result = spadata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_rawData_adj_new.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1 USING RAW DATA: 
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(spadata, 'pairs_' + cluster_method, groups=[str(i[0])], reference=str(i[1]), use_raw=True, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(spadata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, use_raw=True, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = spadata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_rawData_adj_new.xlsx'))

###CALCULATE GENES UPREGULATED IN GROUP 2 USING ONLY HIGHLY VARIABLE GENES:
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(spadata, 'pairs_' + cluster_method, groups=[str(i[1])], reference=str(i[0]), use_raw=False, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(spadata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, use_raw=False, method=method)        
    result = spadata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '_' + cluster_method + '_' + str(resolution) + 'resolution_' + str(n_genes) + 'genes_filteredHVGs_adj.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1 USING ONLY HIGHLY VARIABLE GENES: 
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(spadata, 'pairs_' + cluster_method, groups=[str(i[0])], reference=str(i[1]), use_raw=False, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(spadata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, use_raw=False, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = spadata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '_' + cluster_method + '_' + str(resolution) + 'resolution_' + str(n_genes) + 'genes_filteredHVGs_adj.xlsx'))

#COUNT DEGs:
file='/d1/studies/scanPy/spatial/March2021_FinalAnalysis/Spatial_Stressed_Cropped_DiffExp_UpregulatedStress_t-test_leiden_4192genes_rawData_adj.xlsx'
fcx.countDEGs(file, BaseDirectory, n_genes=4192, pcutoff=.05, plot=True, save=True, imageType='.svg')


###QUERY GENE ONTOLOGY DATABASE FOR ENRICHMENT OF DEGs - (OPTIONAL, DOES NOT AFFECT DOWNSTREAM ANALYSIS):
clus = adata.obs['leiden'].unique().tolist()
for cluster in clus:
    if method != 'logreg':
        sigGenes = cat.loc[cat[r'(0, '+cluster+r')'+'_p']<.05][r'(0, '+cluster+r')'+'_n'].tolist()
    if len(sigGenes)>0:
        enriched = sc.queries.enrich(sigGenes, org='mmusculus')
        if not os.path.exists(os.path.join(BaseDirectory, 'GOenrichment_'+str(g1n)+'/')):
            os.mkdir(os.path.join(BaseDirectory, 'GOenrichment_'+str(g1n)+'/'))
        enriched.to_excel(os.path.join(BaseDirectory, 'GOenrichment_'+str(g1n)+'/' + sampleName + '_GOenrichment_Cluster' + str(cluster) + '_Up' + str(g1n) + '.xlsx'), engine = 'openpyxl')
        categories = enriched['source'].unique().tolist()
        for categ in categories:
            df = enriched.loc[enriched['source']==categ]
            df.to_excel(os.path.join(BaseDirectory, 'GOenrichment_'+str(g1n)+'/' + sampleName + '_Enriched_' + categ + '_Cluster' + str(cluster) + '_Up' + str(g1n) + '.xlsx'), engine='openpyxl')
    else:
        print("No significant genes found in " + str(cluster))
        pass
    
#ENRICHED GROUP 2
clus = adata.obs['leiden'].unique().tolist()
for cluster in clus:
    if method != 'logreg':
        sigGenes = cat.loc[cat[r'(1, '+cluster+r')'+'_p']<.05][r'(1, '+cluster+r')'+'_n'].tolist()
    if len(sigGenes)>0:
        enriched = sc.queries.enrich(sigGenes, org='mmusculus')
        if not os.path.exists(os.path.join(BaseDirectory, 'GOenrichment_'+str(g2n)+'/')):
            os.mkdir(os.path.join(BaseDirectory, 'GOenrichment_'+str(g2n)+'/'))
        enriched.to_excel(os.path.join(BaseDirectory, 'GOenrichment_'+str(g2n)+'/' + sampleName + '_GOenrichment_Cluster' + str(cluster) + '_Up' + str(g2n) + '.xlsx'), engine = 'openpyxl')
        categories = enriched['source'].unique().tolist()
        for categ in categories:
            df = enriched.loc[enriched['source']==categ]
            df.to_excel(os.path.join(BaseDirectory, 'GOenrichment_'+str(g2n)+'/' + sampleName + '_Enriched_' + categ + '_Cluster' + str(cluster) + '_Up' + str(g2n) + '.xlsx'), engine='openpyxl')
    else:
        print("No significant genes found in " + str(cluster))
        pass


###CONGRATULATIONS, THE PRIMARY ANALYSIS IS DONE. THIS IS A GOOD PLACE TO SAVE YOUR RESULTS:
spadata.write(results_file)


###INTEGRATE NDB SINGLE-CELL DATA:
adata = sc.read('/d1/studies/scanPy/spatial/March2021_FinalAnalysis/NDB_scSeq/NDB_scSeq_throughUMAP.h5ad')


###THIS IS THE WAY:
shared_genes = adata.var_names.intersection(spadata.var_names)
adata_ing = adata[:, shared_genes].copy()
spadata_ing = spadata[:, shared_genes].copy()

adata_mapped = sc.tl.ingest(adata_ing, spadata_ing, obs='leiden', inplace=False)
adata_mapped.write('NDB_scSeq_MappedToSpadata.h5ad')

labeled_genes_var = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Pdgfra', 'Flt1']
labeled_genes = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Ntsr2', 'Pdgfra', 'Tmem119', 'Flt1', 'Pcbp3']
labeled_genes_var.insert(0, 'leiden')
labeled_genes.insert(0, 'leiden')
sc.pl.umap(adata_mapped, color=labeled_genes, color_map=color_map, vmax=vmax_raw, wspace=0.5,  save='_NDB_scSeq_MappedToSpadata_raw')
sc.pl.umap(adata_mapped, color=labeled_genes_var, color_map=color_map, vmin=vmin, vmax=vmax_filt, use_raw=False, wspace=0.5, save='_NDB_scSeq_MappedToSpadata_filtered')
sc.pl.umap(adata_mapped, color=['leiden'], use_raw=False, wspace=0.5, save='_NDB_scSeq_MappedToSpadata_leiden')

#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts = adata_mapped.obs['leiden'].value_counts().sort_index()
print(counts)

#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts_leid = spadata.obs['pairs_leiden'].value_counts().sort_index()
print(counts_leid)

###PLOT NUMBEROF CELLS IN EACH CLUSTER:
cdf = pd.DataFrame(counts)
cf = cdf.plot.bar(grid=False)
cf.set_ylabel('# Cells')
cf.set_xlabel('Leiden Cluster')
cf.set_title('Number of cells per cluster')
fig = cf.get_figure()
fig.savefig('figures/NDB_scSeq_adataMapped_CellsPercluster.' + sc.settings.file_format_figs)





###MERGE CLUSTERS 4 & 7:
spadata.obs.loc[spadata.obs['leiden']=='7']='4'
    
#############################################################################
###BELOW HERE ARE OPTIONAL ADDITIONAL ANALYSIS FUNCTIONS & SOME COOL PLOTTING.

###OPTIONALLY RECLUSTER A SINGLE CLUSTER INTO FURTHER SUBTYPES
cluster = 5
sc.tl.leiden(adata, restrict_to=('leiden', [str(cluster)]), resolution=resolution, n_iterations=iterations)
labeled_genes_var.insert(0, 'leiden_R')
labeled_genes.insert(0, 'leiden_R')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_recluster' + str(cluster))
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered_recluster' + str(cluster))


###EXTRACT MARKER GENES FROM CLUSTER TABLE GENERATED ABOVE, GENERATE NEW UMAP PLOT USING THE DATA-DRIVEN MARKERS AS LABELED GENES:

t = table.T
markers = t[0].tolist()
markers.insert(0,cluster_method)

###PLOT MARKERS AS ON UMAP:
#markers.insert(0,cluster_method)
vmin=None
vmax=5
color_map='plasma'
sc.pl.umap(spadata, color=markers, use_raw=False, color_map=color_map, vmin=vmin, vmax=vmax, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_' + cluster_method + '_filtered_MarkerGenes_' + str(method) + '_overlappingMarkers')
adata = sc.read('/d1/studies/scanPy/spatial/March2021_FinalAnalysis/NDB_scSeq_MappedToSpadata.h5ad')
sampleName2='NDB_scSeq_MappedToSpadata'
sc.pl.umap(adata, color=markers, use_raw=False, color_map=color_map, vmin=vmin, vmax=vmax, wspace=0.5, save='_' + str(sampleName2) + '_' + str(resolution) + 'resolution_clusters_labeled_' + cluster_method + '_filtered_MarkerGenes_' + str(method))



#markers = ['Slc17a7', 'Tac1', 'Zic1', 'Nptxr', 'Lhx8', 'Nrgn', 'Cldn11', 'Rxfp2', 'Hmgb2', 'Slc6a13']

markers = ['Vxn', 'Gpr88', 'Zic1', 'Rilpl1', 'Ache', 'Slc17a7', 'Mbp', 'Fam131c', 'Cdca7', 'Slc6a13']

cmap='Blues'
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', dendrogram=True, standard_scale='var', cmap=cmap, show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_PickedGenes7_withDend'+str(cmap))
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', dendrogram=False, use_raw=True, standard_scale='var', cmap=cmap, show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_PickedGenes7_raw_'+str(cmap))
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', vmin=0.3, vmax=1, use_raw=False, standard_scale='obs', cmap=cmap, show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_PickedGenes7_filtered'+str(cmap))

cmaps=['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r']
for cmap in cmaps:
    sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', cmap=cmap, show_gene_labels=True, save='_cmaps/' + str(sampleName) + '_' + method + '_clusters_leiden_PickedGenes6_'+str(cmap))

sc.pl.umap(spadata, color=markers, use_raw=False, wspace=0.5, color_map='plasma', save='_CustomPicked_MarkerGenes')
alpha=1
cmap='plasma'
vmax=None
sc.pl.spatial(stress1, library_id='stress1', img_key="hires", use_raw=False, vmax=vmax, alpha=alpha, color=markers, size=size, color_map=cmap, save='_spatialReconstruction_stress1_CustomPicked_MarkerGenes_alpha'+str(alpha)+'_filtered')


labeled_genes_var = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Pdgfra', 'Flt1']
labeled_genes = ['Drd1', 'Pdyn', 'Drd2', 'Penk', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Ntsr2', 'Pdgfra', 'Tmem119', 'Flt1', 'Pcbp3']
labeled_genes_var.insert(0, 'leiden')
labeled_genes.insert(0, 'leiden')

sc.pl.umap(adata, color=labeled_genes, use_raw=True, wspace=0.5, color_map='plasma', save='_adata_mappedToSpatial_LabeledGenes_raw')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, color_map='plasma', save='_adata_mappedToSpatial_LabeledGenes_raw')
markers = ['Neurod2', 'Gpr88', 'Zic1', 'Nptxr', 'Ache', 'Slc17a7', 'Mog', 'Coro6', 'Sox4', 'Slc6a13']
markers.insert(0,'leiden')
sc.pl.umap(adata, color=markers, use_raw=False, wspace=0.5, color_map='plasma', save='_adata_mappedToSpatial_LabeledGenes_Filtered')
sc.pl.umap(adata, color=['leiden'], use_raw=True, wspace=0.5, color_map='plasma', save='_adata_mappedToSpatial_Leiden')



###PLOT AS HEATMAP:
sc.pl.heatmap(adata, var_names=markers, groupby='louvain', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_louvain')
sc.pl.heatmap(adata, var_names=markers, groupby='pairs', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_pairs')
markers = t[0].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_PickedGenes')
markers = t[1].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden1')
markers = t[2].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden2')
markers = t[3].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden3')
markers = t[4].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden4')
markers = t[5].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden5')
markers = t[6].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden6')
markers = t[7].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden7')
markers = t[8].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden8')
markers = t[9].tolist()
sc.pl.heatmap(spadata, var_names=markers, groupby='leiden', standard_scale='var', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden9')

###HEIRARCHICAL CLUSTERING:
sc.tl.dendrogram(spadata, groupby='louvain')
sc.pl.dendrogram(spadata, groupby='louvain', save='_louvain')
sc.tl.dendrogram(spadata, groupby='pairs')
sc.pl.dendrogram(spadata, groupby='pairs', save='_pairs')
sc.tl.dendrogram(spadata, groupby='leiden')
sc.pl.dendrogram(spadata, groupby='leiden', save='_leiden')

markers_dend = ['Oprm1', 'Slc17a6', 'Xylt1', "Camk2a", 'Grik3', 'Pcbp3', 'Grik1', 'Rgs6', 'Fstl4', 'Luzp2', 'Hs6st3', 'Plp1', 'Pcbp3', 'Grm8', 'Inpp5d', 'Vcan', 'Arhgap6', 'Cpa6', 'Prex2', 'Flt1']


###PLOT CORRELATION MATRIX:
sc.pl.correlation_matrix(spadata, groupby='louvain', save='_louvain')
sc.pl.correlation_matrix(spadata, groupby='leiden', save='_leiden')
sc.pl.correlation_matrix(spadata, groupby='pairs', save='_pairs')
sc.pl.correlation_matrix(spadata, groupby='condition', save='_condition')

###CORRELATION OF EACH GROUP INDIVIDUALLY:
sc.pl.correlation_matrix(spadata_g1, show_correlation_numbers=True, groupby='leiden', save='_leiden_' + g1n)
sc.pl.correlation_matrix(spadata_g2, show_correlation_numbers=True, groupby='leiden', save='_leiden_' + g2n)
mtx_g1 = spadata_g1.uns['dendrogram_leiden']['correlation_matrix']
mtx_g2 = spadata_g2.uns['dendrogram_leiden']['correlation_matrix']
mtx_difference = mtx_g1 - mtx_g2
#PLOT DIFFERENCE BETWEEN TWO HEATMAPS:
order = spadata.uns['dendrogram_leiden']['categories_ordered']
fig, ax = plt.subplots(ncols=3, sharex=True, tight_layout=True, gridspec_kw=dict(width_ratios=[1,1,1]), figsize=(35,10))
fig.suptitle(g1n + ' minus ' + g2n + ' correlation matrix', fontsize=36)
sns.heatmap(mtx_g1, annot=True, square=True, cbar=False, vmin=-1, vmax=1, xticklabels=order, yticklabels=order, linewidths=1, linecolor='black', cmap='bwr', ax=ax[0])
sns.heatmap(mtx_g2, annot=True, square=True, cbar=False, vmin=-1, vmax=1, xticklabels=order, yticklabels=order, linewidths=1, linecolor='black', cmap='bwr', ax=ax[1])
sns.heatmap(mtx_difference, square=True, annot=True, vmin=-1, vmax=1, xticklabels=order, yticklabels=order, linewidths=1, linecolor='black', cmap='bwr', ax=ax[2])
ax[0].set_xticklabels(order, fontsize=24)
ax[0].set_yticklabels(order, fontsize=24)
ax[0].set_title(g1n, fontdict={'fontsize':32})
ax[1].set_xticklabels(order, fontsize=24)
ax[1].set_yticklabels(order, fontsize=24)
ax[1].set_title(g2n, fontdict={'fontsize':32})
ax[2].set_xticklabels(order, fontsize=24)
ax[2].set_yticklabels(order, fontsize=24)
ax[2].set_title('Difference', fontdict={'fontsize':32})
plt.show()
fig.savefig(os.path.join(BaseDirectory, 'figures/CorrelationMatrices_Combined.' + sc.settings.file_format_figs), bbox_inches='tight')

###PLOT EMBEDDING DENSITY:
sc.tl.embedding_density(spadata, basis='umap', groupby='condition', key_added='umap_density')
sc.pl.embedding_density(spadata, basis='umap', key='umap_density', title='UMAP_Density', save='all')
sc.pl.embedding_density(spadata, basis='umap', key='umap_density', title='UMAP_Density_' + g1n, group=[0], save='Group1')
sc.pl.embedding_density(spadata, basis='umap', key='umap_density', title='UMAP_Density_' + g2n, group=[1], save='Group2')


###PLOT PARTITION-BASED GRAPH ABSTRACTIONS
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, save='_paga_louvain')
sc.tl.paga(spadata, groups='leiden')
sc.pl.paga(spadata, save='_paga_leiden')
sc.pl.paga_compare(spadata, basis='umap', save='_paga_compare')

sc.tl.paga(adata, groups='pairs')
sc.pl.paga(adata, single_component=True, save='_paga_pairs_louvain')
sc.tl.paga(adata, groups='pairs_leiden')
sc.pl.paga(adata, single_component=True, save='_paga_pairs_leiden')

###PLOT DOT PLOTS:
adatac5 = adata[adata.obs['leiden']=='5']
mito=['mt-Nd1', 'mt-Nd2', 'mt-Nd3', 'mt-Nd4', 'mt-Nd4l', 'mt-Nd5', 'mt-Nd6', 'mt-Cytb', 'mt-Co1', 'mt-Co2', 'mt-Co3', 'mt-Atp6', 'mt-Atp8', 'Taco1', 'Atg4a']
ieg=['Fos', 'Arc', 'Npas4', 'Cux2', 'Egr1']
sc.pl.stacked_violin(adata, mito, groupby='condition', scale='width', standard_scale='var', use_raw=False, save='_' + str(sampleName) + '_MitoGenes_standard_filt')
sc.pl.dotplot(adatac5, mito, groupby='condition', standard_scale='obs', cmap='Blues', use_raw=False, save='_' + str(sampleName) + '_MitoGenes_filt_obs_c5')
sc.pl.dotplot(adatac5, ieg, groupby='condition', standard_scale='obs', cmap='Blues', use_raw=False, save='_' + str(sampleName) + '_IEGs_filt_obs_c5')

###PLOT STACKED VIOLIN:
sc.pl.stacked_violin(adata, var_names=labeled_genes_var, groupby='leiden', save='_leiden')

###WRITE RESULTS:
spadata.write(results_file)

###IF YOU LATER WANT TO REVERT TO A RESULTS FILE USE:
adata = sc.read(results_file)



