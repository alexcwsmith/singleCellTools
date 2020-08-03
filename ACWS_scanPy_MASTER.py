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
from skimage import io
import os
import matplotlib
import openpyxl

###SET DIRECTORY TO READ/WRITE DATA. SET THE SPYDER WORKING DIRECTORY TO THE SAME PATH (TOP RIGHT OF SPYDER).
#THIS SHOULD BE THE DIRECTORY CONTAINING THE .MTX DATA FILE AND .TSV BARCODES & FEATURE FILES:
BaseDirectory = '/d1/studies/cellranger/ACWS_DP/Batch2_Horde_DP_Sal_Merged/'
sampleName = 'Batch2_Horde_DP_OC_Sal' #This is used for name result output files

###SET SCANPY SETTINGS:
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = os.path.join(BaseDirectory, sampleName + '_scanpy_results.h5ad')  # the file that will store the analysis results
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='png')


###LOAD DATA FROM MTX FILE:
adata = sc.read_10x_mtx(
    BaseDirectory,  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading
adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'


###ADD CONDITION IDs TO ADATA ANNOTATIONS:
#Import merged barcodes from extracted tsv file into pandas dataframe:
cells = pd.read_csv(os.path.join(BaseDirectory, 'barcodes.tsv'), sep='\t', names=['barcode'])
#Preview the data that imported:
cells.head()
#Make a list of barcodes:
barcodes = cells.barcode.tolist()
#Initialize result lists
cellsList = []
groupList=[]
#Split the barcodes from groups and append to their own lists:
for cell in barcodes:
    code, group = cell.split(sep='-')
    cellsList.append(code)
    groupList.append(group)
#Create a new pandas dataframe with the split data:
anno = pd.DataFrame(data=(cellsList, groupList)).T
anno.columns=['barcode', 'sample']
#Split the data into groups & annotate - 
#EDIT GROUP LISTS FOR YOUR SAMPLES. THIS SHOULD BE BASED ON THE ORDER OF THE SAMPLES IN THE CELLRANGER AGGR INDEX CSV, BUT YOU NEED TO DOUBLE CHECK INDIVIDUAL BARCODE TSVs:
g1 = ['1','2', '4']
g2 = ['3','5', '6', '7', '8']
group1 = anno.loc[anno['sample'].isin(g1)]
group1['condition']=1
group1['condition'] = group1['condition'].astype('category')
group2 = anno.loc[anno['sample'].isin(g2)]
group2['condition']=2
group2['condition'] = group2['condition'].astype('category')
#Put it all back together:
anno = pd.concat([group1, group2], axis=0)
#Add the group labels to the adata annotations:
adata.obs['condition'] = anno['condition'].values

#Name the groups, for naming of result output files only:
g1n = 'Control'
g2n = 'Treated'

groupNames = {"1" : g1n, "2" : g2n}

###IF SAMPLES WERE RUN IN MULTIPLE BATCHES (i.e. TAKEN TO THE CORE AT SEPARATE TIMES) ADD BATCH INFO TO ADATA:
###IF NO BATCH CORRECTION, SKIP TO LINE 95
b1 = ['1','2','4']
b2 = ['3','5','6']
b3 = ['7','8']
batch1 = anno.loc[anno['sample'].isin(b1)]
batch1['batch']=1
batch2 = anno.loc[anno['sample'].isin(b2)]
batch2['batch']=2
batch3 = anno.loc[anno['sample'].isin(b3)]
batch3['batch']=3
batches = pd.concat([batch1, batch2, batch3])

adata.obs['batch']=batches['batch'].values
adata.obs['batch']=adata.obs['batch'].astype('category')


###EXPLORE DATA, FILTER HIGHEST EXPRESSING GENES:
sc.pl.highest_expr_genes(adata, n_top=50, save='_' + str(sampleName) + '_highestExpressingGenes')
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata)

###CALCULATE % MITOCHONDRIAL GENES FOR EACH CELL, AND ADD TO ADATA.OBS:
mito_genes = adata.var_names.str.startswith('mt-') 
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


###PLOT % MITO & TOTAL COUNTS:
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save='_' + str(sampleName) + '_plot_percentMito')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_' + str(sampleName) + '_mito_counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_' + str(sampleName) + '_genes_counts')


###SAVE THE CURRENT ADATA STATE, CAN REVERT TO THIS WHILE TUNING PARAMETERS:
adata.raw = adata

###FILTER - TUNE THESE PARAMETERS
max_genes = 3500 #Look at the genes_counts output figure from previous step to decide where you want your cutoff.
max_mito = 0.05 #Look at percentMito output figure from previous step to decide on cutoff.
adata = adata[adata.obs['n_genes'] < max_genes, :]
adata = adata[adata.obs['percent_mito'] < max_mito, :]
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_' + str(sampleName) + '_filtered_mito_counts', title='Filtered < ' + str(max_genes) +  ' total counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_' + str(sampleName) + '_filtered_genes_counts_test', title='Filtered < ' + str(max_mito) + ' mito')

###NORMALIZE & LOG TRANSFORM
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)


###PRE-PROCESS DATA, SELECT HIGHLY VARIABLE GENES. 
###YOU WILL LIKELY WANT TO PLAY WITH THESE AND SEE HOW IT AFFECTS RESULTS.
min_mean = .005
max_mean = 3
min_disp = 0.15
sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
"""This adds following dimensions to the data:
    'highly_variable', boolean vector (adata.var)
    'means', float vector (adata.var)
    'dispersions', float vector (adata.var)
    'dispersions_norm', float vector (adata.var)"""
    
###PLOT DISPERSIONS AND MEANS TO HELP REFINE:
sc.pl.scatter(adata, x='n_cells', y='dispersions', save='_' + str(sampleName) + '_dispersions', title='Dispersions')
sc.pl.scatter(adata, x='n_cells', y='dispersions_norm', save='_' + str(sampleName) + '_dispersions_norm', title='Dispersions Normalized')
sc.pl.scatter(adata, x='n_cells', y='means', save='_' + str(sampleName) + '_means', title='Means')

###WRITE TABLE OF RESULTS BEFORE FILTERING:
adata.var.to_excel(os.path.join(BaseDirectory, 'Adata_var_raw_preFiltering.xlsx'))

###ACTUALLY DO THE FILTERING
sc.pl.highly_variable_genes(adata, save='_' + str(sampleName) + '_highlyVariableGenes')
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

#Batch effect correction:
sc.pp.combat(adata, key='batch')

###WRITE EXCEL FILE OF HIGHLY VARIABLE GENES:
adata.var.to_excel(os.path.join(BaseDirectory, str(sampleName) + 'HighlyVariableGenes_minMean' + str(min_mean) + '_maxMean' + str(max_mean) + '_min_disp' + str(min_disp) + '_batchCorrected.xlsx'))

###SCALE EACH GENE TO UNIT OF VARIANCE, CLIP VALUES EXCEEDING MAX VARIANCE:
sc.pp.scale(adata, max_value=10)
###EXAMPLE HOW TO PLOT EXPRESSION DISTRIBUTION INDIVIDUAL GENES:
sc.pl.violin(adata, 'Oprm1', save='_' + str(sampleName) + '_Oprm1')
   #You can also pass a list of genes here by just converting second argument to a list:
sc.pl.violin(adata, ['Oprm1', 'Penk'], save='_' + str(sampleName) + '_Oprm1_Penk')

###CREATE LIST OF GENES TO LABEL ON PCA/UMAP PLOTS:
labeled_genes = ['Oprm1', 'Slc17a7', 'Slc17a6', 'Slc1a2', 'Aldh1l1', 'Gad1', 'Gad2', 'Tmem119', 'Ache',  'Olig2', 'Dbi', 'Foxj1', 'Cx3cr1']

#RUN PCA
n_comps = 25
sc.tl.pca(adata, n_comps=n_comps, svd_solver='arpack')
sc.pl.pca(adata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_comps) + 'comps_PCA_labeled')
sc.pl.pca_variance_ratio(adata, log=True, save='_' + str(n_comps) + '_ncomps_PCA_VarianceRatio')
###LOOK AT VARIANCE_RATIO OUTPUT FILE BEFORE PROCEEDING TO NEXT STEP
n_pcs = 10
sc.pl.pca_overview(adata, save='_' + str(sampleName) + '_' + str(n_pcs) + 'PCs_PCA_Overview')
sc.pl.pca_loadings(adata, components=list(np.arange(1, n_pcs+1)), save='_' + str(sampleName) + '_' + str(n_pcs) + '_PCs')

#COMPUTING NEIGHBORHOOD GRAPH:
#Uses PCA representation of data matrix
###LOOK AT SAVED FIGURE W/ SUFFIX _PCA_VarianceRatio AND CHOOSE NUMBER OF PCs BEFORE APEX (HERE ~20)
n_neighbors = 50
n_pcs = 15
min_dist = .05
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(adata, min_dist=min_dist)
sc.pl.umap(adata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_neighbors) + 'neighbors_' + str(n_pcs) + 'PCs_umap_' + str(min_dist) + 'minDist_' + 'labeled_raw')
###THIS REQUIRES GENES IN LABELED_GENES TO BE IN HIGHLY VARIABLE LIST:
labeled_genes_var = ['Oprm1', 'Penk', 'Slc17a6', 'Rorb', 'Slc1a2', 'Aldh1l1', 'Gfap', 'Gad1', 'Gad2', 'Olig2', 'Dbi', 'Cx3cr1']
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, save='_' + str(sampleName) + '_nNeighbors' + str(n_neighbors) + '_nPCs' + str(n_pcs) + '_umap_filtered')

###CLUSTER DATA:
sc.tl.louvain(adata)
labeled_genes.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + '_clusters_labeled')
labeled_genes_var.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_clusters_labeled_filtered')



###SEPARATE CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
adata.obs['pairs'] = pairs
adata.obs['louvain'] = adata.obs['louvain'].values.remove_unused_categories()

adata.obs.condition=adata.obs.condition.astype('category')
###RUN STATS, NOTE 't-test' can be changed to 'wilcoxon':
###COMPARE EXPRESSION BY CONDITION (EQUIVALENT TO BULK-SEQ):
sc.settings.verbosity = 2
method = 't-test' #t-test, wilcoxon, or logreg

sc.tl.rank_genes_groups(adata, 'condition', method=method)
sc.pl.rank_genes_groups(adata, groupby='condition', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_conditions')
###FIND UPREGULATED GENES IN EACH CLUSTER COMPARED TO ALL OTHER CLUSTERS:
sc.tl.rank_genes_groups(adata, 'louvain', method=method)
sc.pl.rank_genes_groups(adata, groupby='louvain', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters')
###FIND UPREGULATED GENES IN EACH CLUSTER SEPARATED SEPARATED BY TREATMENT CONDITION:
sc.tl.rank_genes_groups(adata, 'pairs', method=method)
sc.pl.rank_genes_groups(adata, groupby='pairs', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters_grouped')

###WRITE RESULTS:
adata.write(results_file)

###IF YOU LATER WANT TO REVERT TO A RESULTS FILE USE:
adata.read(results_file)

###MAKE TABLES OF GENES IN EACH CLUSTER
method = 'wilcoxon' #t-test, wilcoxon, or logreg

sc.tl.rank_genes_groups(adata, 'louvain', n_genes=500, method=method)
sc.pl.rank_genes_groups(adata, groupby='louvain', n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_clusters')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_cluster_table.xlsx'), engine='openpyxl')
#make table with p-values included
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pval_table = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_pval_table_500genes_clusters.xlsx'), engine='openpyxl')


###DIFFERENTIAL EXPRRESSION OF GENES WITHIN CLUSTERS:
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
pairs_set = list(set(pairs))
s = sorted(pairs_set)
half = int((len(s)/2)+1)
list1 = s[:half]
list2 = s[half-1:]
lz = list(zip(list1, list2))

###CALCULATE GENES UPREGULATED IN GROUP 2:
method = 'wilcoxon' #t-test, wilcoxon, or logreg

cat = pd.DataFrame()
for i in lz:
    sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=500, method=method)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1: 
cat = pd.DataFrame()
for i in lz:
    sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=500, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '.xlsx'))


###PLOT AS HEATMAP:
t = table.T
markers = t[0].tolist()
sc.pl.heatmap(adata, var_names=markers, groupby='louvain', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters')
sc.pl.heatmap(adata, var_names=markers, groupby='pairs', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_pairs')

###HEIRARCHICAL CLUSTERING:
sc.tl.dendrogram(adata, groupby='louvain')
sc.pl.dendrogram(adata, groupby='louvain', orientation='left', save='_louvain')
sc.tl.dendrogram(adata, groupby='pairs')
sc.set_figure_params(fontsize=8, dpi=80, dpi_save=300, format='png')
sc.pl.dendrogram(adata, groupby='pairs', orientation='left', save='_pairs')

###PLOT CORRELATION MATRIX:
sc.pl.correlation_matrix(adata, groupby='louvain', save='_louvain')
sc.pl.correlation_matrix(adata, groupby='pairs', dendrogram=True, save='_pairs')
sc.pl.correlation_matrix(adata, groupby='condition', save='_condition')





