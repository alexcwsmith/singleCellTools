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
BaseDirectory = '/d1/studies/cellranger/ACWS_DP/scanpy_DiffExp_V2/'

###SET SCANPY SETTINGS:
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = os.path.join(BaseDirectory, 'scanpy_results.h5ad')  # the file that will store the analysis results
sc.settings.set_figure_params(dpi=80)


###LOAD DATA FROM MTX FILE:
adata = sc.read_10x_mtx(
    BaseDirectory,  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading
adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'


###ADD CONDITION IDs TO ADATA ANNOTATIONS:
#Import merged barcodes from extracted tsv file into pandas dataframe:
cells = pd.read_csv(os.path.join(BaseDirectory, 'barcodes.tsv', sep='\t', names=['barcode']))
#Preview the data that imported:
cells.head()
#Make a list of barcodes:
barcodes = cells.barcode.tolist()
#Initiate lists
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
#EDIT GROUP LISTS FOR YOUR SAMPLES. THIS WILL BE BASED ON THE ORDER OF THE SAMPLES IN THE CELLRANGER AGGR INDEX CSV:
g1 = ['1','2','3',]
g2 = ['4','5','6',]
group1 = anno.loc[anno['sample'].isin(g1)]
group1['condition']=1
group2 = anno.loc[anno['sample'].isin(g2)]
group2['condition']=2
#Put it all back together:
anno = pd.concat([group1, group2], axis=0)
#Add the group labels to the adata annotations:
adata.obs['condition'] = anno['condition'].values.astype('category')

###EXPLORE DATA, FILTER HIGHEST EXPRESSING GENES:
sc.pl.highest_expr_genes(adata, n_top=50, save='_highestExpressingGenes')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.calculate_qc_metrics(adata)
qc = sc.pp.calculate_qc_metrics(adata)

###CALCULATE % MITOCHONDRIAL GENES FOR EACH CELL, AND ADD TO ADATA.OBS:
mito_genes = adata.var_names.str.startswith('mt-') 
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


###PLOT % MITO & TOTAL COUNTS:
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save='_plot_percentMito.tif')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_mito_counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_genes_counts')

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save='_plot_percentMito.tif')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_mito_counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_genes_counts')


###SAVE THE CURRENT ADATA STATE, CAN REVERT TO THIS WHILE TUNING PARAMETERS:
adata_raw = adata

###FILTER - TUNE THESE PARAMETERS
max_genes = 4000
max_mito = 0.15
adata = adata[adata.obs['n_genes'] < max_genes, :]
adata = adata[adata.obs['percent_mito'] < max_mito, :]
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_filtered_mito_counts', title='Filtered < ' + str(max_genes) +  ' total counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_filtered_genes_counts_test', title='Filtered < ' + str(max_mito) + ' mito')

#NORMALIZE & LOG TRANSFORM
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)


###PRE-PROCESS DATA, SELECT HIGHLY VARIABLE GENES
min_mean = .0125
max_mean = 3
min_disp = 0.25
sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
"""This adds following dimensions to the data:
    'highly_variable', boolean vector (adata.var)
    'means', float vector (adata.var)
    'dispersions', float vector (adata.var)
    'dispersions_norm', float vector (adata.var)"""

sc.pl.highly_variable_genes(adata, save='_highlyVariableGenes')
adata = adata[:, adata.var['highly_variable']]
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

###WRITE EXCEL FILE OF HIGHLY VARIABLE GENES:
adata.var.to_excel(os.path.join(BaseDirectory, 'HighlyVariableGenes_minMean' + str(min_mean) + '_maxMean' + str(max_mean) + '_min_disp' + str(min_disp) + '.xlsx'))

#scale each gene to unit variance. Clip values exceeding standard deviation of 10:
sc.pp.scale(adata, max_value=10)
sc.pl.violin(adata, 'Oprm1', save='_Oprm1')

###CREATE LIST OF GENES TO LABEL ON PCA/UMAP PLOTS:
labeled_genes = ['Oprm1', 'Slc17a7', 'Slc17a6', 'Gad1', 'Slc1a2', 'Olig2', 'Tmem119', 'Dbi', 'Opcml', 'Fos']


###RUN PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color=labeled_genes, save='PCA_labeled')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)
sc.pl.pca_overview(adata, save='PCA_Overview')
sc.pl.pca_loadings(adata, components=list(np.arange(30)), save='_allPCs')

#COMPUTING NEIGHBORHOOD GRAPH:
#Uses PCA representation of data matrix
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=labeled_genes, save='_umap_labeled_raw') #'Slc17a7', 'Slc17a6', 'Cxcr6', 'Ache', 'Olig1'

###THIS REQUIRES GENES IN LABELED_GENES TO BE IN HIGHLY VARIABLE LIST. NOT SUPER USEFUL:
sc.pl.umap(adata, color=labeled_genes, use_raw=False, save='_umap_filtered')

###CLUSTER DATA:
clu = sc.tl.louvain(adata)
labeled_genes.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_clusters_labeled')

###SEPARATE CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
adata.obs['pairs'] = pairs
adata.obs['louvain'] = adata.obs['louvain'].values.remove_unused_categories()


###RUN STATS, NOTE 't-test' can be changed to 'wilcoxon':
###COMPARE EXPRESSION BY CONDITION (EQUIVALENT TO BULK-SEQ):
sc.settings.verbosity = 2
sc.tl.rank_genes_groups(adata, 'condition', method='t-test')
sc.pl.rank_genes_groups(adata, groupby='condition', n_genes=25, sharey=False, save='_t-test_conditions')
###FIND UPREGULATED GENES IN EACH CLUSTER COMPARED TO ALL OTHER CLUSTERS:
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, groupby='louvain', n_genes=25, sharey=False, save='_t-test_clusters')
###FIND UPREGULATED GENES IN EACH CLUSTER SEPARATED BY TREATMENT CONDITION:
sc.tl.rank_genes_groups(adata, 'pairs', method='t-test')
sc.pl.rank_genes_groups(adata, groupby='pairs', n_genes=25, sharey=False, save='_t-test_clusters_grouped')

###WRITE RESULTS:
adata.write(results_file)

###IF YOU LATER WANT TO REVERT TO A RESULTS FILE USE:
adata.read(results_file)

###MAKE TABLES OF GENES IN EACH CLUSTER
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, groupby='louvain', n_genes=25, sharey=False, save='_t-test_clusters')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, '_t-test_condition_table.xlsx'), engine='openpyxl')
#make table with p-values included
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
pval_table = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(500) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
pval_table.to_excel(os.path.join(BaseDirectory, 't-test_pval_table_500genes_conditions.xlsx'), engine='openpyxl')


###DIFFERENTIAL EXPRRESSION OF GENES WITHIN CLUSTERS:
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
pairs_set = list(set(pairs))
s = sorted(pairs_set)
half = int((len(s)/2)+1)
list1 = s[:half]
list2 = s[half-1:]
lz = list(zip(list1, list2))

###CALCULATE GENES UPREGULATED IN GROUP 1:
cat = pd.DataFrame()
for i in lz:
    sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=500, method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel('DiffExp_UpregulatedGroup2.xlsx')
###CALCULATE GENES UPREGULATED IN GROUP 2: 
cat = pd.DataFrame()
for i in lz:
    sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=500, method='wilcoxon')
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(500)
 #   pval_table.columns=[(str(i[0][0]) + '_gene'), (str(i[0][0]) + '_pvalue')]
 #   pval_table.to_excel(os.path.join(BaseDirectory, 't-test_pval_table_500genes_conditions' + str(i) + '.xlsx'), engine='openpyxl')

    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel('DiffExp_UpregulatedGroup1.xlsx')








