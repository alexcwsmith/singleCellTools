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
BaseDirectory = '/d1/studies/cellranger/ACWS_DP/scanpy_DiffExp_V8/'
sampleName = 'DP_OCvsSalineV8' #This is used for name result output files
batches = False #Set to True if you need to do batch correction (i.e. if samples were taken to core and sequenced at different times)

###SET SCANPY SETTINGS:
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
results_file = os.path.join(BaseDirectory, sampleName + 'scanpy_results.h5ad')  # the file that will store the analysis results
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=600, format='svg')


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
#EDIT GROUP LISTS FOR YOUR SAMPLES. THIS SHOULD BE BASED ON THE ORDER OF THE SAMPLES IN THE CELLRANGER AGGR INDEX CSV, BUT YOU NEED TO DOUBLE CHECK INDIVIDUAL BARCODE TSVs:
g1 = ['1','2','3']
g2 = ['4','5','6']
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
g1n = 'Saline'
g2n = 'OC'

groupNames = {"1" : g1n, "2" : g2n}

###IF SAMPLES WERE RUN IN MULTIPLE BATCHES (i.e. TAKEN TO THE CORE AT SEPARATE TIMES) ADD BATCH INFO TO ADATA:
if batches:
    b1 = ['1','3']
    b2 = ['2','4','5','6']
    batch1 = anno.loc[anno['sample'].isin(b1)]
    batch1['batch']=1
    batch2 = anno.loc[anno['sample'].isin(b2)]
    batch2['batch']=2
    batches = pd.concat([batch1, batch2])
    
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

###CALCULATE MEAN & SD OF PERCENT_MITO & N_COUNTS
mito_mean = np.mean(adata.obs.percent_mito)
mito_std = np.std(adata.obs.percent_mito)
mito_suggestedMax = mito_mean + mito_std
print("Suggested max Mito% = " + str(mito_suggestedMax))

n_genes_mean = np.mean(adata.obs.n_genes)
n_genes_std = np.std(adata.obs.n_genes)
n_genes_suggestedMax = n_genes_mean + n_genes_std
print("Suggested max genes = " + str(n_genes_suggestedMax))

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
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_' + str(sampleName) + '_filtered_genes_counts', title='Filtered < ' + str(max_mito) + ' mito')

###NORMALIZE & LOG TRANSFORM
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)


###PRE-PROCESS DATA, SELECT HIGHLY VARIABLE GENES. 
###YOU WILL LIKELY WANT TO PLAY WITH THESE AND SEE HOW IT AFFECTS RESULTS.
sc.pp.highly_variable_genes(adata)

###THE BELOW SUGGESTIONS ARE IN TESTING PHASE. THEY WORKED WELL ON MY DATA, WON'T NECESSARILY FOR YOURS:
genes_min_percentile = 50
genes_min_mean = np.percentile(adata.var.means, genes_min_percentile)
print("Suggested min_mean = " + str(genes_min_mean))

disp_mean = np.mean(adata.var.dispersions)
disp_std = np.std(adata.var.dispersions)
suggestedMinDisp = disp_mean - disp_std
print("Suggested min_disp = " + str(suggestedMinDisp))

min_mean = .01
max_mean = 3.5
min_disp = 0.2
if batches:
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp, batch_key='batch')
elif not batches:
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
if batches:
    sc.pp.combat(adata, key='batch')

###WRITE EXCEL FILE OF HIGHLY VARIABLE GENES:
adata.var.to_excel(os.path.join(BaseDirectory, str(sampleName) + 'HighlyVariableGenes_minMean' + str(min_mean) + '_maxMean' + str(max_mean) + '_min_disp' + str(min_disp) + '.xlsx'))

###SCALE EACH GENE TO UNIT OF VARIANCE, CLIP VALUES EXCEEDING MAX VARIANCE:
sc.pp.scale(adata, max_value=10)
###EXAMPLE HOW TO PLOT EXPRESSION DISTRIBUTION INDIVIDUAL GENES:
sc.pl.violin(adata, 'Oprm1', save='_' + str(sampleName) + '_Oprm1')
#You can also pass a list of genes here by just converting second argument to a list:
sc.pl.violin(adata, ['Oprm1', 'Penk'], save='_' + str(sampleName) + '_Oprm1_Penk')

###CREATE LIST OF GENES TO LABEL ON PCA/UMAP PLOTS:
labeled_genes_var = ['Oprm1', 'Slc17a6', 'Grik3', 'Gad1', 'Gad2', 'Slc4a4', 'Cpa6', 'Ntsr2', 'Pdgfra', 'Luzp2', 'C1qc', 'Flt1', 'Pcbp3', 'Dbi']
labeled_genes = ['Oprm1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Ntsr2', 'Pdgfra', 'Gpr17', 'Tmem119', 'C1qc', 'Cldn5', 'Flt1', 'Dbi']


#RUN PCA
n_comps = 30
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
n_neighbors = 25
n_pcs = 10
min_dist = .1
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(adata, min_dist=min_dist)
sc.pl.umap(adata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_neighbors) + 'neighbors_' + str(n_pcs) + 'PCs_umap_' + str(min_dist) + 'minDist_' + '_labeled_raw')
###THIS REQUIRES GENES IN LABELED_GENES TO BE IN HIGHLY VARIABLE LIST:
labeled_genes_var = ['Oprm1', 'Slc17a6', 'Grik3', 'Gad1', 'Gad2', 'Slc4a4', 'Cpa6', 'Ntsr2', 'Pdgfra', 'Luzp2', 'C1qc', 'Flt1', 'Pcbp3', 'Dbi']
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, save='_' + str(sampleName) + '_nNeighbors' + str(n_neighbors) + '_nPCs' + str(n_pcs) + '_umap_filtered')

###CLUSTER DATA WITH LOUVAIN ALGORITHM:
resolution = 0.9
sc.tl.louvain(adata, resolution=resolution)
labeled_genes.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution' + '_clusters_labeled_louvain')
labeled_genes_var.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution' + '_clusters_labeled_louvain_filtered')

###CLUSTER WITH LEIDEN ALGORITHM:
resolution = 0.9
sc.tl.leiden(adata, resolution=resolution)
labeled_genes.insert(0,'leiden')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution' + '_clusters_labeled_leiden')
labeled_genes_var.insert(0,'leiden')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution' + '_clusters_labeled_leiden_filtered')
sc.pl.umap(adata, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution' + '_clusters_labeled_leiden_only')


###SEPARATE LOUVAIN CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
adata.obs['pairs'] = pairs
adata.obs['louvain'] = adata.obs['louvain'].values.remove_unused_categories()

###SEPARATE LEIDEN CLUSTERS BY CONDITION & APPEND TO ADATA.OBS
pairs_leid = list(zip(adata.obs['condition'], adata.obs['leiden'].astype('int')))
adata.obs['pairs_leiden'] = pairs_leid
adata.obs['leiden'] = adata.obs['leiden'].values.remove_unused_categories()


#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts = adata.obs['pairs'].value_counts().sort_index()
print(counts)

#COUNT NUMBER OF CELLS IN EACH CLUSTER:
counts_leid = adata.obs['pairs_leiden'].value_counts().sort_index()
print(counts_leid)


###RUN STATS, NOTE 't-test' can be changed to 'wilcoxon':
###COMPARE EXPRESSION BY CONDITION (EQUIVALENT TO BULK-SEQ):
adata.obs.condition=adata.obs.condition.astype('category')
sc.settings.verbosity = 2
method = 't-test' #t-test, wilcoxon, or logreg

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




###MAKE TABLES OF GENES IN EACH CLUSTER
method = 't-test' #t-test, wilcoxon, or logreg
cluster_method = 'leiden'
n_genes=3696


sc.tl.rank_genes_groups(adata, cluster_method, n_genes=n_genes, method=method)
sc.pl.rank_genes_groups(adata, groupby=cluster_method, n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_' +  cluster_method + '_clusters')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_' + cluster_method + '_cluster_table_' + str(n_genes) + 'genes.xlsx'), engine='openpyxl')
#make table with p-values included
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(3696) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
pval_table = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(3696) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes.xlsx'), engine='openpyxl')


###DIFFERENTIAL EXPRRESSION OF GENES WITHIN CLUSTERS:
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
pairs_set = list(set(pairs))
s = sorted(pairs_set)
half = int((len(s)/2))
list1 = s[:half]
list2 = s[half:]
lz = list(zip(list1, list2))

pairs_leid = list(zip(adata.obs['condition'], adata.obs['leiden'].astype('int')))
pairs_leid_set = list(set(pairs_leid))
s_leid = sorted(pairs_leid_set)
half_leid = int((len(s_leid)/2)+1)
list1_leid = s_leid[:half_leid-1]
list2_leid = s_leid[half_leid-1:]
#list2_leid.insert(8,(2,8))
#list2.insert(6,(2,6))
#list2.insert(7,(2,7))
lz_leid = list(zip(list1_leid, list2_leid))


###CALCULATE GENES UPREGULATED IN GROUP 2:
method = 't-test' #t-test, wilcoxon, or logreg
cluster_method = 'louvain'
n_genes = 3696

cat = pd.DataFrame()
for i in lz_leid:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, method=method)        
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1: 
cat = pd.DataFrame()
for i in lz_leid:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes.xlsx'))


###PLOT AS HEATMAP:
t = table.T
markers = t[0].tolist()
sc.pl.heatmap(adata, var_names=markers, groupby='louvain', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_louvain')
sc.pl.heatmap(adata, var_names=markers, groupby='pairs', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_pairs')
sc.pl.heatmap(adata, var_names=markers, groupby='leiden', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden')

markers = t[18].tolist()
sc.pl.heatmap(adata, var_names=labeled_genes, groupby='leiden', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_labeledGenes')

labeled_genes = ['Oprm1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Ntsr2', 'Pdgfra', 'Gpr17', 'Tmem119', 'C1qc', 'Cldn5', 'Flt1', 'Dbi']


markers = ['Oprm1', 'Cmss1', 'Xylt1', "Camk2a", 'Grik3', 'Pcbp3', 'Grik1', 'Cux2', 'Rgs20', 'Luzp2', 'Hs6st3', 'Plp1', 'Pcbp3', 'Grm8', 'Inpp5d', 'Vcan', 'Arhgap6', 'Cpa6', 'Prex2', 'Flt1']

###HEIRARCHICAL CLUSTERING:
sc.tl.dendrogram(adata, groupby='louvain')
sc.pl.dendrogram(adata, groupby='louvain', save='_louvain')
sc.tl.dendrogram(adata, groupby='pairs')
sc.pl.dendrogram(adata, groupby='pairs', save='_pairs')
sc.tl.dendrogram(adata, groupby='leiden', var_names=markers_dend)
sc.pl.dendrogram(adata, groupby='leiden', save='_leiden')

markers_dend = ['Oprm1', 'Slc17a6', 'Xylt1', "Camk2a", 'Grik3', 'Pcbp3', 'Grik1', 'Rgs6', 'Fstl4', 'Luzp2', 'Hs6st3', 'Plp1', 'Pcbp3', 'Grm8', 'Inpp5d', 'Vcan', 'Arhgap6', 'Cpa6', 'Prex2', 'Flt1']



###PLOT CORRELATION MATRIX:
sc.pl.correlation_matrix(adata, groupby='louvain', save='_louvain')
sc.pl.correlation_matrix(adata, groupby='leiden', save='_leiden')
sc.pl.correlation_matrix(adata, groupby='pairs', save='_pairs')
sc.pl.correlation_matrix(adata, groupby='condition', save='_condition')

sc.tl.embedding_density(adata, basis='umap', groupby='condition', key_added='umap_density')
sc.pl.embedding_density(adata, basis='umap', key='umap_density', title='UMAP_Density', save='_all')
sc.pl.embedding_density(adata, basis='umap', key='umap_density', group=2, save='_group2')

sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, save='_paga')

sc.tl.paga(adata, groups='pairs_leiden')
sc.pl.paga(adata, single_component=True, save='_paga_pairs_leiden')

sc.pl.stacked_violin(adata, var_names=labeled_genes_var, groupby='leiden', save='_leiden')

###WRITE RESULTS:
adata.write(results_file)

###IF YOU LATER WANT TO REVERT TO A RESULTS FILE USE:
adata = sc.read(results_file)



