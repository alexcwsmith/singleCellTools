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
%logstart -o scanPy_log.txt

###SET DIRECTORY TO READ/WRITE DATA. SET THE SPYDER WORKING DIRECTORY TO THE SAME PATH (TOP RIGHT OF SPYDER).
#THIS SHOULD BE THE DIRECTORY CONTAINING THE .MTX DATA FILE AND .TSV BARCODES & FEATURE FILES:
BaseDirectory = '/d1/studies/cellranger/ACWS_DP/scanPy_V11/'
sampleName = 'DP_OCvsSalineV11' #This is used for name result output files
batches = False #Set to True if you need to do batch correction (i.e. if samples were taken to core and sequenced at different times)

###SET SCANPY SETTINGS:
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.n_jobs=8 #use parallel processing when possible
sc.logging.print_header()
results_file = os.path.join(BaseDirectory, sampleName + 'scanpy_results.h5ad')  # the file that will store the analysis results
results_file_partial = os.path.join(BaseDirectory, sampleName + 'scanpy_adata_prefiltering.h5ad')
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
    batches_combined = pd.concat([batch1, batch2])
    
    adata.obs['batch']=batches_combined['batch'].values
    adata.obs['batch']=adata.obs['batch'].astype('category')


###EXPLORE DATA, PLOT HIGHEST EXPRESSING GENES, FILTER LOW EXPRESSION GENES & CELLS:
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

n_genes_median = np.median(adata.obs.n_genes)
print("Suggested max genes based on median = " + str(n_genes_median + n_genes_std))


###PLOT % MITO & TOTAL COUNTS:
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True, save='_' + str(sampleName) + '_plot_percentMito')
sc.pl.scatter(adata, x='n_counts', y='percent_mito', save='_' + str(sampleName) + '_mito_counts')
sc.pl.scatter(adata, x='n_counts', y='n_genes', save='_' + str(sampleName) + '_genes_counts')


###SAVE THE CURRENT ADATA STATE, CAN REVERT TO THIS WHILE TUNING PARAMETERS:
adata.write(results_file_partial)
###IF YOU LATER WANT TO REVERT TO THE PRE-FILTERED RESULTS FILE USE:
adata = sc.read(results_file_partial)


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

###STORE THE RAW ADATA STATE AS A VARIABLE TO BE ABLE TO LABEL GENES THAT WERE FILTERED OUT ETC:
adata.raw = adata

###PRE-PROCESS DATA, SELECT HIGHLY VARIABLE GENES. 
###YOU WILL LIKELY WANT TO PLAY WITH THESE AND SEE HOW IT AFFECTS RESULTS.
sc.pp.highly_variable_genes(adata)

###THE BELOW SUGGESTIONS ARE IN TESTING PHASE. THEY WORKED WELL ON MY DATA, WON'T NECESSARILY FOR YOURS:
genes_min_percentile = 50
genes_min_mean = np.percentile(adata.var.means, genes_min_percentile)
print("Suggested min_mean = " + str(genes_min_mean))


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
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs=8)

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
ieg=['Fos', 'Arc', 'Npas4', 'Cux2', 'Egr1', 'Oprm1', 'Slc17a6']
sc.pl.stacked_violin(adata, ieg, groupby='condition', num_categories=2, standard_scale='obs', save='_' + str(sampleName) + 'IEGs')


###CREATE LIST OF GENES TO LABEL ON PCA/UMAP PLOTS:
labeled_genes_var = ['Oprm1', 'Slc17a6', 'Grik3', 'Gad1', 'Gad2', 'Slc4a4', 'Cpa6', 'Ntsr2', 'Pdgfra', 'Luzp2', 'C1qc', 'Flt1', 'Pcbp3', 'Dbi']
labeled_genes = ['Oprm1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Slc4a4', 'Ntsr2', 'Pdgfra', 'Gpr17', 'Tmem119', 'C1qc', 'Cldn5', 'Flt1', 'Dbi']


#RUN PCA
n_comps = 30
sc.tl.pca(adata, n_comps=n_comps, svd_solver='arpack')
sc.pl.pca(adata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_comps) + 'comps_PCA_labeled')
sc.pl.pca_variance_ratio(adata, log=True, save='_' + str(n_comps) + '_ncomps_PCA_VarianceRatio')
###LOOK AT VARIANCE_RATIO OUTPUT FILE BEFORE PROCEEDING TO NEXT STEP
n_pcs = 15
sc.pl.pca_overview(adata, save='_' + str(sampleName) + '_' + str(n_pcs) + 'PCs_PCA_Overview')
sc.pl.pca_loadings(adata, components=list(np.arange(1, n_pcs+1)), save='_' + str(sampleName) + '_' + str(n_pcs) + '_PCs')

#COMPUTING NEIGHBORHOOD GRAPH:
#Uses PCA representation of data matrix
###LOOK AT SAVED FIGURE W/ SUFFIX _PCA_VarianceRatio AND CHOOSE NUMBER OF PCs BEFORE APEX (HERE ~20)
n_neighbors = 25
min_dist = .1
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(adata, min_dist=min_dist)
sc.pl.umap(adata, color=labeled_genes, save='_' + str(sampleName) + '_' + str(n_neighbors) + 'neighbors_' + str(n_pcs) + 'PCs_umap_' + str(min_dist) + 'minDist_' + '_labeled_raw')
###THIS REQUIRES GENES IN LABELED_GENES TO BE IN HIGHLY VARIABLE LIST:
labeled_genes_var = ['Oprm1', 'Slc17a6', 'Grik3', 'Gad1', 'Gad2', 'Slc4a4', 'Cpa6', 'Ntsr2', 'Pdgfra', 'Luzp2', 'C1qc', 'Flt1', 'Pcbp3', 'Dbi']
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, save='_' + str(sampleName) + '_nNeighbors' + str(n_neighbors) + '_nPCs' + str(n_pcs) + '_umap_filtered')

###CLUSTER DATA WITH LOUVAIN ALGORITHM:
resolution = 0.5
sc.tl.louvain(adata, resolution=resolution)
labeled_genes.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_louvain')
labeled_genes_var.insert(0,'louvain')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_louvain_filtered')

###CLUSTER WITH LEIDEN ALGORITHM:
resolution = 0.5
iterations = -1 #only change this if you know what you're doing
sc.tl.leiden(adata, resolution=resolution, n_iterations=iterations)
#labeled_genes.insert(0,'leiden')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden')
#labeled_genes_var.insert(0,'leiden')
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered')
sc.pl.umap(adata, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution_clusters_labeled_leiden_only')


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

###SPLIT DATA BY GROUP TO EXAMINE CLUSTERS FOR EACH GROUP INDIVIDUALLY:
adata_g1 = adata[adata.obs['condition']==1]
adata_g2 = adata[adata.obs['condition']==2]
sc.pl.umap(adata_g1, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_filtered_' + str(g1n))
sc.pl.umap(adata_g2, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_filtered_' + str(g2n))
sc.pl.umap(adata_g1, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_only_' + str(g1n))
sc.pl.umap(adata_g2, color=['leiden'], use_raw=False, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_only' + str(g2n))

###RUN STATS, NOTE 't-test' can be changed to 'wilcoxon':
###COMPARE EXPRESSION BY CONDITION (EQUIVALENT TO BULK-SEQ):
adata.obs.condition=adata.obs.condition.astype('category')
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
method = 't-test' #t-test, wilcoxon, or logreg
cluster_method = 'leiden'
n_genes=1000 #set to adata.var.shape[0] to include all genes

sc.tl.rank_genes_groups(adata, cluster_method, use_raw=False, pts=True, n_genes=n_genes, method=method)
sc.pl.rank_genes_groups(adata, groupby=cluster_method, use_raw=False, n_genes=25, sharey=False, save='_' + str(sampleName) + '_' + method + '_' +  cluster_method + '_clusters')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(n_genes) #ALTER NUMBER TO SPECIFIC NUMBER OF GENES TO LIST
table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_' + cluster_method + '_cluster_table_' + str(n_genes) + 'genes.xlsx'), engine='openpyxl')
#make table with p-values included
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pval_table = pd.DataFrame(
        {group + '_' + key[:2]: result[key][group]
        for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join(BaseDirectory, sampleName + '_' + method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_filtered_corrected.xlsx'), engine='openpyxl')


###DIFFERENTIAL EXPRRESSION OF GENES WITHIN CLUSTERS:
pairs = list(zip(adata.obs['condition'], adata.obs['louvain'].astype('int')))
pairs_set = list(set(pairs))
s = sorted(pairs_set)
half = int((len(s)/2))
list1 = s[:half]
list2 = s[half:]
lz_louvain = list(zip(list1, list2))

pairs_leid = list(zip(adata.obs['condition'], adata.obs['leiden'].astype('int')))
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
n_genes = 1000

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
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[1])], reference=str(i[0]), use_raw=True, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, use_raw=True, method=method)        
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_rawData_adj.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1 USING RAW DATA: 
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_rawData_adj.xlsx'))

###CALCULATE GENES UPREGULATED IN GROUP 2 USING ONLY HIGHLY VARIABLE GENES:
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[1])], reference=str(i[0]), use_raw=False, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[1])], reference=str(i[0]), n_genes=n_genes, method=method)        
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g2n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_filteredHVGs_adj.xlsx'))
###CALCULATE GENES UPREGULATED IN GROUP 1 USING ONLY HIGHLY VARIABLE GENES: 
cat = pd.DataFrame()
for i in list2compare:
    if cluster_method=='leiden':
        sc.tl.rank_genes_groups(adata, 'pairs_' + cluster_method, groups=[str(i[0])], reference=str(i[1]), use_raw=False, n_genes=n_genes, method=method)
    elif cluster_method=='louvain':
        sc.tl.rank_genes_groups(adata, 'pairs', groups=[str(i[0])], reference=str(i[1]), n_genes=n_genes, method=method)
    #df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(500)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pval_table = pd.DataFrame(
            {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
    cat = pd.concat([cat, pval_table], axis=1)
    cat.to_excel(os.path.join(BaseDirectory, str(sampleName) + '_DiffExp_Upregulated' + str(g1n) + '_' + method + '_' + cluster_method + '_' + str(n_genes) + 'genes_filteredHVGs_adj.xlsx'))

###OPTIONALLY RECLUSTER A SINGLE CLUSTER INTO FURTHER SUBTYPES
cluster = 5
sc.tl.leiden(adata, restrict_to=('leiden', [str(cluster)]), resolution=resolution, n_iterations=iterations)
labeled_genes_var.insert(0, 'leiden_R')
labeled_genes.insert(0, 'leiden_R')
sc.pl.umap(adata, color=labeled_genes, wspace=0.5, save='_' + str(sampleName) + '_' + str(resolution) + 'resolution' + '_clusters_labeled_leiden_recluster' + str(cluster))
sc.pl.umap(adata, color=labeled_genes_var, use_raw=False, wspace=0.5, save='_' + str(sampleName) + str(resolution) + 'resolution_clusters_labeled_leiden_filtered_recluster' + str(cluster))

###PSEUDOTIME ANALYSIS (STILL IN DEVELOPMENT, CAN SKIP THIS)
adata.uns['iroot'] = np.flatnonzero(adata.obs['condition'] == 1)[0]
sc.tl.diffmap(adata, n_comps=10)
sc.tl.dpt(adata, n_branchings=1)
sc.pl.dpt_groups_pseudotime(adata, save='_pseudotime_groups')
#sc.pl.dpt_timeseries(adata, color_map=None, save='_pseudotime_timeSeries')


###PLOT AS HEATMAP:
t = table.T
markers = t[0].tolist()
sc.pl.heatmap(adata, var_names=markers, groupby='louvain', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_louvain')
sc.pl.heatmap(adata, var_names=markers, groupby='pairs', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_pairs')
sc.pl.heatmap(adata, var_names=markers, groupby='leiden', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden')

markers = t[18].tolist()
sc.pl.heatmap(adata, var_names=markers, groupby='leiden', show_gene_labels=True, save='_' + str(sampleName) + '_' + method + '_clusters_leiden_labeledGenes')

markers = ['Oprm1', 'Cmss1', 'Xylt1', "Camk2a", 'Grik3', 'Pcbp3', 'Grik1', 'Cux2', 'Rgs20', 'Luzp2', 'Hs6st3', 'Plp1', 'Pcbp3', 'Grm8', 'Inpp5d', 'Vcan', 'Arhgap6', 'Cpa6', 'Prex2', 'Flt1']

###HEIRARCHICAL CLUSTERING:
sc.tl.dendrogram(adata, groupby='louvain')
sc.pl.dendrogram(adata, groupby='louvain', save='_louvain')
sc.tl.dendrogram(adata, groupby='pairs')
sc.pl.dendrogram(adata, groupby='pairs', save='_pairs')
sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.dendrogram(adata, groupby='leiden', save='_leiden')

markers_dend = ['Oprm1', 'Slc17a6', 'Xylt1', "Camk2a", 'Grik3', 'Pcbp3', 'Grik1', 'Rgs6', 'Fstl4', 'Luzp2', 'Hs6st3', 'Plp1', 'Pcbp3', 'Grm8', 'Inpp5d', 'Vcan', 'Arhgap6', 'Cpa6', 'Prex2', 'Flt1']



###PLOT CORRELATION MATRIX:
sc.pl.correlation_matrix(adata, groupby='louvain', save='_louvain')
sc.pl.correlation_matrix(adata, groupby='leiden', save='_leiden')
sc.pl.correlation_matrix(adata, groupby='pairs', save='_pairs')
sc.pl.correlation_matrix(adata, groupby='condition', save='_condition')


###PLOT EMBEDDING DENSITY:
sc.tl.embedding_density(adata, basis='umap', groupby='condition', key_added='umap_density')
sc.pl.embedding_density(adata, basis='umap', key='umap_density', title='UMAP_Density', save='all')
sc.pl.embedding_density(adata, basis='umap', key='umap_density', title='UMAP_Density_Group1', group=[1], save='Group1')
sc.pl.embedding_density(adata, basis='umap', key='umap_density', title='UMAP_Density_Group2', group=[2], save='Group2')


###PLOT PARTITION-BASED GRAPH ABSTRACTIONS
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, save='_paga_louvain')
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, save='_paga_leiden')
sc.pl.paga_compare(adata, basis='umap', save='_paga_compare')

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
adata.write(results_file)

###IF YOU LATER WANT TO REVERT TO A RESULTS FILE USE:
adata = sc.read(results_file)



