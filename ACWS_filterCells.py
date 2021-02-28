#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 16:01:27 2021

@author: smith
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy

def findCellsByGeneCoex(adata, gene1, gene2=None, g1thresh=0.6, g2thresh=0.6, gene1up=True, gene2up=True, use_raw=False):
    """Find cells based on co-expression of two genes.
    
    Parameters
    ----------
    adata : AnnData object
        Dataset to filter
    gene1 : string
        First gene to filter by.
    gene2 : string, optional (default None)
        Second gene to filter by.
    g1thresh : float, optional (default 0.6)
        Expression theshold for gene1. Be careful to know the data type of your input (e.g. has the data been log transformed?)
    g2thresh : float, optional (default 0.6)
        Expression theshold for gene1. Be careful to know the data type of your input (e.g. has the data been log transformed?)       
    gene1Up : bool, optional (default True)
        Whether you want to filter gene1 above threshold (set False to keep cells below threshold)
    gene2Up : bool, optional (default True)
        Whether you want to filter gene2 above threshold (set False to keep cells below threshold)
    use_raw : bool, optional (default False)
        Whether to use the adata.X matrix containing only highly variable genes (if True), or adata.raw.X 
        Note that using raw data takes significantly longer than only HVGs so should only be used when necessary.
        
    Returns
    -------
    Filtered dataframe.
    
    Examples
    --------
    >>> result = findCellsByGeneCoex(adata, 'Sst', 'Gad1', 0.6, 1.0, True, True, False)
    >>> result = findCellsByGeneCoex(adata, gene1='Slc17a6', gene2='Emx1', g1thresh=1.5, g2thresh=0.6, gene1up=False, gene2up=True)
    """
    if type(adata.X)==scipy.sparse.csr.csr_matrix:
        if use_raw:
            mtx = pd.DataFrame(adata.raw.X.toarray())
            mtx.columns=adata.raw.var_names
            mtx.index=adata.raw.obs_names
        elif not use_raw:
            mtx = pd.DataFrame(adata.X.toarray())
            mtx.columns = adata.var.index
            mtx.index = adata.obs.index
    elif type(adata.X)==np.ndarray:
        if use_raw:
            mtx = pd.DataFrame(adata.raw.X.toarray())
            mtx.columns=adata.raw.var_names
            mtx.index=adata.raw.obs_names            
        elif not use_raw:
            mtx = pd.DataFrame(adata.X)
            mtx.columns = adata.var.index
            mtx.index = adata.obs.index
    else:
        raise TypeError(str(type(adata.X)) + " is not a valid data type. Must be a scipy sparse matrix or numpy ndarray.")
    try:
        if gene1up:
            df = mtx.loc[mtx[gene1]>=g1thresh]
        elif not gene1up:
            df = mtx.loc[mtx[gene1]<g1thresh]
    except KeyError:
        if not use_raw:
            raise KeyError(str(gene1) + " not found, you may need to set use_raw=True.")
        elif use_raw:
            raise KeyError(str(gene1) + " not found in dataset")
    try:
        if gene2:
            if gene2up:
                df1 = df.loc[df[gene2]>=g2thresh]
            elif not gene2up:
                df1 = df.loc[df[gene2]<g2thresh]
        elif not gene2:
            print(str(df.shape[0]) + " cells stored with mean " + str(gene1) + "expression" + str(np.mean(df[gene1])))
            return df
    except KeyError:
        if not use_raw:
            raise KeyError(str(gene2) + " not found, you may need to set use_raw=True.")
        elif use_raw:
            raise KeyError(str(gene2) + " not found in dataset")
    print(str(df1.shape[0]) + " cells stored with mean " + str(gene1) + " expression " + str(round(np.mean(df1[gene1]), 4)) + 
          " and mean " + str(gene2) + " expression " + str(round(np.mean(df1[gene2]), 4)))
    return df1


def filterCellsByCoex(adata, df):
    """Removes cells from adata if they are contained in df.
    
    Parameters
    ----------
    adata : AnnData object.
        Data to filter
    df : pandas DataFrame
        Dataframe of cells to remove (result of findCellsByGeneCoex)
    """
    cells = adata.obs.index.tolist()
    exCells = df.index.tolist()
    filteredCells = [c for c in cells if c not in exCells]
    return(adata[adata.obs.index.isin(filteredCells)].copy())
    
def subsampleData(adata, df):
    """Returns a full AnnData object containing only cells identified in findCellsByGeneCoex.
    
    Parameters
    ----------
    adata : AnnData object.
        The scanPy AnnData object results were originally extracted from.
    df : Pandas DataFrame.
        The result dataframe from findCellsByGeneCoex
        
    Example:
    >>> result = findCellsByGeneCoex(adata, 'Gad1', 'Sst', 0.6, gene1up=True, gene2up=True, use_raw=True)
    >>> adata=subsampleAdata(result,adata)

    """
    return(adata[adata.obs.index.isin(df.index.tolist())])
    
def countDEGs(file, directory, n_genes=1000, pcutoff=.05, plot=True, save=False):
    """Count number of differentially expressed genes in scanpy result file.
    
    Parameters
    ----------
    file : string
        Path to saved .xlsx or .csv file containing differential expression data.
    directory : string
        Directory to save results.
    n_genes : int, (optional, default 1000)
        Number of genes used in original data analysis.
    pcutoff : float (optional, default .05)
        Alpha value for significance.
    save : bool (optional, default False)
        Whether to save or only return result.
        
    Returns
    -------
    Pandas DataFrame with # of DEGs for each cluster.
    """
    fname, ext = os.path.splitext(os.path.basename(file))
    if file.endswith('.xlsx'):
        df = pd.read_excel(file, index_col=0, engine='openpyxl')
        df = df[:n_genes]
    elif file.endswith('.csv'):
        df = pd.read_csv(file, index_col=0)
        df = df[:n_genes]
    clusters=[]
    degs=[]
    for col in df.columns:
        if col.endswith('_p'):
            count = (df[col]<pcutoff).value_counts()
            try:
                count = count.loc[count.index==True].values[0]
            except IndexError:
                count=0
            clu = int(col.strip('_p').split(' ')[-1].strip(')'))
            clusters.append(clu)
            degs.append(count)
    lz = list(zip(clusters,degs))
    res = pd.DataFrame(lz)
    res.columns=['Cluster', 'DEGs']
    res.set_index('Cluster', inplace=True, drop=True)
    if plot:
        fig = res.plot(kind='bar', grid=False)
        ax = fig.get_figure()
        ax.savefig(os.path.join(directory, 'figures/' + fname + '_DEG_Counts.png'))
    if save:
        res.to_excel(os.path.join(directory, fname + '_DEG_Counts.xlsx'))
    return res

def meanMito():
    meanMito = adata.obs['percent_mito'].mean()
    return(meanMito)
    
