#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 16:01:27 2021

@author: smith
"""
import numpy as np
import pandas as pd
import scanpy as sc

#adata=sc.read('/home/smith/FinalAnalysis.h5ad')

def findCellsByGeneCoex(adata, gene1, gene2=None, thresh=0.6, gene1up=True, gene2up=True):
    """Find cells based on co-expression of two genes. Currently only works if both are in "highly variable genes" list.
    
    Parameters
    ----------
    adata : AnnData object
        Dataset to filter
    gene1 : string
        First gene to filter by.
    gene2 : string, optional (default None)
        Second gene to filter by.
    thresh : float, optional (default 0.6)
        Expression theshold to filter. Be careful to know the data type of your input (e.g. has the data been log transformed?)
    gene1Up : bool, optional (default True)
        Whether you want to filter gene1 above threshold (set False to keep cells below threshold)
    gene2Up : bool, optional (default True)
        Whether you want to filter gene2 above threshold (set False to keep cells below threshold)
        
    Returns
    -------
    Filtered dataframe.
    
    Examples
    --------
    result = findCellsByGeneCoex(adata, 'Sst', 'Gad1', 0.6, True, True)
    result = findCellsByGeneCoex(adata, gene1='Slc17a6', gene2='Emx1', thresh=1.5, gene1up=False, gene2up=True)
    """
    mtx = pd.DataFrame(adata.X)
    mtx.columns = adata.var.index
    mtx.index = adata.obs.index
    try:
        if gene1up:
            df = mtx.loc[mtx[gene1]>thresh]
        elif not gene1up:
            df = mtx.loc[mtx[gene1]<thresh]
    except KeyError:
        raise KeyError(str(gene1) + " not found, is it included in highly variable genes?")
    try:
        if gene2:
            if gene2up:
                df1 = df.loc[df[gene2]>thresh]
            elif not gene2up:
                df1 = df.loc[df[gene2]<thresh]
        elif not gene2:
            print(str(df.shape[0] + " cells stored with mean " + str(gene1) + "expression" + str(round(np.mean(df[gene1],4)))))
            return df
    except KeyError:
        raise KeyError(str(gene2) + " not found. Is it included in highly variable genes")
    print(str(df1.shape[0]) + " cells stored with mean " + str(gene1) + " expression " + str(round(np.mean(df1[gene1]), 4)) + 
          " and mean " + str(gene2) + " expression " + str(round(np.mean(df1[gene2]), 4)))
    return df1



#result = findCellsByGeneCoex(adata, 'Sst', 'Gad1', 0.6, gene1up=True, gene2up=True)
