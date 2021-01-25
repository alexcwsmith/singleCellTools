#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:17:11 2021

@author: smith
"""
import numpy as np
import pandas as pd
import os

def countDEGs(file, directory, pcutoff=.05, plot=True, save=False):
    """Count number of differentially expressed genes in scanpy result file.
    
    Parameters
    ----------
    file : string
        Path to saved .xlsx or .csv file containing differential expression data.
    directory : string
        Directory to save results.
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
    elif file.endswith('.csv'):
        df = pd.read_csv(file, index_col=0)
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
        fig = res.plot(kind='bar')
        ax = fig.get_figure()
        ax.savefig(os.path.join(directory, fname + '_DEG_Counts.png'))
    if save:
        res.to_excel(os.path.join(directory, fname + '_DEG_Counts.xlsx'))
    return res
    
