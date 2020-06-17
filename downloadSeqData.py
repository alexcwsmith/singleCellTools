#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:53:45 2020

@author: smith
"""

import urllib
import os

def downloadSeqData():
    sampleName = input("Enter Sample ID: ")
    baseUrl = input("Enter directory to download from (should end in /outs/): " ) #'https://wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01155_StephanieCaligiuri_PaulKenney/710_3GEX_N1_Chp_3rd/outs/'
    path = input("Enter save directory: " ) #'/d1/studies/cellranger/ACWS_DP/SeqDownloader/'
    p = os.path.expanduser(path)
    flist = ['molecule_info.h5', 'filtered_feature_bc_matrix.h5', 'metrics_summary.csv']
    for file in flist:
        urllib.request.urlretrieve(os.path.join(baseUrl, file), filename=os.path.join(p, sampleName + '_' + file))
             
    dlist = ['matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz']
    url2 = os.path.join(baseUrl, 'filtered_feature_bc_matrix/')
    p = os.path.join(p, 'filtered_feature_bc_matrix')
    if not os.path.exists(p):
        os.mkdir(p)
    for file in dlist:
        urllib.request.urlretrieve(os.path.join(url2, file), filename=os.path.join(p, sampleName + '_' + file))
    

downloadSeqData()



