#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:53:45 2020

@author: smith

Script to download 10X Genomics sequencing data off of Minerva @ Mount Sinai


"""

import urllib.request
import os

def downloadSpatialSeqData():
    sampleName = input("Enter Sample ID: ")
    baseUrl = input("Enter directory to download from (should end in /outs/): " ) #'https://wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01155_StephanieCaligiuri_PaulKenney/710_3GEX_N1_Chp_3rd/outs/'
    path = input("Enter save directory: " ) #'/d1/studies/cellranger/ACWS_DP/SeqDownloader/'
    p = os.path.expanduser(path)
    flist = ['molecule_info.h5', 'filtered_feature_bc_matrix.h5', 'metrics_summary.csv']
    for file in flist:
        urllib.request.urlretrieve(os.path.join(baseUrl, file), filename=os.path.join(p, sampleName + '_' + file))
             
    dlist = ['matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz']
    url2 = os.path.join(baseUrl, 'filtered_feature_bc_matrix/')
    p2 = os.path.join(p, 'filtered_feature_bc_matrix')
    if not os.path.exists(p):
        os.mkdir(p2)
    for file in dlist:
        urllib.request.urlretrieve(os.path.join(url2, file), filename=os.path.join(p2, sampleName + '_' + file))
    imList = ['aligned_fiducials.jpg', 'detected_tissue_image.jpg', 'scalefactors_json.json', 'tissue_hires_image.png', 'tissue_lowres_image.png', 'tissue_positions_list.csv']
    spatialUrl = os.path.join(baseUrl, 'spatial/')
    imDir = os.path.join(p, 'images/')
    if not os.path.exists(imDir):
        os.mkdir(imDir)
    for file in imList:
        try:
            urllib.request.urlretrieve(os.path.join(spatialUrl, file), filename=os.path.join(imDir, sampleName + '_' + file))
        except:
            print(file + ' not found')

if __name__ == '__main__':
    downloadSpatialSeqData()



