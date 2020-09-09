#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 20:57:34 2020


@author: smith
"""

import ACWS_GOAtools as go

"""
Before starting anything in this script, open ACWS_GOAtools.py and change variables:
    goaResultDir - change to path you want to save results to
    genesDf - change to path to .xlsx file of differential expression data from scanPy
    comparison - change to descriptor of data you are analyzing, example 'UpregulatedSaline'
"""

comparison = 'Clem_Oxy'

#Run for a single cluster
go._runGOanalysis('Cluster8', n_genes=75)

#Make a list of clusters and run for all. Change integer to number of clusters.
clusters = ['Cluster' + str(x) for x in range(16)]
newGeneIndex = go.runGOanalysis(clusters, processes=8)

#Combine individual cluster data files into one large data file
GOcategory = 'CC' #Can be 'MF' for molecular  functions, 'BP' for biogical processes, or 'CC' for cellular components
go.combineGOresults(clusters, GOcategory, comparison)



