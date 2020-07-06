#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 20:57:34 2020

From example template at:
    http://localhost:8888/notebooks/goatools/notebooks/goea_nbt3102.ipynb

@author: smith
"""
from __future__ import print_function
import pandas as pd
import numpy as np
import os
from Bio import Entrez
import ACWS_GOAtools as go

"""
Before starting anything in this script, open ACWS_GOAtools.py and change variables:
    goaResultDir - change to path you want to save results to
    genesDf - change to path to .xlsx file of differential expression data from scanPy
    comparison - change to descriptor of data you are analyzing, example 'UpregulatedSaline'
"""

comparison = 'DiffExp_UpOC'

#Run for a single cluster
go.runGOanalysis('Cluster3')

#Make a list of clusters and run for all. Change integer to number of clusters.
clusters = ['Cluster' + str(x) for x in range(20)]
for cluster in clusters:
    go.runGOanalysis(cluster)


GOcategory = 'MF' #Can be 'MF' for molecular  functions, 'BP' for biogical processes, or 'CC' for cellular components
comparison = 'DiffExp_UpOC'
go.combineGOresults(clusters, GOcategory, comparison)



