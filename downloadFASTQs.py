#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 21:10:43 2021

@author: smith
"""
import requests
import urllib
import os
from bs4 import BeautifulSoup

def downloadFASTQs():
    directory = input("Enter URL to FASTQ directory: ")
    saveDir = input("Enter directory to save to: ")
    prefix = input("Enter prefix of samples to download: ")
    r = requests.get(directory)
    soup = BeautifulSoup(r.text)
    table = soup.findAll('table')[0]
    links = table.findAll('a')
    sampleNames=[]
    for link in links:
        if link.text.startswith(prefix):
            sampleNames.append(link.text)
    DLpaths = []
    for sample in sampleNames:
        url = os.path.join(directory, sample)
        sample = sample.strip('.fastq.gz')
        DLpaths.append(url)
    
    for path in DLpaths:
        name = path.split('/')[-1]
        urllib.request.urlretrieve(path, filename=os.path.join(saveDir, name))

    
if __name__ == '__main__':
    downloadFASTQs()









