#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:09:58 2020

@author: smith
"""
import os
import pandas as pd
import numpy as np
import collections as cx
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.godag_plot import plot_gos, plot_results
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from genes_ncbi_10090_proteincoding import GENEID2NT as GeneID2nt_mus
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from Bio import Entrez
Entrez.email = 'alexander.smith@mssm.edu'
Entrez.api_key = '8f015d17da5ac2df68b84d558b46e5150508'

#Change this to true before first import:
download = False

comparison = 'DiffExp_DownOC'
goaResultDir = '/home/smith/Smith_Scripts/NLP_GeneExpression/spaCy_model062920/Results/GOA_Results/'
genesDf = pd.read_excel('/d1/studies/cellranger/ACWS_DP/scanpy_DiffExp_V2/results_maxGenes3000_maxMito.05_MinDisp0.2/DP_OC_Saline_Merged_DiffExp_UpregulatedSaline_t-test.xlsx', index_col=0)


if download:
    obo_fname = download_go_basic_obo()
    fin_gene2go = download_ncbi_associations()

if not download:
    obo_fname = '/d1/studies/singleCellTools/go-basic.obo'
    fin_gene2go = '/d1/studies/singleCellTools/gene2go'
    
obodag = GODag("go-basic.obo")
objanno = Gene2GoReader(fin_gene2go, taxids=[10090])
ns2assoc = objanno.get_ns2assc()
for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))
goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_mus.keys(), # List of mouse protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

indexDirectory = '/home/smith/Smith_Scripts/AllenSDK/CircuitFinder/Indexes/'
geneIndex = pd.read_excel(os.path.join(indexDirectory, 'ISH_experiment_gene_index.xlsx'), index_col=0)

def getEntrezIDs(genes, setName):
    ID = []
    IDs = []
    nullGenes = []
    for gene in genes:
        try:
            gene_info = geneIndex.loc[geneIndex.acronym==str(gene)]
            entrez_id = gene_info['entrez_id'].astype(int)
            eid = entrez_id.unique()
            if len(eid) == 0:
                nullGenes.append(gene)
            ID.append(eid.tolist())
        except ValueError:
            nullGenes.append(gene)
            pass
    for sublist in ID:
        for item in sublist:
            IDs.append(item)
    return(IDs, nullGenes)

def runGOanalysis(cluster):
    clusterNum=cluster.replace('Cluster', '')
    genesList = []
    if comparison == 'DiffExp_DownOC':
        sigGenes = genesDf.loc[genesDf['(1, ' + str(clusterNum) + ')_p'] < .05]
        genesList = sigGenes['(1, ' + str(clusterNum) + ')_n'].tolist()
    elif comparison == 'DiffExp_UpOC':
        sigGenes = genesDf.loc[genesDf['(2, ' + str(clusterNum) + ')_p'] < .05]
        genesList = sigGenes['(2, ' + str(clusterNum) + ')_n'].tolist()        
    elif comparison =='MarkerGenes':
        genesList = genesDf[str(clusterNum) + '_n']

    setName = cluster + '_' + comparison
    entrezIDlist, nullGenes = getEntrezIDs(genesList, setName)
    
    print(str(len(entrezIDlist)) + " Entrez IDs retrieved, searching NCBI for " + str(len(nullGenes)) + " genes")
    
    newEntrezList = []
    unfoundList = []
    for gene in nullGenes:
        search = Entrez.esearch(db='gene', term=gene, retmax=1, start=0, sort='relevance')
        record = Entrez.read(search)
        ID = record.get('IdList')
        if not len(ID) == 0:
            newEntrezList.append(ID[0])
        else:
            unfoundList.append(gene)
            
    print(str(len(newEntrezList)) + " additional Entrez IDs retrieved, " + str(len(unfoundList)) + " genes with no results")
    
    entrezIDlist.extend(newEntrezList)
    
    print("Proceeding with GO analysis of " + str(len(entrezIDlist)) + " genes")
    
    goea_results_all = goeaobj.run_study(entrezIDlist)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    
    print('{N} of {M:,} results were significant'.format(
        N=len(goea_results_sig),
        M=len(goea_results_all)))
    
    print('Significant results: {E} enriched, {P} purified'.format(
        E=sum(1 for r in goea_results_sig if r.enrichment=='e'),
        P=sum(1 for r in goea_results_sig if r.enrichment=='p')))
    
    ctr = cx.Counter([r.NS for r in goea_results_sig])
    print('Significant results[{TOTAL}] = {BP} BP + {MF} MF + {CC} CC'.format(
        TOTAL=len(goea_results_all),
        BP=ctr['BP'],  # biological_process
        MF=ctr['MF'],  # molecular_function
        CC=ctr['CC'])) # cellular_component
    
    
    goeaobj.wr_xlsx(os.path.join(goaResultDir, cluster + '_' + comparison + "_GOA_Results.xlsx"), goea_results_sig)
    goeaobj.wr_txt(os.path.join(goaResultDir, cluster + '_' + comparison + "_GOA_Results.txt"), goea_results_sig)
    
    plot_results(os.path.join(goaResultDir, cluster + '_' + comparison + "_GOA.png"), goea_results_sig)
    
    goid_subset = [
        'GO:0098984', # CC neuron to neuron synapse
        'GO:0099699', # CC integral component of synaptic membrane
        'GO:2000474', # BP regulation of opioid receptor signaling
        'GO:0038003', # BP opioid receptor signaling pathway
        'GO:0048167', # BP regulation of synaptic plasticity
        'GO:0021897', # BP forebrain astrocyte development
        'GO:0014003', # BP oligodendroyte development
        'GO:0035249', # BP synaptic transmission, glutamatergic
        'GO:0051932', # BP synaptic transmission, GABAergic
        'GO:0014004', # BP microglia differentiation
    ]
    
    
    """
    This plot contains GOEA results:
    
        GO terms colored by P-value:
            pval < 0.005 (light red)
            pval < 0.01 (light orange)
            pval < 0.05 (yellow)
            pval > 0.05 (grey) Study terms that are not statistically significant
        GO terms with study gene counts printed. e.g., "32 genes"
    
    """
    
    
    plot_gos(os.path.join(goaResultDir, cluster + '_' + comparison + "_GOA_Subset_Plot.png"), 
        goid_subset, # Source GO ids
        obodag,
        goea_results=goea_results_all, # use pvals for coloring
        # We can further configure the plot...
        id2symbol=geneid2symbol, # Print study gene Symbols, not Entrez GeneIDs
        study_items=6, # Only only 6 gene Symbols max on GO terms
        items_p_line=3, # Print 3 genes per line
        )

def combineGOresults(clusters, GOcategory, comparison):
    """

    Parameters
    ----------
    clusters : list
        List of clusters to combine
    GOcategory : string
        GO namespace to combine. Options are 'BP', 'MF', 'CC'
    comparison : string
        Descriptor for results being compared, for reading/writing result files.
        

    Returns
    -------
    None.

    """
    catdf1 = pd.DataFrame()
    for cluster in clusters:
        resDf = pd.read_excel(os.path.join(goaResultDir, cluster + '_' + comparison + '_GOA_Results.xlsx'), index_col=0)
        res = resDf.loc[resDf['NS']==GOcategory]
        res = res.reset_index()
        res = res.iloc[:,[0,3,9]]
        res['Cluster'] = cluster
        catdf1 = pd.concat([catdf1, res])
    catdf1.to_excel(os.path.join(goaResultDir, GOcategory + '_GOA_Enriched_' + comparison + '.xlsx'))
    names = catdf1.name.unique().tolist()
    
    df2 = pd.DataFrame()
    catdf2 = pd.DataFrame()
    for process in names:
        resdf = catdf1.loc[catdf1['name']==process]
        resdf = resdf.reset_index()
        count = resdf.shape[0]
        clus = resdf['Cluster'].unique().tolist()
        df2['GO'] = resdf['GO']
        df2['name'] = process
        df2['count'] = count
        df2['ave_corrected_pval'] = np.mean(resdf['p_fdr_bh'])
        df2['clusters'] = pd.Series([clus])
        catdf2 = pd.concat([catdf2, df2], axis=0)   
    catdf2 = catdf2.dropna(how='any', subset=['clusters'])
    catdf2 = catdf2.sort_values(by='count', ascending=False)
    catdf2.to_excel(os.path.join(goaResultDir, GOcategory + '_GOA_Enriched_' + comparison + '_Combined.xlsx'))




