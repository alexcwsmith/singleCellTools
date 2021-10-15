# singleCellTools
Tools for analysis of single-cell sequencing data, by Alex Smith in the lab of Paul Kenny @ Mount Sinai.

[![DOI](https://zenodo.org/badge/273064116.svg)](https://zenodo.org/badge/latestdoi/273064116)

#Overview:

* ACWS_scanPy_MASTER.py - contains a full workflow for analysis of single-cell sequencing data. Currently uses scanpy 1.7 and AnnData 0.7.5

* ACWS_filterCells.py - contains functions for filtering and extracting cells from adata objects based on gene co-expression.
This can be good for identifying cells of interest, or for filtering of "cells" that express high levels of multiple genes that
are believed to be mutually exclusive, and thus are likely dead cells or ambient RNA.

* ACWS_GOAtools.py & goaTools_MASTER.py are deprecated now that scanpy has the integrated gprofiler ontology queries, but I will
leave the scripts here.

* countDEGs.py contains a function for counting # of differentially expressed genes in each cluster.

* downloadSeqData.py & downloadFASTQs.py are for retrieving data from Mount Sinai's Minerva HPC.

## Cellranger tools:

* cellrangerCount.sh - shell script for running cellranger count. Users will need to make their own copy, this file is read-only for everyone but root.

* cellrangerAggr.sh - shell script for running cellranger aggr. Again, make a copy to edit.

## CellBender tools:

* cellbender_conda_env.yml - conda environment for running cellbender, a package for de-noising data by removing ambient RNA that I find quite useful.
  Details on cellbender are here: https://github.com/broadinstitute/CellBender

* cellbender.sh - shell script for running cellbender (inside activated conda env above).
