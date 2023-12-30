#!/bin/bash

#PMA Workflow

#process the 16S and 18S amplicons via command line 
#qiime2-pma

#prep the amplicon data for analysis and create phyloseq objects
Rscript pma-data-prep.R

#explore and remove contaminants
Rscript pma-contamination.R

#normalise data
Rscript pma-normalisation.R

#get community profiles
Rscript pma-community-profiles.R

#prep the qPCR data for analysis
#python PMA-qPCR.py

#plot and test for significant differences
#Rscript qPCR-plot.R

#METATRANSCRIPTOMICS workflow

#process the 16S and 18S amplicons via command line 
#qiime2-metatranscriptomics 

#prep the amplicon data for analysis and create phyloseq objects
Rscript mt-dna-data-prep.R