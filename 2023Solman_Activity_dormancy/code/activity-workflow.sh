#!/bin/bash

#PMA Workflow

#process the 16S and 18S amplicons via qiime2 in command line 
#qiime2-pma

#prep the amplicon data for analysis and create phyloseq objects
Rscript pma-data-prep.R

#explore and remove contaminants
Rscript pma-contamination.R

#normalise data
Rscript pma-normalisation.R

#get community profiles
Rscript pma-community-profiles.R

#get figures about community structure
Rscript pma-community-figures.R

#alpha diversity analysis
Rscript pma-alpha-diversity.R

#beta diversity analysis with unconstrained ordination
Rscript pma-unconstrained-ordination.R

#test for significant differences in structure between viable and total communities
Rscript pma-PERMANOVA.R

#ASV-level viability ratios
Rscript pma-idna-tdna-ratios.R

#phylum-level viability ratios
Rscript pma-phylum-idna-tdna-ratios.R

#prep the qPCR data for analysis
python pma-qPCR-data-prep.py

#plot and test for significant differences
Rscript pma-qPCR-plot.R

#pma qPCR inhibition test
Rscript pma-qPCR-inhibition-test.R

#pma protocol test analysis script
Rscript pma-qPCR-test.R

#METATRANSCRIPTOMICS workflow

#process the 16S and 18S amplicons via command line 
#qiime2-metatranscriptomics 

#prep the amplicon data for analysis and create phyloseq objects
Rscript mt-dna-data-prep.R