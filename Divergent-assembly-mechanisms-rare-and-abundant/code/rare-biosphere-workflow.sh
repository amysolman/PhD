#!/bin/bash

#prepare the data for analysis by making phyloseq objects 
Rscript data-prep.R

#normalise the data
#Rscript repreated-rarefying.R

#MultiCoLA analysis
#Rscript MultiCoLA.R

#create subcommunities
#Rscript define-subcommunities.R