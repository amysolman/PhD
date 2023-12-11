#!/bin/bash

#prepare the data for analysis by making phyloseq objects 
Rscript data-prep.R

#normalise the data
Rscript repeated-rarefying.R

#MultiCoLA analysis
Rscript MultiCoLA.R

#create subcommunities
Rscript define-subcommunities.R

#explore community profiles
Rscript community-profiles.R

#unconstrained ordination PCoA, with Un-/Weighted Unifrac and Bray-Curtis Dissimilarities 
Rscript unconstrained-ordination-pcoa.R

#unconstrained ordination NMDS and Bray-Curtis Dissimilarities
Rscript unconstrained-ordination-nmds.R

#PERMANOVA for differences between community structure
Rscript PERMANOVA.R

#beta diversity partitioning analysis
Rscript beta-diversity-partitioning.R

#percentage of ASVs from each subcommunity in 25%, 75%, 99% of sites
Rscript percentage-in-sites.R

#correlations between subcommunity alpha diversity and environmental variables
Rscript alpha-diversity-correlations.R

#distance-decay relationship
Rscript DDR.R

#taxa-area relationship
Rscript TAR.R

#Mantel tests for single variables
Rscript mantel-single-variable.R

#Z-transform the environmental variables
Rscript z-transform.R

#Canonical correspondence analysis
Rscript CCA.R

#Distance-based redundancy analysis
Rscript db-RDA.R

#Variation partitioning analysis
Rscript VPA.R

#Mantel and Partial Mantel Tests
Rscript mantel-multiple-variables.R

#Neutral model - remember to make the inset plots manually using powerpoint!
Rscript neutral-model.R

#abundance-occupancy plot
Rscript abundance-occupancy.R

#Phylogenetic null model analysis - don't run this locally as it's too computationally expensive!
#Rscript null-model.R