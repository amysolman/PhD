#!/bin/bash

#analyse physio-chemical profiles of our sites
Rscript environmental-analysis-stats.R
Rscript environmental-analysis-plots.R

#prepare the data for analysis by making phyloseq objects 
Rscript data-prep.R

#explore control samples
Rscript controls-lib-size.R #library size of controls vs environmental samples
Rscript controls-profiles.R #what do our control profiles look like?
Rscript controls-high-abundance.R #which ASVs in controls are abundant in our environmental samples?
Rscript controls-remove.R #remove potentially contaminant ASVs identified from abundance and taxonomy
Rscript controls-NMDS.R #check for differences between communities before/after removing "contaminants" using NMDS and ANOVA

#normalise the data via proportional transformation and get basic community stats
Rscript normalise.R

#community composition plots
Rscript community-profiles.R

#uncostrained ordination (NMDS with Bray-Curtis Dissimilarities)
Rscript unconstrained-ordination.R

#hierarchical clustering analysis
Rscript hierarchical-clustering.R

#PERMANOVA
Rscript PERMANOVA.R

#Differential abundance analysis
Rscript differential-abundance.R

#Distance-decay analysis
Rscript DDR.R

#Explore differences in communities over different spatial scales
Rscript spatial-scales.R
