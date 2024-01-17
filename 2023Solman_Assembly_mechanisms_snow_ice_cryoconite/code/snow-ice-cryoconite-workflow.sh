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

#Explore differences in communities and environment over different spatial scales
Rscript spatial-scales.R

#Single variable mantel tests
Rscript mantel-single-var.R

#CCA and db-RDA 
Rscript constrained-ordination.R

#mantel and db-RDA plots were combined manually using powerpoint because I struggled to do this with code

#VPA and (partial) Mantel Tests
Rscript vpa-mantel.R

#Mantel correlograms - don't run this locally as it's too computationally expensive! Run on HPC.
#Rscript mantel-correlogram.R #full script
#model for each group to make it run quicker
#Rscript mantel-correlogram-pro.R 
#Rscript mantel-correlogram-euk.R 
#Rscript mantel-correlogram-mm.R 

#Mantel correlogram plots - plot the data from the HPC
Rscript mantel-correlogram-plots.R

#Null model - don't run this locally as it's too computationally expensive! Run on HPC.
#Rscript null-model.R #full model
#model in parts to make it quicker to run
#Rscript null-model-pro-bnti.R
#Rscript null-model-euk-bnti.R
#Rscript null-model-pro-rcbray.R
#Rscript null-model-euk-rcbray.R

#Neutral model and null model plot
Rscript neutral-model-with-null-model-plot.R

#Network analysis
Rscript network-analysis.R

#Network properties - stats for the whole networks
Rscript network-properties.R

#Node properties - stats for ASVs in the networks
Rscript node-properties.R