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

#unconstrained ordination PCoA
Rscript unconstrained-ordination-pcoa.R

#unconstrained ordination NMDS
#Rscript unconstrained-ordination-nmds.R