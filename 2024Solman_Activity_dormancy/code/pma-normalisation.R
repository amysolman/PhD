#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)

#2. Import data
pro <- readRDS("../results/pma-16S-phylo-object-no-controls.rds")
euk <- readRDS("../results/pma-18S-phylo-object-no-controls.rds")

x = data.frame(tax_table(euk))

#Normalise data

sort(sample_sums(pro), decreasing = FALSE)
sort(sample_sums(euk), decreasing = FALSE)

#as the minimum number of reads per environmental sample is 2300+ I'm not going to remove low abundance reads
#I'm going to use proportional transformation to normalise all data
pro.rel = transform_sample_counts(pro, function(x) x/sum(x))
#sanity check
colSums(data.frame(otu_table(pro.rel)))
euk.rel = transform_sample_counts(euk, function(x) x/sum(x))
#sanity check
colSums(data.frame(otu_table(euk.rel)))

#alternative method: Rarefy data
pro.rare = rarefy_even_depth(pro, sample.size = min(sample_sums(pro)))
euk.rare = rarefy_even_depth(euk, sample.size = min(sample_sums(euk)))

#save phyloseq objects
saveRDS(pro.rare, "../results/pma-16S-phylo-object-norm-rare.rds")
saveRDS(euk.rare, "../results/pma-18S-phylo-object-norm-rare.rds")
saveRDS(pro.rel, "../results/pma-16S-phylo-object-norm-rel.rds")
saveRDS(euk.rel, "../results/pma-18S-phylo-object-norm-rel.rds")

####################################################################################################################

#output results to text document
sink("../results/pma-asv-read-numbers.txt", type="output")
writeLines("===============================================================
NUMBER OF ASVS AND READS IN OUR PHYLOSEQ OBJECTS
===============================================================")
writeLines("Prokaryote dataset rarefied to:")
min(sample_sums(pro))
writeLines("Eukaryote dataset rarefied to:")
min(sample_sums(euk))
writeLines("Number of ASVs in prokaryote dataset before rarefying:")
ntaxa(pro)
writeLines("Number of ASVs in prokaryote dataset after rarefying:")
ntaxa(pro.rare)
writeLines("Number of reads in prokaryote dataset before rarefying:")
sum(sample_sums(pro))
writeLines("Number of reads in prokaryote dataset after rarefying:")
sum(sample_sums(pro.rare))
writeLines("Number of ASVs in eukaryote dataset before rarefying:")
ntaxa(euk)
writeLines("Number of ASVs in eukaryote dataset after rarefying:")
ntaxa(euk.rare)
writeLines("Number of reads in eukaryote dataset before rarefying:")
sum(sample_sums(euk))
writeLines("Number of reads in eukaryote dataset after rarefying:")
sum(sample_sums(euk.rare))
sink()
