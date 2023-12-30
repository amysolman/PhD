
# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
# library(RColorBrewer)
# library(tidyr) #trans data between wide and long format
# library(ggplot2)
library(vegan)
# source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
# library(dplyr)
# library(cowplot)
# library(reshape) 
# library(funrar) #for make relative
# library(stringr) #mutate function, to split columns
# library(gridExtra) #for exporting as pdf
# library(scales)
library(microbiome) #for summarize_phyloseq

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-controls-removed.rds")
#eukaryotes without micrometazoans
ps.euk.rm <- readRDS("../results/18S-ps-no-mm-controls-removed.rds")
#eukaryotes with micrometazoans
ps.euk <- readRDS("../results/18S-ps-mm-keep-controls-removed.rds")
#micrometazoans only
ps.mm <- readRDS("../results/18S-ps-mm-only-controls-removed.rds")

#Let's plot rarefaction curves to explore where to remove samples

asv.tab.pro = data.frame(t(otu_table(ps.pro)))
pdf("../results/prokaryote-rarefaction-curve.pdf")
rarecurve(asv.tab.pro, step=100, cex=0.5)
dev.off()
sort(sample_sums(ps.pro), decreasing = FALSE) #remove samples with less than 1900 reads

asv.tab.euk = data.frame(t(otu_table(ps.euk.rm)))
pdf("../results/total-eukaryote-rarefaction-curve.pdf")
rarecurve(asv.tab.euk, step=100, cex=0.5)
dev.off()
sort(sample_sums(ps.euk.rm), decreasing = FALSE) #remove samples with less than 1900 reads

asv.tab.euk.full = data.frame(t(otu_table(ps.euk)))
pdf("../results/micro-eukaryote-rarefaction-curve.pdf")
rarecurve(asv.tab.euk.full, step=100, cex=0.5)
dev.off()
sort(sample_sums(ps.euk), decreasing = FALSE) #remove samples with less than 1900 reads

asv.tab.mm = data.frame(t(otu_table(ps.mm)))
pdf("../results/microfauna-rarefaction-curve.pdf")
rarecurve(asv.tab.mm, step=100, cex=0.5)
dev.off()
sort(sample_sums(ps.mm), decreasing = FALSE) #remove samples with less than 500 reads


remove_proportion <- function(ps, num){
  
  sort(sample_sums(ps))
  
  #retain samples with >= num counts
  ps.p = prune_samples(sample_sums(ps)>= num, ps)
  
  #remove ASVs with zero counts
  ps.f = filter_taxa(ps.p, function(x) sum(x) > 0, TRUE)
  
  #proportional transformation
  ps.ra <- transform_sample_counts(ps.f, function(x) x / sum(x) )
  
  return(ps.ra)
  
}

pro.p <- remove_proportion(ps.pro, 1900)
euk.p <- remove_proportion(ps.euk.rm, 1900)
euk.p.full <- remove_proportion(ps.euk, 1900)
mm.p <- remove_proportion(ps.mm, 500) #this at 500 because it's a subset of the euk community

saveRDS(pro.p, "../results/16S-ps-norm.rds")
saveRDS(euk.p, "../results/18S-no-mm-ps-norm.rds")
saveRDS(euk.p.full, "../results/18S-mm-included-ps-norm.rds")
saveRDS(mm.p, "../results/18S-mm-only-ps-norm.rds")

###################################################################################################
###################################################################################################

#Report normalisation Results
sink("../results/normalisation.txt", type="output")
writeLines("===============================================================
NORMALISATION RESULTS
===============================================================")
writeLines("After removing samples with low read depth...") 
writeLines("Number of prokaryote samples remaining:") 
nsamples(pro.p)
writeLines("Number of prokaryote ASVs remaining:") 
ntaxa(pro.p)
writeLines("Number of total eukaryote samples remaining:") 
nsamples(euk.p.full)
writeLines("Number of total eukaryote ASVs remaining:") 
ntaxa(euk.p.full)
writeLines("Number of micro-eukaryote samples remaining:") 
nsamples(euk.p)
writeLines("Number of micro-eukaryote ASVs remaining:") 
ntaxa(euk.p)
writeLines("Number of micro-fauna samples remaining:") 
nsamples(mm.p)
writeLines("Number of micro-fauna ASVs remaining:") 
ntaxa(mm.p)
sink()
###################################################################################################
###################################################################################################
