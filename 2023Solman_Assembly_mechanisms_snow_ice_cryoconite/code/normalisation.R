
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
  
  return(list(ps.f, ps.ra))
  
}

pro.p <- remove_proportion(ps.pro, 1900)
euk.p <- remove_proportion(ps.euk.rm, 1900)
euk.p.full <- remove_proportion(ps.euk, 1900)
mm.p <- remove_proportion(ps.mm, 500) #this at 500 because it's a subset of the euk community

saveRDS(pro.p[[2]], "../results/16S-ps-norm.rds")
saveRDS(euk.p[[2]], "../results/18S-no-mm-ps-norm.rds")
saveRDS(euk.p.full[[2]], "../results/18S-mm-included-ps-norm.rds")
saveRDS(mm.p[[2]], "../results/18S-mm-only-ps-norm.rds")

###################################################################################################
###################################################################################################

#Report normalisation Results
sink("../results/normalisation.txt", type="output")
writeLines("===============================================================
NORMALISATION RESULTS
===============================================================")
writeLines("After removing samples with low read depth...") 
writeLines("Number of prokaryote samples remaining:") 
nsamples(pro.p[[2]])
writeLines("Number of prokaryote ASVs remaining:") 
ntaxa(pro.p[[2]])
writeLines("Number of total eukaryote samples remaining:") 
nsamples(euk.p.full[[2]])
writeLines("Number of total eukaryote ASVs remaining:") 
ntaxa(euk.p.full[[2]])
writeLines("Number of micro-eukaryote samples remaining:") 
nsamples(euk.p[[2]])
writeLines("Number of micro-eukaryote ASVs remaining:") 
ntaxa(euk.p[[2]])
writeLines("Number of micro-fauna samples remaining:") 
nsamples(mm.p[[2]])
writeLines("Number of micro-fauna ASVs remaining:") 
ntaxa(mm.p[[2]])
sink()
###################################################################################################
###################################################################################################

####################################################################################################################

#output results to text document
sink("../results/community-profiles.txt", type="output")
writeLines("===============================================================
COMMUNITY PROFILES
===============================================================")

writeLines("Percentage of reads per phylum in snow prokaryote subcommunity:")
sub = subset_samples(pro.p[[1]], Habitat == "Snow")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in spring ice prokaryote subcommunity:")
sub = subset_samples(pro.p[[1]], Habitat == "Spring Ice")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in summer ice prokaryote subcommunity:")
sub = subset_samples(pro.p[[1]], Habitat == "Summer Ice")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in cryoconite prokaryote subcommunity:")
sub = subset_samples(pro.p[[1]], Habitat == "Cryoconite")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("BREAK DOWN OF GENERA IN CRYOCONITE:")
sub = subset_samples(pro.p[[1]], Habitat == "Cryoconite")
glom <- tax_glom(sub, taxrank = 'Genus')
t.glom = data.frame(otu_table(glom))
x = data.frame(tax_table(glom))$Genus
x2 = paste0(x, 1:length(x))
rownames(t.glom) = x2
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in snow eukaryote subcommunity:")
sub = subset_samples(euk.p[[1]], Habitat == "Snow")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in spring ice eukaryote subcommunity:")
sub = subset_samples(euk.p[[1]], Habitat == "Spring Ice")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in summer ice eukaryote subcommunity:")
sub = subset_samples(euk.p[[1]], Habitat == "Summer Ice")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in cryoconite eukaryote subcommunity:")
sub = subset_samples(euk.p[[1]], Habitat == "Cryoconite")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

#do we have any archaea in our datasets?
#remember microfauna have alread been removed

#Archaea in our datasets

#Snow
sub = subset_samples(pro.p[[1]], Habitat == "Snow") #subset to snow
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
xi = data.frame(tax_table(sub2))
xia = xi[xi$Domain == "Archaea",]

#Spring Ice
sub = subset_samples(pro.p[[1]], Habitat == "Spring Ice") #subset to spring ice
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
xii = data.frame(tax_table(sub2))
xiia = xii[xii$Domain == "Archaea",]

#Summer Ice
sub = subset_samples(pro.p[[1]], Habitat == "Summer Ice") #subset to summer ice
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
xiii = data.frame(tax_table(sub2))
xiiia = xiii[xiii$Domain == "Archaea",]

#Cryoconite
sub = subset_samples(pro.p[[1]], Habitat == "Cryoconite") #subset to cryoconite
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
xiiii = data.frame(tax_table(sub2))
xiiiia = xiiii[xiiii$Domain == "Archaea",]

writeLines("There are") 
nrow(xia) 
writeLines("archaeal ASVs in snow.") 
xia
writeLines("There are") 
nrow(xiia) 
writeLines("archaeal ASVs in spring ice.") 
xiia
writeLines("There are") 
nrow(xiiia) 
writeLines("archaeal ASVs in summer ice.") 
xiiia
writeLines("There are") 
nrow(xiiii) 
writeLines("archaeal ASVs in cryoconite.") 
xiiiia

writeLines("MICROFAUNA") 

#Snow
sub = subset_samples(mm.p[[2]], Habitat == "Snow") #subset to snow
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
yi = data.frame(tax_table(sub2))

writeLines("Number of microfauna ASVs in snow") 
nrow(yi)

writeLines("Taxonomy of microfauna ASVs in snow") 
yi

#Spring Ice
sub = subset_samples(mm.p[[2]], Habitat == "Spring Ice") #subset
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
yii = data.frame(tax_table(sub2))

writeLines("Number of microfauna ASVs in spring ice") 
nrow(yii)

writeLines("Taxonomy of microfauna ASVs in spring ice") 
yii

#Summer Ice
sub = subset_samples(mm.p[[2]], Habitat == "Summer Ice") #subset
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
yiii = data.frame(tax_table(sub2))

writeLines("Number of microfauna ASVs in summer ice") 
nrow(yiii)

writeLines("Taxonomy of microfauna ASVs in summer ice") 
yiii

#Cryoconite
sub = subset_samples(mm.p[[2]], Habitat == "Cryoconite") #subset
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
yiiii = data.frame(tax_table(sub2))

writeLines("Number of microfauna ASVs in cryoconite") 
nrow(yiiii)

writeLines("Taxonomy of microfauna ASVs in cryoconite") 
yiiii

  
  writeLines("OVERALL TAXONOMIC DIVERSITY") 
  
#Snow
sub = subset_samples(pro.p[[1]], Habitat == "Snow")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Prokaryote Phyla in snow") 
length(unique(x$Phylum))

writeLines("Prokaryote Genera in snow") 
length(unique(x$Genus))

#Spring Ice
sub = subset_samples(pro.p[[1]], Habitat == "Spring Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Prokaryote Phyla in spring ice") 
length(unique(x$Phylum))

writeLines("Prokaryote Genera in spring ice") 
length(unique(x$Genus))

#Summer Ice
sub = subset_samples(pro.p[[1]], Habitat == "Summer Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Prokaryote Phyla in summer ice") 
length(unique(x$Phylum))

writeLines("Prokaryote Genera in summer ice") 
length(unique(x$Genus))

#Cryoconite
sub = subset_samples(pro.p[[1]], Habitat == "Cryoconite")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Prokaryote Phyla in cryoconite") 
length(unique(x$Phylum))

writeLines("Prokaryote Genera in cryoconite") 
length(unique(x$Genus))

#Snow
sub = subset_samples(euk.p[[1]], Habitat == "Snow")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Eukaryote Phyla in snow") 
length(unique(x$Phylum))

writeLines("Eukaryote Genera in snow") 
length(unique(x$Genus))

#Spring Ice
sub = subset_samples(euk.p[[1]], Habitat == "Spring Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Eukaryote Phyla in spring ice") 
length(unique(x$Phylum))

writeLines("Eukaryote Genera in spring ice") 
length(unique(x$Genus))

#Summer Ice
sub = subset_samples(euk.p[[1]], Habitat == "Summer Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Eukaryote Phyla in summer ice") 
length(unique(x$Phylum))

writeLines("Eukaryote Genera in summer ice") 
length(unique(x$Genus))

#Cryoconite
sub = subset_samples(euk.p[[1]], Habitat == "Cryoconite")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Eukaryote Phyla in cryoconite") 
length(unique(x$Phylum))

writeLines("Eukaryote Genera in cryoconite") 
length(unique(x$Genus))


#Snow
sub = subset_samples(mm.p[[1]], Habitat == "Snow")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Microfauna Phyla in snow") 
length(unique(x$Phylum))

writeLines("Microfauna Genera in snow") 
length(unique(x$Genus))

#Spring Ice
sub = subset_samples(mm.p[[1]], Habitat == "Spring Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Microfauna Phyla in spring ice") 
length(unique(x$Phylum))

writeLines("Microfauna Genera in spring ice") 
length(unique(x$Genus))

#Summer Ice
sub = subset_samples(mm.p[[1]], Habitat == "Summer Ice")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Microfauna Phyla in summer ice") 
length(unique(x$Phylum))

writeLines("Microfauna Genera in summer ice") 
length(unique(x$Genus))

#Cryoconite
sub = subset_samples(mm.p[[1]], Habitat == "Cryoconite")
sub2 = filter_taxa(sub, function(x) sum(x) > 0, TRUE) #remove ASVs with no counts
x = data.frame(tax_table(sub2))

writeLines("Microfauna Phyla in cryoconite") 
length(unique(x$Phylum))

writeLines("Microfauna Genera in cryoconite") 
length(unique(x$Genus))
sink()

