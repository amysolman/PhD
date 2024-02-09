# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
#library(vegan)
#source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
library(dplyr)
library(cowplot)
library(reshape) 
library(funrar) #for make relative
library(stringr) #mutate function
library(gridExtra) #for exporting as pdf
library(scales)
# devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic")
library(tidyverse)
library(patchwork)
library(ggpubr)

#prokaryotes
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
#eukaryotes without micrometazoans
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 

#Community structure numbers

####################################################################################################################

#output results to text document
sink("../results/pma-community-profiles.txt", type="output")
writeLines("===============================================================
COMMUNITY PROFILES
===============================================================")
#Snow
writeLines("Percentage of reads per phylum total prokaryote community in snow:")
sub = subset_samples(pro, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in total snow community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the total snow community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

writeLines("Percentage of reads per phylum viable prokaryote community in snow:")
sub = subset_samples(pro, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in viable snow community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the viable snow community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

#Spring Ice
writeLines("Percentage of reads per phylum total prokaryote community in spring ice:")
sub = subset_samples(pro, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in total spring ice community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the total spring ice community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

writeLines("Percentage of reads per phylum viable prokaryote community in spring ice:")
sub = subset_samples(pro, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in viable spring ice community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the viable spring ice community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

#Summer Ice
writeLines("Percentage of reads per phylum total prokaryote community in summer ice:")
sub = subset_samples(pro, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in total summer ice community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the total summer ice community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

writeLines("Percentage of reads per phylum viable prokaryote community in summer ice:")
sub = subset_samples(pro, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in viable summer ice community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the viable summer ice community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

#Cryoconite
writeLines("Percentage of reads per phylum total prokaryote community in cryoconite:")
sub = subset_samples(pro, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in total cryoconite community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the total cryoconite community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

writeLines("Percentage of reads per phylum viable prokaryote community in cryoconite:")
sub = subset_samples(pro, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Archaea in viable cryoconite community:")
arch.com = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
snow.t.a = data.frame(tax_table(arch.com))[data.frame(tax_table(arch.com))$Domain == "Archaea",]
writeLines("There are") 
nrow(snow.t.a) 
writeLines("archaeal ASVs in the viable cryoconite community. These belong to:") 
rownames(snow.t.a) <- NULL
unique(snow.t.a[,c(2:6)])

#EUKARYOTES

#Snow
writeLines("Percentage of reads per phylum total eukaryote community in snow:")
sub = subset_samples(euk, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum viable eukaryote community in snow:")
sub = subset_samples(euk, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

#Spring Ice
writeLines("Percentage of reads per phylum total eukaryote community in spring ice:")
sub = subset_samples(euk, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum viable eukaryote community in spring ice:")
sub = subset_samples(euk, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

#Summer Ice
writeLines("Percentage of reads per phylum total eukaryote community in summer ice:")
sub = subset_samples(euk, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum viable eukaryote community in summer ice:")
sub = subset_samples(euk, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

#Cryoconite
writeLines("Percentage of reads per phylum total eukaryote community in cryoconite:")
sub = subset_samples(euk, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "tDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum viable eukaryote community in cryoconite:")
sub = subset_samples(euk, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "iDNA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

sink()