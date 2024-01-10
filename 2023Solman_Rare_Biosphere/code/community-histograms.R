# 1. Clear workspace and load packages

rm(list=ls())
graphics.off()

source("00-abundance-functions.R")
source("00-solman-functions.R")
library(dplyr) #for manipulating data where you see %>% and computing summary stats
#install.packages("naniar")
library(naniar) #for replace NA 
library(ggplot2)
library(ggpubr)
library(stringr)
library(phyloseq)
library(tidyr)
library(RColorBrewer)
#install.packages("ggvenn")             
library(ggvenn)
library(cowplot) #multiple plots
library(grid) #for removing white space
library(dplyr)
library(tibble)
#install.packages("betapart")
library(betapart)
library(reshape2)
library(vegan)
library(tidyverse)
#install.packages("glue")
library(glue)
#library(scales) #for percentages
# install.packages("indicspecies")
library(indicspecies)
library(ggrepel)
library(fantaxtic)
library(patchwork)
library(zoo) #load zoo package


# Method

#Using the appropriate abundance cut offs defined by MultiCoLA, the average rarefied datasets were divided into abundant, intermediate and rare subcommunities. Abundant porkaryote ASVs were those with relative abundance above `r pro_abun_cut*100`%. Rare prokaryote ASVs were those with relative abundance below `r pro_rare_cut*100`%. Abundant eukaryote ASVs were those with relative abundance above `r euk_abun_cut*100`%. Rare eukaryote ASVs were those with relative abundance below `r euk_rare_cut*100`%. Initial exploration of the subcommunties was carried out by calculating and visualizing phylum- and order-level relative abundances. Species richness and Shannon's diversity (to look at evenness) were calculated for each subcommunity within each glacier. Alpha diversities were plotted for subcommunitiy on each glacier. Bray-Curtis dissimilarities (chosen as it is not effected by high numbers of null values e.g. zeros) with NMDS ordination (chosen because it can be used with any dissimilarity metric, unlike PCA) was performed. Samples were normalised via proportional prior to beta-diversity analysis, as the Bray-Curtis dissimilarity metric assumes samples are of equal size. PCoA was also performed as this technique is believed to be more capable of dealing with a large number of zeros than NMDS (such as our rare communities). ANOSIM was carried out to test for significant differences in community composition between regions. To further explore the relative contributions of nestedness and turnover, β-diversity values were split into two components, namely the balanced variation (richness) and the unidirectional abundance gradients (turnover) within the distinct subcommunities, using the “bray.part” function in the “betapart” R package Balance variation in abundance (i.e. the turnover component) means to what extent differences in community are due to one species being replaced by another new species. Abundance gradients (i.e. the richness component) means to what extent differences in community are due to poorer sites being subsets of richer sites (i.e. individuals are lost from one site to another and no species is replaced by another species, it is simply missing). These are nested assemblages.  

# Results

# 3. Import data and separate into Arctic and Antarctic samples
# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds") 

#Arctic prokaryotes
p.ar.a <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
p.ar.i <- readRDS("../results/16S-phylo-object-arc-int.rds") 
p.ar.r <- readRDS("../results/16S-phylo-object-arc-rare.rds") 

#Antarctic prokaryotes
p.an.a <- readRDS("../results/16S-phylo-object-ant-abun.rds") 
p.an.i <- readRDS("../results/16S-phylo-object-ant-int.rds") 
p.an.r <- readRDS("../results/16S-phylo-object-ant-rare.rds") 

#Arctic eukaryotes
e.ar.a <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
e.ar.i <- readRDS("../results/18S-phylo-object-arc-int.rds") 
e.ar.r <- readRDS("../results/18S-phylo-object-arc-rare.rds") 

#Antarctic eukaryotes
e.an.a <- readRDS("../results/18S-phylo-object-ant-abun.rds") 
e.an.i <- readRDS("../results/18S-phylo-object-ant-int.rds") 
e.an.r <- readRDS("../results/18S-phylo-object-ant-rare.rds") 

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#How many reads to each of these subcommunities have? Do they have similar numbers of reads? Will they need to be normalised again?

#Arctic Abundant Prokaryotes
p.ar.a.df = data.table::data.table(as(sample_data(p.ar.a), "data.frame"),
                                   TotalReads = sample_sums(p.ar.a), keep.rownames = TRUE)

p.ar.a_p<-ggplot(data=p.ar.a.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.ar.a_p


#Arctic Intermediate  Prokaryotes - Very different library depth per sample!

p.ar.i.df = data.table::data.table(as(sample_data(p.ar.i), "data.frame"),
                                   TotalReads = sample_sums(p.ar.i), keep.rownames = TRUE)

p.ar.i_p<-ggplot(data=p.ar.i.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.ar.i_p

#Arctic Rare  Prokaryotes - Very different library depth per sample!

p.ar.r.df = data.table::data.table(as(sample_data(p.ar.r), "data.frame"),
                                   TotalReads = sample_sums(p.ar.r), keep.rownames = TRUE)

p.ar.r_p<-ggplot(data=p.ar.r.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.ar.r_p


#Antarctic Abundant Prokaryotes
p.an.a.df = data.table::data.table(as(sample_data(p.an.a), "data.frame"),
                                   TotalReads = sample_sums(p.an.a), keep.rownames = TRUE)

p.an.a_p<-ggplot(data=p.an.a.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.an.a_p


#Antarctic Intermediate  Prokaryotes - Very different library depth per sample!

p.an.i.df = data.table::data.table(as(sample_data(p.an.i), "data.frame"),
                                   TotalReads = sample_sums(p.an.i), keep.rownames = TRUE)

p.an.i_p<-ggplot(data=p.an.i.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.an.i_p

#Antarctic Rare  Prokaryotes - Very different library depth per sample!

p.an.r.df = data.table::data.table(as(sample_data(p.an.r), "data.frame"),
                                   TotalReads = sample_sums(p.an.r), keep.rownames = TRUE)

p.an.r_p<-ggplot(data=p.an.r.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

p.an.r_p


#Arctic Abundant Eukaryotes
e.ar.a.df = data.table::data.table(as(sample_data(e.ar.a), "data.frame"),
                                   TotalReads = sample_sums(e.ar.a), keep.rownames = TRUE)

e.ar.a_p<-ggplot(data=e.ar.a.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.ar.a_p


#Arctic Intermediate  Eukaryotes - Very different library depth per sample!

e.ar.i.df = data.table::data.table(as(sample_data(e.ar.i), "data.frame"),
                                   TotalReads = sample_sums(e.ar.i), keep.rownames = TRUE)

e.ar.i_p<-ggplot(data=e.ar.i.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.ar.i_p

#Arctic Rare  Eukaryotes - Very different library depth per sample!

e.ar.r.df = data.table::data.table(as(sample_data(e.ar.r), "data.frame"),
                                   TotalReads = sample_sums(e.ar.r), keep.rownames = TRUE)

e.ar.r_p<-ggplot(data=e.ar.r.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.ar.r_p


#Antarctic Abundant Eukaryotes
e.an.a.df = data.table::data.table(as(sample_data(e.an.a), "data.frame"),
                                   TotalReads = sample_sums(e.an.a), keep.rownames = TRUE)

e.an.a_p<-ggplot(data=e.an.a.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.an.a_p


#Antarctic Intermediate  Eukaryotes - Very different library depth per sample!

e.an.i.df = data.table::data.table(as(sample_data(e.an.i), "data.frame"),
                                   TotalReads = sample_sums(e.an.i), keep.rownames = TRUE)

e.an.i_p<-ggplot(data=e.an.i.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.an.i_p

#Antarctic Rare  Eukaryotes - Very different library depth per sample!

e.an.r.df = data.table::data.table(as(sample_data(e.an.r), "data.frame"),
                                   TotalReads = sample_sums(e.an.r), keep.rownames = TRUE)

e.an.r_p<-ggplot(data=e.an.r.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

e.an.r_p