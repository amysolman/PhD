# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Define the abundance cut offs
# 3. Import data
# 4. Use the previously written abundance functions to assign each ASV to an abundance category and create new phyloseq objects 
# 5. Export the phyloseq objects
# 6. Define function to merge subcommunities into one phyloseq object
# 7. Merge subcommunities
# 8. Define relative abundances function
# 9. Plot Phylum level abundances
# 9. Plot Class level abundances
# 10. Plot Order level abundances.
# 12. Visualise PCoA plots of based on Bray-Curtis dissimilarities between subcommunities
# 13. Test for differences in community composition between groups (glaciers and poles).

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

# What do the different communities look like in terms of taxa?

#what are the proportions of each community that are rare, abundant and intermediate?

#Prokaryote Antarctic 

####get points for putting labels
pie.df1 = data.frame(Num = c(ntaxa(p.an.r), ntaxa(p.an.i), ntaxa(p.an.a)),
                     Group = c("Rare", "Intermediate", "Abundant"))
pie.df1$value = round(pie.df1$Num/sum(pie.df1$Num)*100, 1)

pie.df1b <- pie.df1 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df1b <- pie.df1b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie1 = ggplot(pie.df1, aes(x = "" , y = value, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_fill_brewer(palette="Set1")+
  geom_label_repel(data = pie.df1b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 15, nudge_x = 1, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie1

#Prokaryote Arctic 

  ####get points for putting labels
  pie.df2 = data.frame(Num = c(ntaxa(p.ar.r), ntaxa(p.ar.i), ntaxa(p.ar.a)),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df2$value = round(pie.df2$Num/sum(pie.df2$Num)*100, 1)
  
  pie.df2b <- pie.df2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
  pie.df2b <- pie.df2b %>%
         mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie2 = ggplot(pie.df2, aes(x = "" , y = value, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_fill_brewer(palette="Set1")+
  geom_label_repel(data = pie.df2b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 15, nudge_x = 1, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
  theme_void()
  pie2
  
  #Eukaryote Antarctic 
  
  ####get points for putting labels
  pie.df3 = data.frame(Num = c(ntaxa(e.an.r), ntaxa(e.an.i), ntaxa(e.an.a)),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df3$value = round(pie.df3$Num/sum(pie.df3$Num)*100, 1)
  
  pie.df3b <- pie.df3 %>% 
    mutate(csum = rev(cumsum(rev(value))), 
           pos = value/2 + lead(csum, 1))
  pie.df3b <- pie.df3b %>%
    mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))
  
  pie3 = ggplot(pie.df3, aes(x = "" , y = value, fill = fct_inorder(Group))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    #scale_fill_brewer(palette = "Pastel1") +
    scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    #scale_fill_brewer(palette="Set1")+
    geom_label_repel(data = pie.df3b,
                     aes(y = pos, label = paste0(value, "%")),
                     size = 15, nudge_x = 1, show.legend = FALSE) +
    guides(fill = "none") +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
    theme_void()
  pie3
  
    #Eukaryote Arctic 
  
    ####get points for putting labels
  pie.df4 = data.frame(Num = c(ntaxa(e.ar.r), ntaxa(e.ar.i), ntaxa(e.ar.a)),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df4$value = round(pie.df4$Num/sum(pie.df4$Num)*100, 1)
  
  pie.df4b <- pie.df4 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
  pie.df4b <- pie.df4b %>%
         mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))
  
  pie4 = ggplot(pie.df4, aes(x = "" , y = value, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_fill_brewer(palette="Set1")+
  geom_label_repel(data = pie.df4b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 15, nudge_x = 1, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
  theme_void()
  pie4

# What do the different communities look like in terms of number of reads?

#what are the proportions of each community that are rare, abundant and intermediate?
  
  #Prokaryote Antarctic 
  
  ####get points for putting labels
  pie.df5 = data.frame(Num = c(sum(sample_sums(p.an.r)), sum(sample_sums(p.an.i)), sum(sample_sums((p.an.a)))),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df5$value = round(pie.df5$Num/sum(pie.df5$Num)*100, 1)
  
  pie.df5b <- pie.df5 %>% 
    mutate(csum = rev(cumsum(rev(value))), 
           pos = value/2 + lead(csum, 1))
  pie.df5b <- pie.df5b %>%
    mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))
  
  pie5 = ggplot(pie.df5, aes(x = "" , y = value, fill = fct_inorder(Group))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    #scale_fill_brewer(palette = "Pastel1") +
    scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    #scale_fill_brewer(palette="Set1")+
    geom_label_repel(data = pie.df5b,
                     aes(y = pos, label = paste0(value, "%")),
                     size = 15, nudge_x = 1, show.legend = FALSE) +
    guides(fill = "none") +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
    theme_void()
  pie5

#Prokaryote Arctic 

  ####get points for putting labels
  pie.df6 = data.frame(Num = c(sum(sample_sums(p.ar.r)), sum(sample_sums(p.ar.i)), sum(sample_sums((p.ar.a)))),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df6$value = round(pie.df6$Num/sum(pie.df6$Num)*100, 1)
  
  pie.df6b <- pie.df6 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
  pie.df6b <- pie.df6b %>%
         mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie6 = ggplot(pie.df6, aes(x = "" , y = value, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_fill_brewer(palette="Set1")+
  geom_label_repel(data = pie.df6b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 15, nudge_x = 1, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
  theme_void()
  pie6
  
  #Eukaryote Antarctic 
  
  ####get points for putting labels
  pie.df7 = data.frame(Num = c(sum(sample_sums(e.an.r)), sum(sample_sums(e.an.i)), sum(sample_sums((e.an.a)))),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df7$value = round(pie.df7$Num/sum(pie.df7$Num)*100, 1)
  
  pie.df7b <- pie.df7 %>% 
    mutate(csum = rev(cumsum(rev(value))), 
           pos = value/2 + lead(csum, 1))
  pie.df7b <- pie.df7b %>%
    mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))
  
  pie7 = ggplot(pie.df7, aes(x = "" , y = value, fill = fct_inorder(Group))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    #scale_fill_brewer(palette = "Pastel1") +
    scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    #scale_fill_brewer(palette="Set1")+
    geom_label_repel(data = pie.df7b,
                     aes(y = pos, label = paste0(value, "%")),
                     size = 15, nudge_x = 1, show.legend = FALSE) +
    guides(fill = "none") +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
    theme_void()
  pie7
  
    #Eukaryote Arctic 
  
    ####get points for putting labels
  pie.df8 = data.frame(Num = c(sum(sample_sums(e.ar.r)), sum(sample_sums(e.ar.i)), sum(sample_sums((e.ar.a)))),
                       Group = c("Rare", "Intermediate", "Abundant"))
  pie.df8$value = round(pie.df8$Num/sum(pie.df8$Num)*100, 1)
  
  pie.df8b <- pie.df8 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
  pie.df8b <- pie.df8b %>%
         mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))
  
  pie8 = ggplot(pie.df8, aes(x = "" , y = value, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_fill_brewer(palette="Set1")+
  geom_label_repel(data = pie.df8b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 15, nudge_x = 1, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank())+
  theme_void()
  pie8
  

# 6. Define function to merge subcommunities into one phyloseq object

#Function to seperate ASVs into subcommunity samples and merge as phyloseq object. Requires original phyloseq object and list of ASV IDs and their subcommunities.

#What we want to end up with is pro_merge object with 3 versions of each sample (e.g. 150A, 150M, 150R)

# phylo = pro #full phylo object
# res = pro_res #subcommunities results

subcommunity_merge <- function(phylo, res){
  
    #data check
    #divide into Arctic and Antarctic samples
    ps_arc <- subset_samples(phylo, Pole=="Arctic")
    ps_arc = filter_taxa(ps_arc, function(x) sum(x) >= 1, TRUE)
    ps_ant <- subset_samples(phylo, Pole=="Antarctic")
    ps_ant = filter_taxa(ps_ant, function(x) sum(x) >= 1, TRUE)

    #do our subsets capture the whole community
    ## Get vectors of numbered OTUs/ASVs
    treatment1 <- colnames(t(otu_table(ps_arc)))
    treatment2 <- colnames(t(otu_table(ps_ant)))

    ## Get the intersection
   shared <- intersect(treatment1, treatment2)

    #Should equal TRUE
    ntaxa(res[[1]]) + ntaxa(res[[2]]) + ntaxa(res[[3]]) + ntaxa(res[[4]]) + ntaxa(res[[5]]) + ntaxa(res[[6]]) - length(shared) == ntaxa(phylo)
  
    #change samples names to reflect subcommunity
    sample_names(res[[1]]) = paste0(sample_names(res[[1]]), "A")
    sample_names(res[[2]]) = paste0(sample_names(res[[2]]), "M")
    sample_names(res[[3]]) = paste0(sample_names(res[[3]]), "R")
    sample_names(res[[4]]) = paste0(sample_names(res[[4]]), "A")
    sample_names(res[[5]]) = paste0(sample_names(res[[5]]), "M")
    sample_names(res[[6]]) = paste0(sample_names(res[[6]]), "R")
    
    #extract ASV tables and merge into one count table
    x1 = data.frame(otu_table(res[[1]]), check.names=FALSE)
    x2 = data.frame(otu_table(res[[2]]), check.names=FALSE)
    x3 = data.frame(otu_table(res[[3]]), check.names=FALSE)
    x4 = data.frame(otu_table(res[[4]]), check.names=FALSE)
    x5 = data.frame(otu_table(res[[5]]), check.names=FALSE)
    x6 = data.frame(otu_table(res[[6]]), check.names=FALSE)
    
    #create empty dataframe with ASVs as rows and sample names as columns
    nams = c(colnames(x1), colnames(x2), colnames(x3), colnames(x4), colnames(x5), colnames(x6))
    row.nams = taxa_names(phylo)
    final.df <- data.frame(matrix(ncol = length(nams), nrow = length(row.nams)))
    colnames(final.df) = nams
    row.names(final.df) = row.nams
    final.df[is.na(final.df)] <- 0
    
    #add in the data 
    #make rownames a column
    library(tibble)
    final.df <- tibble::rownames_to_column(final.df, "ASV")
    x1 <- tibble::rownames_to_column(x1, "ASV")
    x2 <- tibble::rownames_to_column(x2, "ASV")
    x3 <- tibble::rownames_to_column(x3, "ASV")
    x4 <- tibble::rownames_to_column(x4, "ASV")
    x5 <- tibble::rownames_to_column(x5, "ASV")
    x6 <- tibble::rownames_to_column(x6, "ASV")


  final.df.plus = final.df %>% 
  rows_update(x1, by = "ASV")
  
  final.df.plus = final.df.plus %>% 
  rows_update(x2, by = "ASV")
    final.df.plus = final.df.plus %>% 
  rows_update(x3, by = "ASV")
      final.df.plus = final.df.plus %>% 
  rows_update(x4, by = "ASV")
        final.df.plus = final.df.plus %>% 
  rows_update(x5, by = "ASV")
          final.df.plus = final.df.plus %>% 
  rows_update(x6, by = "ASV")
          
    final.df.plus2 <- final.df.plus[,-1]
    rownames(final.df.plus2) <- final.df.plus[,1]
    
    sum(colSums(final.df.plus2))
    
    # METADATA
    #prep metadata table
    samp.data1 = data.frame(sample_data(res[[1]]))
    samp.data1$SampleID = rownames(samp.data1)
    samp.data1$Subcommunity = "Abundant"
    samp.data2 = data.frame(sample_data(res[[2]]))
    samp.data2$SampleID = rownames(samp.data2)
    samp.data2$Subcommunity = "Intermediate"
    samp.data3 = data.frame(sample_data(res[[3]]))
    samp.data3$SampleID = rownames(samp.data3)
    samp.data3$Subcommunity = "Rare"
    samp.data4 = data.frame(sample_data(res[[4]]))
    samp.data4$SampleID = rownames(samp.data4)
    samp.data4$Subcommunity = "Abundant"
    samp.data5 = data.frame(sample_data(res[[5]]))
    samp.data5$SampleID = rownames(samp.data5)
    samp.data5$Subcommunity = "Intermediate"
    samp.data6 = data.frame(sample_data(res[[6]]))
    samp.data6$SampleID = rownames(samp.data6)
    samp.data6$Subcommunity = "Rare"
    
    final.samp.data = rbind(samp.data1, samp.data2, samp.data3, samp.data4, samp.data5, samp.data6)
    
    new.phylo = phyloseq(otu_table(final.df.plus2, taxa_are_rows = TRUE), sample_data(final.samp.data), tax_table(phylo), phy_tree(phylo))
    
    #remove samples and taxa with zero counts
    new.phylo.edit = filter_taxa(new.phylo, function(x) sum(x) > 0, TRUE)
    new.phylo.edit = prune_samples(sample_sums(new.phylo.edit) > 0, new.phylo.edit)
    
    new.phylo
    new.phylo.edit #they're the same

#check it has worked
# #how many ASVs in the Arctic dataset
#divide into Arctic and Antarctic samples
ps_arc.new <- subset_samples(new.phylo, Pole=="Arctic")
ps_arc.new = filter_taxa(ps_arc.new, function(x) sum(x) >= 1, TRUE)
ps_ant.new <- subset_samples(new.phylo, Pole=="Antarctic")
ps_ant.new = filter_taxa(ps_ant.new, function(x) sum(x) >= 1, TRUE)

ps_arc
ps_arc.new
ps_ant
ps_ant.new

return(new.phylo)
}


# 7. Merge subcommunities
pro_res = list(p.ar.a, p.ar.i, p.ar.r, p.an.a, p.an.i, p.an.r)
pro_sub_merge = subcommunity_merge(pro, pro_res)
ntaxa(pro_sub_merge) == ntaxa(pro) #should be TRUE - does the new objects have the same num ASVs as the original?

euk_res = list(e.ar.a, e.ar.i, e.ar.r, e.an.a, e.an.i, e.an.r)
euk_sub_merge = subcommunity_merge(euk, euk_res)
ntaxa(euk_sub_merge) == ntaxa(euk) #should be TRUE - does the new objects have the same num ASVs as the original?

#save merged phyloseq objects
saveRDS(pro_sub_merge, "../results/16S-phylo-object-sub-coms-merged.rds")
saveRDS(euk_sub_merge, "../results/18S-phylo-object-sub-coms-merged.rds")


#PLOT PROKARYOTE COMMUNITIES

#STEP ONE: find the 7 phyla with the highest mean relative abundance for each subcommunity in each pole (ignoring NA phyla)
top.phy <- top_taxa(pro_sub_merge, 
                    tax_level = "Phylum", 
                    n_taxa = 7,
                    grouping = c("Pole", "Subcommunity"))
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
pro.m = pro_sub_merge
x = data.frame(tax_table(pro.m)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))

#replace anything that says NA with Genus Unknown
y = x
y$Phylum <- gsub("P.NA", "(Genus Unknown)", y$Phylum)
y$Class <- gsub("P.NA C.NA", "(Genus Unknown)", y$Class)
y$Class <- gsub("C.NA", "(Genus Unknown)", y$Class)
y$Order <- gsub("P.NA C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("O.NA", "(Genus Unknown)", y$Order)
y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("F.NA", "(Genus Unknown)", y$Family)
y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)

tax_table(pro.m) = as.matrix(y)


#STEP THREE: find the three most abundant genus within each phyla by agglomerating data at the genus level then finding the genera with the highest mean abundance in each phylum

#get our data
data <-
  pro.m %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum - we'll use this dataframe to assign colours
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = pro.cols[pro.cols$Phylum == mini.df$Phylum,]$Colour1
  col2 = pro.cols[pro.cols$Phylum == mini.df$Phylum,]$Colour2
  
  #get colour function for that phylum
  cols.fun = colorRampPalette(c(col1, col2))
  
  #get colours for each class
  save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#add phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
#number of all phyla in our dataframe
length(unique(data2plot$Phylum))

#number of all genera in our dataframe
length(unique(data2plot$Genus))

#remove unwanted genera and replace with other
# for (i in 1:nrow(data2plot)){
#   if(!data2plot$Genus[i] %in% out$Genus){ 
#     #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
#     data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
#   } else if (data2plot$Genus[i] %in% out$Genus){
#     #else they are just given the phyla name with the genus name
#     data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
#   }
# }

#combine genus and phylum names for those we will keep
#out$Genus2 = paste0(out$Phylum, ": ", out$Genus)

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Genus2 = paste0(data2plot$Phylum, ": ", data2plot$Genus)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Genus2[i] %in% out.col$Genus){
    #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Genus2[i] %in% out.col$Genus){ 
    #else they are just given the phyla name with the genus name
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    
  }
}

#sanity check
#are these genera the same as those in my keep dataframe? Assuming Proteobacteria are in our keep phyla
unique(data2plot[data2plot$Phylum == "Proteobacteria",]$Genus)
df2keep[df2keep$Phylum == "Proteobacteria",]$Genus

length(unique(data2plot$Phylum)) # - our number of phyla stay the same
length(unique(data2plot$Genus)) # - our number of genus' has decreased because many low abundance genera have become "other"
sort(unique(data2plot$Genus))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "genus" with Other
    data2plot$Genus[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Genus))
unique(data2plot$Genus)
sort(unique(data2plot$Genus))

#add black to genera classified as "Other"
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))

#add colours to plotting dataframe
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
#set data as factors to keep our plotting in the right order
df$Subcommunity <- factor(df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))
df$Pole <- factor(df$Pole, levels = c("Antarctic", "Arctic"))
df$SubPole <- paste0(df$Pole, " (", df$Subcommunity, ")")
df$SubPole <- factor(df$SubPole, levels = c("Antarctic (Rare)", "Antarctic (Intermediate)", "Antarctic (Abundant)",
                                            "Arctic (Rare)", "Arctic (Intermediate)", "Arctic (Abundant)"))

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Genus as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Genus)), ]
levs = gen.df$Genus
levs = levs[levs != "Other"]
df.plot$Genus <- factor(df.plot$Genus, levels = out.col2$Genus)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Genus)

df.plot2 = 
  df.plot %>%
  group_by(Genus, SubPole, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p1 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol = 5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p1

pdf("../results/prokaryote-cryoconite-per-sample-community-plot.pdf", width=35, height=10)
print(p1)
dev.off()

df.plot3 = 
  df.plot %>%
  group_by(Genus, SubPole, Glacier) %>%
  dplyr::summarise(across(c(Abundance), sum))


p2 = ggplot(df.plot3, aes(fill=Genus, y=Abundance, x=Glacier)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(nrow = 8, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p2

pdf("../results/prokaryote-cryoconite-per-glacier-community-plot.pdf", width=35, height=10)
print(p2)
dev.off()

df.plot4 = 
  df.plot %>%
  group_by(Genus, SubPole) %>%
  dplyr::summarise(across(c(Abundance), sum))

p3 = ggplot(df.plot4, aes(fill=Genus, y=Abundance, x=SubPole)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(nrow=8, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p3

pdf("../results/prokaryote-cryoconite-per-subcom-community-plot.pdf", width=20, height=10)
print(p3)
dev.off()

#PLOT EUKARYOTE COMMUNITIES

#STEP ONE: find the 10 phyla with the highest mean relative abundance for each subcommunity in each pole (ignoring NA phyla)
top.phy <- top_taxa(euk_sub_merge,
                    tax_level = "Phylum",
                    n_taxa = 10,
                    grouping = c("Pole", "Subcommunity"))
# 
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep

#KEEP ALL PHYLA AS THERE ARE ONLY 10 IN THIS DATASET
# x = data.frame(tax_table(euk_sub_merge))
# length(unique(x$Phylum)) #10 unique phylum
# phy.keep = unique(x$Phylum)
# phy.keep = phy.keep[!is.na(phy.keep)] #remove NA
# phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
euk.m = euk_sub_merge
x = data.frame(tax_table(euk.m)) %>%
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))

#replace anything that says NA with Genus Unknown
y = x
y$Phylum <- gsub("P.NA", "(Genus Unknown)", y$Phylum)
y$Class <- gsub("P.NA C.NA", "(Genus Unknown)", y$Class)
y$Class <- gsub("C.NA", "(Genus Unknown)", y$Class)
y$Order <- gsub("P.NA C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("O.NA", "(Genus Unknown)", y$Order)
y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("F.NA", "(Genus Unknown)", y$Family)
y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)

tax_table(euk.m) = as.matrix(y)

#STEP THREE: find the three most abundant genus within each phyla by agglomerating data at the genus level then finding the genera with the highest mean abundance in each phylum

#get our data
data <-
  euk.m %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum - we'll use this dataframe to assign colours 
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)


#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = euk.cols[euk.cols$Phylum == mini.df$Phylum,]$Colour1
  col2 = euk.cols[euk.cols$Phylum == mini.df$Phylum,]$Colour2
  
  #get colour function for that phylum
  cols.fun = colorRampPalette(c(col1, col2))
  
  #get colours for each class
  save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#add phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
#number of all phyla in our dataframe
length(unique(data2plot$Phylum))
#number of all genera in our dataframe
length(unique(data2plot$Genus))

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Genus2 = paste0(data2plot$Phylum, ": ", data2plot$Genus)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Genus2[i] %in% out.col$Genus){
    #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Genus2[i] %in% out.col$Genus){ 
    #else they are just given the phyla name with the genus name
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    
  }
}

#sanity check
length(unique(data2plot$Phylum)) #- our number of phya stay the same
length(unique(data2plot$Genus)) #- our number of genus' has decreased because many low abundance genera have become "other"
sort(unique(data2plot$Genus))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "genus" with Other
    data2plot$Genus[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Genus))
unique(data2plot$Genus)
sort(unique(data2plot$Genus))

#add black to genera classified as "Other"
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))

#add colours to plotting dataframe
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
#set data as factors to keep our plotting in the right order
df$Subcommunity <- factor(df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))
df$Pole <- factor(df$Pole, levels = c("Antarctic", "Arctic"))
df$SubPole <- paste0(df$Pole, " (", df$Subcommunity, ")")
df$SubPole <- factor(df$SubPole, levels = c("Antarctic (Rare)", "Antarctic (Intermediate)", "Antarctic (Abundant)",
                                            "Arctic (Rare)", "Arctic (Intermediate)", "Arctic (Abundant)"))

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Genus as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Genus)), ]
levs = gen.df$Genus
levs = levs[levs != "Other"]
df.plot$Genus <- factor(df.plot$Genus, levels = out.col2$Genus)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Genus)

df.plot5 = 
  df.plot %>%
  group_by(Genus, SubPole, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p4 = ggplot(df.plot5, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p4

pdf("../results/eukaryote-cryoconite-per-sample-community-plot.pdf", width=35, height=10)
print(p4)
dev.off()

df.plot6 = 
  df.plot %>%
  group_by(Genus, SubPole, Glacier) %>%
  dplyr::summarise(across(c(Abundance), sum))


p5 = ggplot(df.plot6, aes(fill=Genus, y=Abundance, x=Glacier)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(nrow=8, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size = 20))

p5

pdf("../results/eukaryote-cryoconite-per-glacier-community-plot.pdf", width=35, height=10)
print(p5)
dev.off()

df.plot7 = 
  df.plot %>%
  group_by(Genus, SubPole) %>%
  dplyr::summarise(across(c(Abundance), sum))

p6 = ggplot(df.plot7, aes(fill=Genus, y=Abundance, x=SubPole)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(nrow=8, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20))

p6

pdf("../results/eukaryote-cryoconite-per-subcom-community-plot.pdf", width=20, height=10)
print(p6)
dev.off()

#COMBINE RELATIVE ABUNDANCE BARPLOTS AND PIECHARTS
# p1 + p3

###PUT THE TWO PLOTS TOGETHER
#put our pie charts together
# all.pie = pie1 | pie2 | pie3 | pie4

#get legend
legend <- get_legend(
  pie1  +
    guides(fill = guide_legend(nrow = 3, override.aes = list(size = 20))) +
    theme(legend.position = "bottom",
          legend.box="horizontal", legend.margin=margin(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text = element_text(size=50),
          legend.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
)



pie1 = pie1 + ggtitle("Antarctic Prokaryote") + theme(plot.title = element_text(size=45), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie2 = pie2 + ggtitle("Arctic Prokaryote") + theme(plot.title = element_text(size=45), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie3 = pie3 + ggtitle("Antarctic Eukaryote") + theme(plot.title = element_text(size=45), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie4 = pie4 + ggtitle("Arctic Eukaryote") + theme(plot.title = element_text(size=45), plot.margin = unit(c(0, 0, 0, 0), "cm"))

pie5 = pie5 + theme(plot.title = element_text(size=0), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie6 = pie6 + theme(plot.title = element_text(size=0), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie7 = pie7 + theme(plot.title = element_text(size=0), plot.margin = unit(c(0, 0, 0, 0), "cm"))
pie8 = pie8 + theme(plot.title = element_text(size=0), plot.margin = unit(c(0, 0, 0, 0), "cm"))

# all.pie1
# all.pie2 = plot_grid(all.pie1,
# legend,
# nrow=2)
# all.pie2

#combine with barcharts
# all.p = p1 / p3 / all.pie2 + plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 25))

#put pie charts and the legend together
# all.pie = pie1 | pie2 | pie3 | pie4 | legend
all.pie1 = plot_grid(pie1, pie2, pie3, pie4, legend, nrow=1)
# title1 <- ggdraw() + draw_label("Proportion of ASVs",fontface = 'bold', x = 0,hjust = 0, size = 50) +theme(
#     # add margin on the left of the drawing canvas,
#     # so title is aligned with left edge of first plot
#     plot.margin = margin(0, 0, 0, 7))
# 
# all.pie1 = plot_grid(
#   title, all.pie1,
#   ncol = 1,
#   # rel_heights values control vertical title margins
#   rel_heights = c(0.1, 1)
# )

all.pie2 = plot_grid(pie5, pie6, pie7, pie8, NULL, nrow=1)
all.pie2

# all.p = (all.pie) / p1 / p4 + plot_annotation(tag_levels = 'A') & 
#  theme(plot.tag = element_text(size = 55, face="bold"))

all.p = (all.pie1) / (all.pie2) / p1 / p4 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 55, face="bold"))

# all.p = plot_grid(all.pie,
#                   p1,
#                   p4,
#                   nrow=3,
#                   ncol=1,
#                   lables="AUTO")

pdf("../results/all-community-plot.pdf", width=35, height=30)
print(all.p)
dev.off()

####################################################################################################################

#output results to text document
sink("../results/community-profiles.txt", type="output")
writeLines("===============================================================
COMMUNITY PROFILES
===============================================================")
#Arctic abundant
writeLines("Percentage of reads per phylum in abundant Arctic prokaryote subcommunity:")
glom <- tax_glom(p.ar.a, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in intermediate Arctic prokaryote subcommunity:")
#Arctic intermediate
glom <- tax_glom(p.ar.i, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)


writeLines("Percentage of reads per phylum in rare Arctic prokaryote subcommunity:")
#Arctic rare
glom <- tax_glom(p.ar.r, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in abundant Antarctic prokaryote subcommunity:")
#Antarctic abundant
glom <- tax_glom(p.an.a, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in intermediate Antarctic prokaryote subcommunity:")
#Antarctic intermediate
glom <- tax_glom(p.an.i, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in rare Antarctic prokaryote subcommunity:")
#Antarctic rare
glom <- tax_glom(p.an.r, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in abundant Arctic eukaryote subcommunity:")
#what percentage of reads are chlorophyta in our subcommunities?
#Arctic abundant
glom <- tax_glom(e.ar.a, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in intermediate Arctic eukaryote subcommunity:")
#Arctic intermediate
glom <- tax_glom(e.ar.i, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in rare Arctic eukaryote subcommunity:")
#Arctic rare
glom <- tax_glom(e.ar.r, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in abundant Antarctic eukaryote subcommunity:")
#Antarctic abundant
glom <- tax_glom(e.an.a, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in intermediate Antarctic eukaryote subcommunity:")
#Antarctic intermediate
glom <- tax_glom(e.an.i, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

writeLines("Percentage of reads per phylum in rare Antarctic eukaryote subcommunity:")
#Antarctic rare
glom <- tax_glom(e.an.r, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort((rowSums(t.glom))/sum(t.glom), decreasing=TRUE)

#do we have any archaea in our datasets?
#remember microfauna have alread been removed

#Archaea in our Arctic dataset
#subset to archaea only  
xi = data.frame(tax_table(p.ar.a))[data.frame(tax_table(p.ar.a))$Domain == "Archaea",]
xii = data.frame(tax_table(p.ar.i))[data.frame(tax_table(p.ar.i))$Domain == "Archaea",]
xiii = data.frame(tax_table(p.ar.r))[data.frame(tax_table(p.ar.r))$Domain == "Archaea",]

#Archaea in our Antarctic dataset
yi = data.frame(tax_table(p.an.a))[data.frame(tax_table(p.an.a))$Domain == "Archaea",]
yii = data.frame(tax_table(p.an.i))[data.frame(tax_table(p.an.i))$Domain == "Archaea",]
yiii = data.frame(tax_table(p.an.r))[data.frame(tax_table(p.an.r))$Domain == "Archaea",]

writeLines("There are") 
nrow(xi) 
writeLines("archaeal ASVs in the abundant Arctic subcommunity.") 
xi
writeLines("There are") 
nrow(xii) 
writeLines("archaeal ASVs in the intermediate Arctic subcommunity.") 
xii
writeLines("There are") 
nrow(xiii) 
writeLines("archaeal ASVs in the rare Arctic subcommunity.") 
xiii
writeLines("There are") 
nrow(yi) 
writeLines("archaeal ASVs in the abundant Antarctic subcommunity.") 
yi
writeLines("There are") 
nrow(yii) 
writeLines("archaeal ASVs in the intermediate Antarctic subcommunity.") 
yii
writeLines("There are") 
nrow(yiii) 
writeLines("archaeal ASVs in the rare Antarctic subcommunity.") 
yiii


sink()





