# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Define the abundance cut offs
# 3. Import data
# 4. Use the previously written abundance functions to assign each ASV to an abundance category and create new phyloseq objects 
# 5. Export the phyloseq objects

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
library(scales) #for percentages
# install.packages("indicspecies")
library(indicspecies)
library(ggrepel)

# 2. Define the abundance cut offs
# Prokaryote – abundant >0.1% RA, rare <0.01% RA
# Eukaryote – abundant >0.1% RA, rare <0.05% RA

pro_abun_cut = 0.001
pro_rare_cut = 0.0001 
euk_abun_cut = 0.001
euk_rare_cut = 0.0005

# Method

#Using the appropriate abundance cut offs defined by MultiCoLA, the average rarefied datasets were divided into abundant, 
#intermediate and rare subcommunities. Abundant porkaryote ASVs were those with relative abundance above `r pro_abun_cut*100`%. 
#Rare prokaryote ASVs were those with relative abundance below `r pro_rare_cut*100`%. 
#Abundant eukaryote ASVs were those with relative abundance above `r euk_abun_cut*100`%. 
#Rare eukaryote ASVs were those with relative abundance below `r euk_rare_cut*100`%. 

# Results

# 3. Import data and separate into Arctic and Antarctic samples

pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds") 

# 4. Use the previously written abundance functions to assign each ASV to an abundance category and create new phyloseq objects 

get_coms <- function(phylo, abun_cut, rare_cut){
  
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
  ntaxa(ps_arc) + ntaxa(ps_ant) - length(shared) == ntaxa(phylo)

  #Arctic
  #classify each ASV as abundant, intermediate or rare
  res_arc <- as.data.frame(three_cats(ps_arc, rare_cut, abun_cut))
  nrow(res_arc) == ntaxa(ps_arc) #this should equal true
  
  #Antarctic
  #classify each ASV as abundant, intermediate or rare
  res_ant <- as.data.frame(three_cats(ps_ant, rare_cut, abun_cut))
  nrow(res_ant) == ntaxa(ps_ant) #this should equal true

  x1 = data.frame(table(res_arc$V1))
  x1$Data = "Arctic"
  x2 = data.frame(table(res_ant$V1))
  x2$Data = "Antarctic"

  #number of ASVs in each subcommunity
  df = rbind(x1, x2)

  #format table for putting into Markdown document
  df2 = reshape(df, idvar = "Data", timevar = "Var1", direction = "wide")
  df2[is.na(df2)] <- 0
  df3 = df2[, c(1, ncol(df2)-2, ncol(df2)-1, ncol(df2))]
  names(df3) = c("Data", "Abundant", "Intermediate", "Rare")

  #add percentages to df
  df4 = data.frame(df3$Data, df3$Abundant, "Abundant (%)" = c(round(df3$Abundant[1]/sum(df3$Abundant[1], df3$Intermediate[1], df3$Rare[1]), 3)*100, round(df3$Abundant[2]/sum(df3$Abundant[2], df3$Intermediate[2], df3$Rare[2]), 3)*100),
                   df3$Intermediate, "Intermediate (%)" = c(round(df3$Intermediate[1]/sum(df3$Abundant[1], df3$Intermediate[1], df3$Rare[1]), 3)*100, round(df3$Intermediate[2]/sum(df3$Abundant[2], df3$Intermediate[2], df3$Rare[2]), 3)*100),
                   df3$Rare, "Rare (%)" = c(round(df3$Rare[1]/sum(df3$Abundant[1], df3$Intermediate[1], df3$Rare[1]), 3)*100, round(df3$Rare[2]/sum(df3$Abundant[2], df3$Intermediate[2], df3$Rare[2]), 3)*100)
                 )
  names(df4) = c("Data", "Abundant", "Abundant (%)", "Intermediate", "Intermediate (%)", "Rare", "Rare (%)")


#Divide the phyloseq object by into abundance subcommuntiies by pole
  
my_subset_arc_abun <- subset(otu_table(ps_arc), rownames(otu_table(ps_arc)) %in% res_arc[res_arc$V1 == "MAT",]$MAT)
new_physeq_arc_abun <- merge_phyloseq(my_subset_arc_abun, tax_table(ps_arc), sample_data(ps_arc), phy_tree(ps_arc))
ntaxa(new_physeq_arc_abun) == df3$Abundant[1] #this should be true

my_subset_arc_int <- subset(otu_table(ps_arc), rownames(otu_table(ps_arc)) %in% res_arc[res_arc$V1 == "MMT",]$MAT)
new_physeq_arc_int <- merge_phyloseq(my_subset_arc_int, tax_table(ps_arc), sample_data(ps_arc), phy_tree(ps_arc))
ntaxa(new_physeq_arc_int) == df3$Intermediate[1] #this should be true

my_subset_arc_rare <- subset(otu_table(ps_arc), rownames(otu_table(ps_arc)) %in% res_arc[res_arc$V1 == "MRT",]$MAT)
new_physeq_arc_rare <- merge_phyloseq(my_subset_arc_rare, tax_table(ps_arc), sample_data(ps_arc), phy_tree(ps_arc))
ntaxa(new_physeq_arc_rare) == df3$Rare[1] #this should be true
#ntaxa(new_physeq_arc_rare)

#do all samples have a rare community?
test = data.frame(otu_table(new_physeq_arc_rare))
min(colSums(test)) == 0 #if this equals true then some samples don'e have a "rare" community

#Divide the phyloseq object by into abundance subcommuntiies by pole
my_subset_ant_abun <- subset(otu_table(ps_ant), rownames(otu_table(ps_ant)) %in% res_ant[res_ant$V1 == "MAT",]$MAT)
new_physeq_ant_abun <- merge_phyloseq(my_subset_ant_abun, tax_table(ps_ant), sample_data(ps_ant), phy_tree(ps_ant))
ntaxa(new_physeq_ant_abun) == df3$Abundant[2] #this should be true

my_subset_ant_int <- subset(otu_table(ps_ant), rownames(otu_table(ps_ant)) %in% res_ant[res_ant$V1 == "MMT",]$MAT)
new_physeq_ant_int <- merge_phyloseq(my_subset_ant_int, tax_table(ps_ant), sample_data(ps_ant), phy_tree(ps_ant))
ntaxa(new_physeq_ant_int) == df3$Intermediate[2] #this should be true

my_subset_ant_rare <- subset(otu_table(ps_ant), rownames(otu_table(ps_ant)) %in% res_ant[res_ant$V1 == "MRT",]$MAT)
new_physeq_ant_rare <- merge_phyloseq(my_subset_ant_rare, tax_table(ps_ant), sample_data(ps_ant), phy_tree(ps_ant))
ntaxa(new_physeq_ant_rare) == df3$Rare[2] #this should be true

#table with number and percentage of sequences
arc_abun_seq = sum(data.frame(otu_table(new_physeq_arc_abun)))
arc_int_seq = sum(data.frame(otu_table(new_physeq_arc_int)))
arc_rare_seq = sum(data.frame(otu_table(new_physeq_arc_rare)))
ant_abun_seq = sum(data.frame(otu_table(new_physeq_ant_abun)))
ant_int_seq = sum(data.frame(otu_table(new_physeq_ant_int)))
ant_rare_seq = sum(data.frame(otu_table(new_physeq_ant_rare)))

#table
df5 = data.frame(df3$Data, "Abundant" = c(arc_abun_seq, ant_abun_seq), "Abundant (%)" = c(round(arc_abun_seq/sum(arc_abun_seq, arc_int_seq, arc_rare_seq), 3)*100, round(ant_abun_seq/sum(ant_abun_seq, ant_int_seq, ant_rare_seq), 3)*100),
                 
                 "Intermediate" = c(arc_int_seq, ant_int_seq), "Intermediate (%)" = c(round(arc_int_seq/sum(arc_abun_seq, arc_int_seq, arc_rare_seq), 3)*100, round(ant_int_seq/sum(ant_abun_seq, ant_int_seq, ant_rare_seq), 3)*100),
                 
                 "Rare" = c(arc_rare_seq, ant_rare_seq), "rare (%)" = c(round(arc_rare_seq/sum(arc_abun_seq, arc_int_seq, arc_rare_seq), 3)*100, round(ant_rare_seq/sum(ant_abun_seq, ant_int_seq, ant_rare_seq), 3)*100) )

names(df5) = c("Data", "Abundant", "Abundant (%)", "Intermediate", "Intermediate (%)", "Rare", "Rare (%)")



res.list = list(df4, df5, res_arc, res_ant, new_physeq_arc_abun, new_physeq_arc_int, new_physeq_arc_rare, new_physeq_ant_abun, new_physeq_ant_int, new_physeq_ant_rare)

return(res.list)
}

pro_res = get_coms(pro, pro_abun_cut, pro_rare_cut)
euk_res = get_coms(euk, euk_abun_cut, euk_rare_cut)

# 5. Export the phyloseq objects

saveRDS(pro_res[[5]], "../results/16S-phylo-object-arc-abun.rds")
saveRDS(pro_res[[6]], "../results/16S-phylo-object-arc-int.rds")
saveRDS(pro_res[[7]], "../results/16S-phylo-object-arc-rare.rds")
saveRDS(pro_res[[8]], "../results/16S-phylo-object-ant-abun.rds")
saveRDS(pro_res[[9]], "../results/16S-phylo-object-ant-int.rds")
saveRDS(pro_res[[10]], "../results/16S-phylo-object-ant-rare.rds")
saveRDS(euk_res[[5]], "../results/18S-phylo-object-arc-abun.rds")
saveRDS(euk_res[[6]], "../results/18S-phylo-object-arc-int.rds")
saveRDS(euk_res[[7]], "../results/18S-phylo-object-arc-rare.rds")
saveRDS(euk_res[[8]], "../results/18S-phylo-object-ant-abun.rds")
saveRDS(euk_res[[9]], "../results/18S-phylo-object-ant-int.rds")
saveRDS(euk_res[[10]], "../results/18S-phylo-object-ant-rare.rds")