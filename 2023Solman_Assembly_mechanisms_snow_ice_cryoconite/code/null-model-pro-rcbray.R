# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(scales) #for percentages
library(picante) #for faith's pd and comdistnt functions
library(ecodist) #for distance() function
library(parallel)
library(svglite)

ps.pro <- readRDS("../results/16S-ps-controls-removed.rds") 

#retain samples with >= num counts
ps.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro)
#remove ASVs with zero counts
ps.pro = filter_taxa(ps.pro, function(x) sum(x) >= 1, TRUE)

# C) Compute RC_Bray Model

# 1 Define function for calculating null values

# Function for calculating the distances in the null communities
RCbray_null_func <- function(i, freq.abd.df, alpha1, alpha2, N){
  # Get simulated communities and distance
  ## initally select OTUs weighted by their frequency. The number of OTUs selected should equal the richness of the samples.
  simcom1 = data.frame(table(sample(freq.abd.df$OTU, size=alpha1, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
  colnames(simcom1) = c("OTU","simcom1")
  simcom1$OTU = as.character(simcom1$OTU)
  simcom1 = inner_join(simcom1, freq.abd.df, by="OTU")
  simcom2 = data.frame(table(sample(freq.abd.df$OTU, size=alpha2, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
  colnames(simcom2) = c("OTU","simcom2")
  simcom2$OTU = as.character(simcom2$OTU)
  simcom2 = inner_join(simcom2, freq.abd.df, by="OTU")
  
  ## Now recruit OTUs based on their abundance in the metacommunity
  simcom1.abd = data.frame(table(sample(simcom1$OTU, size=N-alpha1, replace=T, prob=simcom1$p)), stringsAsFactors = F)
  colnames(simcom1.abd) = c("OTU","simcom1.abd")
  simcom1.abd$OTU = as.character(simcom1.abd$OTU)
  simcom1 = full_join(simcom1, simcom1.abd, by="OTU") %>%
    mutate(simcom1.abd = ifelse(is.na(simcom1.abd), 1, simcom1.abd)) %>%
    select(OTU, simcom1.abd)
  
  simcom2.abd = data.frame(table(sample(simcom2$OTU, size=N-alpha2, replace=T, prob=simcom2$p)), stringsAsFactors = F)
  colnames(simcom2.abd) = c("OTU","simcom2.abd")
  simcom2.abd$OTU = as.character(simcom2.abd$OTU)
  simcom2 = full_join(simcom2, simcom2.abd, by="OTU") %>%
    mutate(simcom2.abd = ifelse(is.na(simcom2.abd), 1, simcom2.abd)) %>%
    select(OTU, simcom2.abd)
  
  
  simcom = full_join(simcom1, simcom2, by="OTU")
  simcom[is.na(simcom)] = 0
  rownames(simcom) = simcom$OTU
  simcom$OTU = NULL
  
  #relative abundance communities
  simcom.rel = simcom
  simcom.rel$simcom1.abd = simcom$simcom1.abd/sum(simcom$simcom1.abd)
  simcom.rel$simcom2.abd = simcom$simcom2.abd/sum(simcom$simcom2.abd)
  
  null.dist.rel = vegdist(t(simcom.rel), method="bray")[1]
  return(null.dist.rel)
}


#2 Define main function for calculating RC Bray


# Main function for calculating RCbray
Calc_RCbray <- function(physeq, reps, nproc){
  
  # Get OTU table from phyloseq object
  otu.table = otu_table(physeq)
  
  #check ASV table
  df.count = data.frame(otu.table)
  
  # Get alpha diversity for each sample
  otu.PA.table = otu.table
  otu.PA.table[otu.PA.table > 0] = 1
  alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)
  
  #transform data for Bray-Curtis calculation
  physeq.rel = transform_sample_counts(physeq, function(x) x/sum(x))
  
  #get community data table
  rel.tab = data.frame(t(otu_table(physeq.rel)))
  
  # Get beta diversity matrix
  beta.table = as.matrix(vegdist(rel.tab), method="bray", diag=TRUE, upper=TRUE)
  
  ## Get metacommunity
  # Calculate the number of individuals in the meta community (Average read depth)
  N <- mean(apply(t(otu.table), 1, sum))
  
  # Calculate the average relative abundance of each taxa across communities
  p.m <- apply(t(otu.table), 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  
  # Calculate the occurrence frequency of each taxa across communities
  otu.table.bi <- 1*(t(otu.table)>0)
  freq <- apply(otu.table.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # Combine
  freq.abd.df = data.frame(p=p, freq=freq) %>%
    tibble::rownames_to_column(var="OTU") %>%
    filter(p != 0, freq != 0) %>%
    arrange(p)
  
  # For each pair of samples run the RCbray analysis
  comps = combn(alpha.df$Sample_ID, m=2, simplify = F)
  RCb.df = data.frame(Site1 = character(), Site2 = character(), RCb = numeric(), stringsAsFactors = F)
  for (j in seq(1, length(comps))){
    sam = comps[[j]]
    alpha1 = alpha.df[alpha.df$Sample_ID == sam[1],]$OTU.n
    alpha2 = alpha.df[alpha.df$Sample_ID == sam[2],]$OTU.n
    # Permute "reps" many times
    rep.list = seq(1, reps)
    null.list = mclapply(rep.list, RCbray_null_func, freq.abd.df=freq.abd.df, alpha1=alpha1, alpha2=alpha2, N=N, mc.cores=nproc)
    
    RCb = (length(null.list[null.list > beta.table[sam[1], sam[2]]]) + (0.5*length(null.list[null.list == beta.table[sam[1], sam[2]]])))/reps
    RCb = (RCb - 0.5)*2
    
    RCb.df = rbind(RCb.df, data.frame(Site1=sam[1], Site2=sam[2], RCb=RCb, stringsAsFactors = F))
  }
  
  RCb.df
  return(RCb.df)
}


get_rcbray_res <- function(ps, hab, group, habitat){
  
  phylo = prune_samples(hab, ps)
  #phylo
  df = Calc_RCbray(phylo, 999, 20) #should be 999 20
  df$Group = group
  df$Habitat = habitat
  return(df)
  
}

pro.sn.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
pro.sp.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice")
pro.sm.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice")
pro.cr.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite")

full.rcbray.df = rbind(pro.sn.rcbray.df, pro.sp.rcbray.df, pro.sm.rcbray.df, pro.cr.rcbray.df)

write.csv(full.rcbray.df, "../results/16S-RCbray-results-table.csv")
