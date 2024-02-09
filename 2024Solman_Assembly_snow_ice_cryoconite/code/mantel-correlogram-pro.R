# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(ecodist) #for distance() function
library(parallel)
library(svglite)

ps.pro <- readRDS("../results/16S-ps-controls-removed.rds") 

#retain samples with >= num counts
ps.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro)
#remove ASVs with zero counts
ps.pro = filter_taxa(ps.pro, function(x) sum(x) >= 1, TRUE)


mantel_correlogram <- function(ps, hab){
  
  set.seed(666)
  #Separate habitat data
  h <- prune_samples(hab, ps)
  
  #PROPORTIONAL TRANSFORMATION
  ps.rel <- transform_sample_counts(h, function(x) x/sum(x))
  
  #Get metadata of samples
  meta <- data.frame(sample_data(ps.rel))
  
  #list of variables we're interested in 
  meta = meta[,-c(1:8, 11:13)]
  #remove columns with all NA values
  meta.trim <- meta[ , colSums(is.na(meta)) < nrow(meta)] 
  #remove columns with less than 75% values > 0
  c <- meta.trim
  c[c > 0 | c < 0] <- 1
  c[is.na(c)] <- 0
  meta.trim <- meta.trim[,colSums(c) >= round(nrow(meta.trim)*0.75)]
  
  #Only keep complete cases
  meta.final <- meta.trim[complete.cases(meta.trim),]
  
  #remove samples from our phyloseq object that aren't in our meta data
  ps.smol = prune_samples(rownames(meta.final), ps.rel)
  
  #remove ASVs with 0 relative abundance
  trim.filter <- filter_taxa(ps.smol, function(x) sum(x) > 0, TRUE)
  
  #Get phylogenetic distances
  #pull out phylogenetic tree
  tree <- phy_tree(trim.filter) 
  #compute phylogenetic tree distances between each ASV
  tree.dist <- cophenetic(tree)
  # #give the tree labels (ASV IDs) to the vector asv
  asvs <- tree$tip.label
  #make sure they have the same ASVs as rows and coloums
  phylo.dists <- tree.dist[asvs, asvs]
  #make the upper triangle of the phylogeneic distance object NA
  phylo.dists[upper.tri(phylo.dists, diag=TRUE)] = NA
  
  #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  #turn the metadata into a matrix
  niches <- as.matrix(meta.final) 
  #get our relative abundance ASV table
  asv.table <- data.frame(t(otu_table(trim.filter)), check.names = FALSE)
  #this will find the weighted mean environmental parameters for each asv
  asv.niche <- data.frame(wascores(niches, asv.table), check.names = FALSE) 
  #generate euclidean distance matrix for each ASV using combined environmental parameters
  dist.out = as.matrix(dist(asv.niche), labels=TRUE) 
  
  #calculate mantel correlogram
  corlg <- mantel.correlog(dist.out, phylo.dists, r.type="spearman")
  
  # Prep data
  crlg <- data.frame(corlg$mantel.res)
  
  crlg <- crlg %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  
  
  crlg$D.cl = rownames(crlg)
  
  return(crlg)
}


pro.sn <- mantel_correlogram(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID)
pro.sp <- mantel_correlogram(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID)
pro.sm <- mantel_correlogram(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID)
pro.cr <- mantel_correlogram(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID)

#combine into single dataframe
pro.sn$Group = "Prokaryote"
pro.sn$Habitat = "Snow"
pro.sp$Group = "Prokaryote"
pro.sp$Habitat = "Spring Ice"
pro.sm$Group = "Prokaryote"
pro.sm$Habitat = "Summer Ice"
pro.cr$Group = "Prokaryote"
pro.cr$Habitat = "Cryoconite"

plot.df = rbind(pro.sn, pro.sp, pro.sm, pro.cr)

write.csv(plot.df, "../results/prokaryote-mantel-results.csv")