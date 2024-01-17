# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(parallel)
library(svglite)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds")
#eukaryotes without microfauna
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds")
#microfauna only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds")

#Mantel function
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

  #MANTEL CORRELOGRAM ONE
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

euk.sn <- mantel_correlogram(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID)
euk.sp <- mantel_correlogram(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID)
euk.sm <- mantel_correlogram(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID)
euk.cr <- mantel_correlogram(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID)

mm.sn <- mantel_correlogram(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID)
mm.sp <- mantel_correlogram(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID)
mm.sm <- mantel_correlogram(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID)
mm.cr <- mantel_correlogram(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID)

#save results
pro.sn$Group = "Prokaryote"
pro.sn$Habitat = "Snow"
pro.sp$Group = "Prokaryote"
pro.sp$Habitat = "Spring Ice"
pro.sm$Group = "Prokaryote"
pro.sm$Habitat = "Summer Ice"
pro.cr$Group = "Prokaryote"
pro.cr$Habitat = "Cryoconite"
euk.sn$Group = "Eukaryote"
euk.sn$Habitat = "Snow"
euk.sp$Group = "Eukaryote"
euk.sp$Habitat = "Spring Ice"
euk.sm$Group = "Eukaryote"
euk.sm$Habitat = "Summer Ice"
euk.cr$Group = "Eukaryote"
euk.cr$Habitat = "Cryoconite"
mm.sn$Group = "Microfauna"
mm.sn$Habitat = "Snow"
mm.sp$Group = "Microfauna"
mm.sp$Habitat = "Spring Ice"
mm.sm$Group = "Microfauna"
mm.sm$Habitat = "Summer Ice"
mm.cr$Group = "Microfauna"
mm.cr$Habitat = "Cryoconite"

man.pro = rbind(pro.sn, pro.sp, pro.sm, pro.cr)
man.euk = rbind(euk.sn, euk.sp, euk.sm, euk.cr)
man.mm = rbind(mm.sn, mm.sp, mm.sm, mm.cr)

#export results
write.csv(man.pro, "../results/prokaryote-mantel-results.csv")
write.csv(man.euk, "../results/eukaryote-mantel-results.csv")
write.csv(man.mm, "../results/microfauna-mantel-results.csv")