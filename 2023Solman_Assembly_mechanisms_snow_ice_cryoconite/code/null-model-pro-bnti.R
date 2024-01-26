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

#import data
ps.pro <- readRDS("../results/16S-ps-controls-removed.rds") 

#retain samples with >= num counts
ps.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro)
#remove ASVs with zero counts
ps.pro = filter_taxa(ps.pro, function(x) sum(x) >= 1, TRUE)

# #rarefy data to an even depth - possibly remove this?
# set.seed(72)  # setting seed for reproducibility
# ps.pro = rarefy_even_depth(ps.pro)

#For testing
# #reduced size dataset for checked code
# ps.pro = rarefy_even_depth(ps.pro, sample.size = 50)
# #proportionally transform to look at the difference
# ps.pro = transform_sample_counts(ps.pro, function(x) x/sum(x)) 


#Code taken from Barnett et al., (2020) was used for the model fitting (https://github.com/seb369/landuse_comm_assembly).  

#Compute Beta-nearest Taxon Index (βNTI)

#Define function for calculating the βMNTD for each random null community

# Function for calculating the βMNTD for each random null community
bMNTD_null_func <- function(i, OTU.table, tree){
  
  #randomly re-shuffle the tip lables of the phylogenetic tree
  tree$tip.label = sample(tree$tip.label)
  
  #get the phylogenetic distances between the reshuffled communities
  bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  
  #collects the same data as the previous function
  A <- attr(bMNTD_s, "Size")
  B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
  if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
  bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_s),
                          rep=i)
  return(bMNTD_s.df)
}


#Define main function for calculating BNTI

# #testing the function
# physeq = prune_samples(data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, ps.pro)
# #remove taxa with zero counts
# physeq = filter_taxa(physeq, function(x) sum(x) > 0, TRUE)
# #remove samples with zero counts
# physeq = prune_samples(sample_sums(physeq) > 0, physeq)
# reps = 10
# nproc = 2

# The main function for calculating βNTI
Phylo_turnover <- function(physeq, reps, nproc){
  
  # Extract OTU table
  OTU.table = t(otu_table(physeq))
  
  # Extract phylogenetic tree
  tree = phy_tree(physeq)
  
  # Get βMNTD between all communities
  #so for each pair of communities, for each taxa in those communities
  #this function finds the most closely related taxa and finds the phylogenetic distance to that taxa
  #this is weighted by the abundance of taxa
  #and repeated for each ASV
  #the mean of all ASVs is then calculated 
  #to give a metric of phylogenetic relatedness between two communities
  
  #WHY WE DON'T USE PROPORTIONAL DATA
  #When abundance.weighted = TRUE the metric changes from mean phylogenetic distance among taxa
  #to the mean phylogenetic distance among individuals
  #if we use relative abundance data, the difference in abundance of an individual from sample 1
  #to sample 2 is lost and reduced to a decimal number
  #this essentially obscures the abudance.weighted aspect of the function
  bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  
  #lets have a look
  print(bMNTD_o)
  
  #A = number of samples
  A <- attr(bMNTD_o, "Size")
  
  #B = sample names
  B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
  
  #remove diagonal values of distance matrix
  if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
  
  #lets have a look
  print(bMNTD_o)
  
  #remove upper triangle of distance matrix
  if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
  
  #lets have a look
  print(bMNTD_o)
  
  #put the distance matrix values into a long format dataframe
  bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_o))
  
  # Get βMNTD for randomized null communities
  
  #make a vector of number for each randomized null community we will make
  rep.list = seq(1, reps)
  
  #for each rep indicated by rep.list, run the function bMNTD_null_func
  #with the variables OTU.table from our environmental samples
  #the phylogenetic tree from our environmental samples
  #and the assigned number of cores (processing units) to use to run this analysis 
  bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)
  
  # Combine all data together and calculate βNTI for each sample pair
  bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
  bMNTD_s.means.df = bMNTD_s.df %>%
    group_by(Sample_1, Sample_2) %>%
    dplyr::summarize(mean_bMNTD = mean(bMNTD),
                     sd_bMNTD = sd(bMNTD))
  
  bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
    mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)
  return(bMNTD_o.df)
}

get_bnti_res <- function(ps, hab, group, habitat){

  #subset to habitat
  phylo = prune_samples(hab, ps)
  #rarefy data to an even depth
  set.seed(72)  # setting seed for reproducibility
  phylo = rarefy_even_depth(phylo)
  #remove taxa with zero counts
  phylo = filter_taxa(phylo, function(x) sum(x) > 0, TRUE)
  #remove samples with zero counts
  phylo = prune_samples(sample_sums(phylo) > 0, phylo)
  
  #run function on rarefied data
  df = Phylo_turnover(phylo, 1000, 10) #should be 1000 10
  
  df$Group = group
  df$Habitat = habitat
  return(df)
}

#test
# bnti.test.rel = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
# bnti.test.rare = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")

pro.sn.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
pro.sp.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice")
pro.sm.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice")
pro.cr.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite")



full.bNTI.df = rbind(pro.sn.bNTI.df, pro.sp.bNTI.df, pro.sm.bNTI.df, pro.cr.bNTI.df)

write.csv(full.bNTI.df, "../results/16S-bNTI-results-table.csv")