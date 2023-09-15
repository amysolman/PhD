# rm(list=ls())
# my_phylo <- readRDS("../results/18S-phylo-object.rds")
# abun_cut = 0.001
# rare_cut = 0.0005

nine_cats <- function(my_phylo, rare_cut, abun_cut){
  
  my_phylo_trim <- prune_samples(sample_sums(my_phylo) >= 1, my_phylo)
  
  data = data.frame(otu_table(transform_sample_counts(my_phylo_trim, function(x) x / sum(x) )))
  
  #mean relative abundance for all ASVs
  mean_rel_abun <- rowMeans(data)
  
  #Create vectors to save our ASV names
  MAT <- vector()
  MRT <- vector()
  MMT <- vector()
  AAT <- vector()
  ART <- vector()
  AMT <- vector()
  CRT <- vector()
  CAT <- vector()
  CRAT <- vector()
  
  for (i in 1:nrow(data)){
    
    #if every value in the row is over 0.01 then save the row name to Always Abundant]
    x <- data[i,] #data for our ASV
    n <- length(x) #our number of samples
    a <- length(x[x >= abun_cut*10]) #number of samples the ASV is over or equal to 1% 
    b <- length(x[x <= rare_cut*10]) # less than or equal to 0.01%
    c <- length(x[x > rare_cut*10 & x < abun_cut*10]) # between 0.01% and 1%
    
    if (a == n){
      AAT <- c(AAT,row.names(data)[i])
    } else if (b == n){
      ART <- c(ART, row.names(data)[i])
    } else if (c == n){
      AMT <- c(AMT, row.names(data)[i])
    } else if (a == 0 & b > 0 & b < n){
      CRT <- c(CRT, row.names(data)[i])
    } else if (b == 0 & a > 0 & a < n){
      CAT <- c(CAT, row.names(data)[i])
    } else if (a > 0 & a < n & b > 0 & b < n){
      CRAT <- c(CRAT, row.names(data)[i])
    }
  }
  
  for (i in 1:length(mean_rel_abun)){
    
    if (mean_rel_abun[[i]] >= abun_cut){
      MAT <- c(MAT,names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] <= rare_cut){
      MRT <- c(MRT, names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] < abun_cut & mean_rel_abun[[i]] > rare_cut){
      MMT <- c(MMT, names(mean_rel_abun)[i])
    } 
  }
  
  #bind and return our results
  MAT_res <- cbind(rep("MAT", length(MAT)), MAT)
  MRT_res <- cbind(rep("MRT", length(MRT)), MRT)
  MMT_res <- cbind(rep("MMT", length(MMT)), MMT)
  AAT_res <- cbind(rep("AAT", length(AAT)), AAT)
  ART_res <- cbind(rep("ART", length(ART)), ART)
  AMT_res <- cbind(rep("AMT", length(AMT)), AMT)
  CRT_res <- cbind(rep("CRT", length(CRT)), CRT)
  CAT_res <- cbind(rep("CAT", length(CAT)), CAT)
  CRAT_res <- cbind(rep("CRAT", length(CRAT)), CRAT)
  
  nine_cat_res <- rbind(AAT_res, ART_res, AMT_res, CRT_res, CAT_res, CRAT_res, MAT_res, MRT_res, MMT_res)
  
  return(nine_cat_res)
  
}

# my_phylo = ps_arc
# abun_cut = euk_abun_cut
# rare_cut = euk_rare_cut

three_cats <- function(my_phylo, rare_cut, abun_cut){
  
  nsamples(my_phylo)
  
  my_phylo_trim <- prune_samples(sample_sums(my_phylo) >= 1, my_phylo)
  
  nsamples(my_phylo_trim)
  
  data = data.frame(otu_table(transform_sample_counts(my_phylo_trim, function(x) x / sum(x) )))
  
  nrow(data) == ntaxa(my_phylo)
  colSums(data)
  
  #mean relative abundance for all ASVs
  mean_rel_abun <- rowMeans(data)
  
  #Create vectors to save our ASV names
  MAT <- vector()
  MRT <- vector()
  MMT <- vector()
  
  for (i in 1:length(mean_rel_abun)){
    
    if (mean_rel_abun[[i]] >= abun_cut){
      MAT <- c(MAT,names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] <= rare_cut){
      MRT <- c(MRT, names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] < abun_cut & mean_rel_abun[[i]] > rare_cut){
      MMT <- c(MMT, names(mean_rel_abun)[i])
    } 
  }
  
  length(MAT) + length(MMT) + length(MRT) == ntaxa(my_phylo) #we've classified all ASVs
  
  #bind and return our results
  MAT_res <- cbind(rep("MAT", length(MAT)), MAT)
  MRT_res <- cbind(rep("MRT", length(MRT)), MRT)
  MMT_res <- cbind(rep("MMT", length(MMT)), MMT)
  
  three_cat_res <- rbind(MAT_res, MRT_res, MMT_res)
  
  nrow(three_cat_res) == ntaxa(my_phylo)
  
  return(three_cat_res)
  
}