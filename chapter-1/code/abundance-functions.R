six_cats <- function(my_phylo, rare_cut, abun_cut){
  
  data = data.frame(otu_table(transform_sample_counts(my_phylo, function(x) x / sum(x) )))
  
  #Create vectors to save our ASV names
  AT <- vector()
  RT <- vector()
  MT <- vector()
  CRT <- vector()
  CAT <- vector()
  CRAT <- vector()
  
  for (i in 1:nrow(data)){
    
    #if every value in the row is over 0.01 then save the row name to Always Abundant]
    x <- data[i,] #data for our ASV
    n <- length(x) #our number of samples
    a <- length(x[x >= abun_cut]) #number of samples the ASV is over or equal to 1% 
    b <- length(x[x <= rare_cut]) # less than or equal to 0.01%
    c <- length(x[x > rare_cut & x < abun_cut]) # between 0.01% and 1%
    
    if (a == n){
      AT <- c(AT,row.names(data)[i])
    } else if (b == n){
      RT <- c(RT, row.names(data)[i])
    } else if (c == n){
      MT <- c(MT, row.names(data)[i])
    } else if (a == 0 & b > 0 & b < n){
      CRT <- c(CRT, row.names(data)[i])
    } else if (b == 0 & a > 0 & a < n){
      CAT <- c(CAT, row.names(data)[i])
    } else if (a > 0 & a < n & b > 0 & b < n){
      CRAT <- c(CRAT, row.names(data)[i])
    }
  }
  
  #bind and return our results
  AT_res <- cbind(rep("AT", length(AT)), AT)
  RT_res <- cbind(rep("RT", length(RT)), RT)
  MT_res <- cbind(rep("MT", length(MT)), MT)
  CRT_res <- cbind(rep("CRT", length(CRT)), CRT)
  CAT_res <- cbind(rep("CAT", length(CAT)), CAT)
  CRAT_res <- cbind(rep("CRAT", length(CRAT)), CRAT)
  
  six_cat_res <- rbind(AT_res, RT_res, MT_res, CRT_res, CAT_res, CRAT_res)
  
  return(six_cat_res)
  
}

sample_three_cats <- function(my_phylo, rare_cut, abun_cut){
  
  data = data.frame(otu_table(transform_sample_counts(my_phylo, function(x) x / sum(x) )))
  
  #mean relative abundance for all ASVs
  mean_rel_abun <- rowMeans(data)
  
  #Create vectors to save our ASV names
  AT <- vector()
  RT <- vector()
  MT <- vector()
  
  for (i in 1:length(mean_rel_abun)){
    
    if (mean_rel_abun[[i]] >= abun_cut){
      AT <- c(AT,names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] <= rare_cut){
      RT <- c(RT, names(mean_rel_abun)[i])
    } else if (mean_rel_abun[[i]] < abun_cut & mean_rel_abun[[i]] > rare_cut){
      MT <- c(MT, names(mean_rel_abun)[i])
    } 
  }
  
  #bind and return our results
  AT_res <- cbind(rep("AT", length(AT)), AT)
  RT_res <- cbind(rep("RT", length(RT)), RT)
  MT_res <- cbind(rep("MT", length(MT)), MT)
  
  sample_three_cat_res <- rbind(AT_res, RT_res, MT_res)
  
  return(sample_three_cat_res)
  
}

data_three_cats <- function(my_phylo, rare_cut, abun_cut){
  
  #relative abundance for each ASV across the whole dataset
  data <- data.frame(otu_table(my_phylo))
  ASV_rel_abun <- rowSums(data)/sum(colSums(data))
  
  #Create vectors to save our ASV names
  AT <- vector()
  RT <- vector()
  MT <- vector()
  
  for (i in 1:length(ASV_rel_abun)){
    
    if (ASV_rel_abun[[i]] >= abun_cut){
      AT <- c(AT,names(ASV_rel_abun)[i])
    } else if (ASV_rel_abun[[i]] <= rare_cut){
      RT <- c(RT, names(ASV_rel_abun)[i])
    } else if (ASV_rel_abun[[i]] < abun_cut & ASV_rel_abun[[i]] > rare_cut){
      MT <- c(MT, names(ASV_rel_abun)[i])
    } 
  }
  
  #bind and return our results
  AT_res <- cbind(rep("AT", length(AT)), AT)
  RT_res <- cbind(rep("RT", length(RT)), RT)
  MT_res <- cbind(rep("MT", length(MT)), MT)
  
  data_three_cat_res <- rbind(AT_res, RT_res, MT_res)
  
  return(data_three_cat_res)
  
}