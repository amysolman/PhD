# Method 

# The relative contribution of stochastic and deterministic processes in shaping cryoconite communities was quantified using a phylogenetic null model (Stegen et al., 2013). β mean nearest taxon distance metric (βMNTD) was calculated using the comdistnt function of the picante package in R (Kembel et al., 2010). For each pair of communities, the null βMNTD was calculated by randomly shuffling the tips of the phylogeny and recalculating the βMNTD 999 times to get distribution of null values.
# 
# Observed βMNTD and the null βMNTD were used to calculate β nearest taxon index (βNTI). Lower values of βNTI ( < -2) indicates ASVs are phylogenetically clustered due to homogenous selection. Higher values of βNTI ( > +2) indicates ASVs are variable selection. Where |βNTI| < +2 indicates stochastic processes are predominantly driving community assembly. 
# 
# As stochastic processes include homogenizing dispersal, disprsal limitation and ecological drift, the Bray–Curtis-based Raup–Crick metric (RCbray) is calculated. Observed Bray-Curtis values are compared to a null distribution to give the RCbray metric. Dispersal limitation is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and RCbray > +0.95. Homogenizing dispersal is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and RCbray < -0.95. Ecological drift is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and |RCbray| < +0.95.
# 
# Phylogenetic null models rely on significant phylogenetic signal within the data. Taxa that are more closely related must exhibit niche similarities. To test this, Mantel correlation plots between Cophenetic distances and Eucldean distances between niche preferences were generated with 999 random permutations for test of significance. Cophenetic distances between each ASV were calculated. The abundance-weighted mean of each environmental variable for each ASV was calculated using function wascores in ‘vegan’ package (Oksanen et al., 2022). Euclidean distances between niche preferences for each ASV were calculated. Mantel correlation plots could not be produced for Arctic dataset as environmental data was not available. 
# 
# Variation in community assembly processes along each environmental gradient was assessed using regression analysis of BNTI values and Euclidean distances of major environmental variables. Mantel tests with 9999 permutations were used to test statistical significance. This method was also used to assess the relationship between phylogenetic turnover and environmental factors after controlling for geographic distance. This was carried out using the mantel function in the ecodist package. Code adapted from 
# #https://github.com/seb369/landuse_comm_assembly.  

rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan) #wascores
library(ggplot2)
library(picante)
library(ecodist) #for distance() function
library(parallel)
library(svglite)
library(dplyr) #for %>%
library(scales) #for percentages

# Step 2: Read in data
pro <- readRDS("../results/16S-phylo-object-rarefied-var-trans.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied-var-trans.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) > 0, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) > 0, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) > 0, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) > 0, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.int <- prune_samples(sample_sums(euk.ant.int)>0, euk.ant.int)
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")

pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds") 
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds") 
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds") 
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds")

#Make each dataset mini so we can locally test the whole script
# pro = rarefy_even_depth(pro, sample.size = min(sample_sums(pro))/100)
# euk = rarefy_even_depth(euk, sample.size = min(sample_sums(euk))/100)
# 
# pro.arc = rarefy_even_depth(pro.arc, sample.size = min(sample_sums(pro.arc))/100)
# pro.ant = rarefy_even_depth(pro.ant, sample.size = min(sample_sums(pro.ant))/1000)
# pro.ant <- filter_taxa(pro.ant, function(x) sum(x) > 0, TRUE)
# euk.arc = rarefy_even_depth(euk.arc, sample.size = min(sample_sums(euk.arc))/100)
# euk.ant = rarefy_even_depth(euk.ant, sample.size = min(sample_sums(euk.ant))/200)
# 
# pro.arc.abun = rarefy_even_depth(pro.arc.abun, sample.size = min(sample_sums(pro.arc.abun))/100)
# pro.arc.int = rarefy_even_depth(pro.arc.int, sample.size = min(sample_sums(pro.arc.int))/10)
# pro.arc.rare = rarefy_even_depth(pro.arc.rare, sample.size = min(sample_sums(pro.arc.rare))/10)
# pro.ant.abun = rarefy_even_depth(pro.ant.abun, sample.size = min(sample_sums(pro.ant.abun))/100)
# pro.ant.int = rarefy_even_depth(pro.ant.int, sample.size = min(sample_sums(pro.ant.int))/100)
# pro.ant.rare = rarefy_even_depth(pro.ant.rare, sample.size = min(sample_sums(pro.ant.rare))/5)
# 
# euk.arc.abun = rarefy_even_depth(euk.arc.abun, sample.size = min(sample_sums(euk.arc.abun))/100)
# euk.arc.int = rarefy_even_depth(euk.arc.int, sample.size = min(sample_sums(euk.arc.int))/10)
# euk.arc.rare = rarefy_even_depth(euk.arc.rare, sample.size = min(sample_sums(euk.arc.rare))/10)
# euk.ant.abun = rarefy_even_depth(euk.ant.abun, sample.size = min(sample_sums(euk.ant.abun))/100)
# euk.ant.int = rarefy_even_depth(euk.ant.int, sample.size = min(sample_sums(euk.ant.int))*10)
# euk.ant.rare = rarefy_even_depth(euk.ant.rare, sample.size = min(sample_sums(euk.ant.rare))*10)

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



get_rcbray_res <- function(phylo, subcom, group, pole){
  
  
  df = Calc_RCbray(phylo, 999, 20) #should be 999 20
  df$Group = group
  df$Subcommunity = subcom
  df$Pole = pole
  return(df)
  
}

pro.ant.rcbray.df = get_rcbray_res(pro.ant, "Full", "Prokaryote", "Antarctic")
pro.ant.abun.rcbray.df = get_rcbray_res(pro.ant.abun, "Abundant", "Prokaryote", "Antarctic")
pro.ant.int.rcbray.df = get_rcbray_res(pro.ant.int, "Intermediate", "Prokaryote", "Antarctic")
pro.ant.rare.rcbray.df = get_rcbray_res(pro.ant.rare, "Rare", "Prokaryote", "Antarctic")
pro.arc.rcbray.df = get_rcbray_res(pro.arc, "Full", "Prokaryote", "Arctic")
pro.arc.abun.rcbray.df = get_rcbray_res(pro.arc.abun, "Abundant", "Prokaryote", "Arctic")
pro.arc.int.rcbray.df = get_rcbray_res(pro.arc.int, "Intermediate", "Prokaryote", "Arctic")
pro.arc.rare.rcbray.df = get_rcbray_res(pro.arc.rare, "Rare", "Prokaryote", "Arctic")
euk.ant.rcbray.df = get_rcbray_res(euk.ant, "Full", "Eukaryote", "Antarctic")
euk.ant.abun.rcbray.df = get_rcbray_res(euk.ant.abun, "Abundant", "Eukaryote", "Antarctic")
euk.ant.int.rcbray.df = get_rcbray_res(euk.ant.int, "Intermediate", "Eukaryote", "Antarctic")
euk.ant.rare.rcbray.df = get_rcbray_res(euk.ant.rare, "Rare", "Eukaryote", "Antarctic")
euk.arc.rcbray.df = get_rcbray_res(euk.arc, "Full", "Eukaryote", "Arctic")
euk.arc.abun.rcbray.df = get_rcbray_res(euk.arc.abun, "Abundant", "Eukaryote", "Arctic")
euk.arc.int.rcbray.df = get_rcbray_res(euk.arc.int, "Intermediate", "Eukaryote", "Arctic")
euk.arc.rare.rcbray.df = get_rcbray_res(euk.arc.rare, "Rare", "Eukaryote", "Arctic")

full.rcbray.df = rbind(pro.ant.rcbray.df, pro.ant.abun.rcbray.df, pro.ant.int.rcbray.df, pro.ant.rare.rcbray.df, 
                     pro.arc.rcbray.df, pro.arc.abun.rcbray.df, pro.arc.int.rcbray.df, pro.arc.rare.rcbray.df,
                     euk.ant.rcbray.df, euk.ant.abun.rcbray.df, euk.ant.int.rcbray.df, euk.ant.rare.rcbray.df, 
                     euk.arc.rcbray.df,euk.arc.abun.rcbray.df, euk.arc.int.rcbray.df, euk.arc.rare.rcbray.df)

write.csv(full.rcbray.df, "../results/RCbray-results-table.csv")