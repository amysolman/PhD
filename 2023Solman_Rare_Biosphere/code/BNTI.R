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


#Compute Beta-nearest Taxon Index (βNTI)

#Define function for calculating the βMNTD for each random null community

# Function for calculating the βMNTD for each random null community
bMNTD_null_func <- function(i, OTU.table, tree){
  tree$tip.label = sample(tree$tip.label)
  bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
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

# The main function for calculating βNTI
Phylo_turnover <- function(physeq, reps, nproc){
  # Extract OTU table
  OTU.table = t(otu_table(physeq))
  # Extract phylogenetic tree
  tree = phy_tree(physeq)
  # Get βMNTD between all communities
  bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  A <- attr(bMNTD_o, "Size")
  B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
  if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
  bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_o))
  
  # Get βMNTD for randomized null communities
  rep.list = seq(1, reps)
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

##RUN THE MODEL

get_bnti_res <- function(phylo, subcom, group, pole){
  
  #remove taxa with zero counts
  phylo = filter_taxa(phylo, function(x) sum(x) > 0, TRUE)
  #remove samples with zero counts
  phylo = prune_samples(sample_sums(phylo) > 0, phylo)
  
  df = Phylo_turnover(phylo, 1000, 10) #should be 1000 10
  df$Group = group
  df$Subcommunity = subcom
  df$Pole = pole
  return(df)
}

pro.ant.bNTI.df = get_bnti_res(pro.ant, "Full", "Prokaryote", "Antarctic")
pro.ant.abun.bNTI.df = get_bnti_res(pro.ant.abun, "Abundant", "Prokaryote", "Antarctic")
pro.ant.int.bNTI.df = get_bnti_res(pro.ant.int, "Intermediate", "Prokaryote", "Antarctic")
pro.ant.rare.bNTI.df = get_bnti_res(pro.ant.rare, "Rare", "Prokaryote", "Antarctic")
pro.arc.bNTI.df = get_bnti_res(pro.arc, "Full", "Prokaryote", "Arctic")
pro.arc.abun.bNTI.df = get_bnti_res(pro.arc.abun, "Abundant", "Prokaryote", "Arctic")
pro.arc.int.bNTI.df = get_bnti_res(pro.arc.int, "Intermediate", "Prokaryote", "Arctic")
pro.arc.rare.bNTI.df = get_bnti_res(pro.arc.rare, "Rare", "Prokaryote", "Arctic")
euk.ant.bNTI.df = get_bnti_res(euk.ant, "Full", "Eukaryote", "Antarctic")
euk.ant.abun.bNTI.df = get_bnti_res(euk.ant.abun, "Abundant", "Eukaryote", "Antarctic")
euk.ant.int.bNTI.df = get_bnti_res(euk.ant.int, "Intermediate", "Eukaryote", "Antarctic")
euk.ant.rare.bNTI.df = get_bnti_res(euk.ant.rare, "Rare", "Eukaryote", "Antarctic")
euk.arc.bNTI.df = get_bnti_res(euk.arc, "Full", "Eukaryote", "Arctic")
euk.arc.abun.bNTI.df = get_bnti_res(euk.arc.abun, "Abundant", "Eukaryote", "Arctic")
euk.arc.int.bNTI.df = get_bnti_res(euk.arc.int, "Intermediate", "Eukaryote", "Arctic")
euk.arc.rare.bNTI.df = get_bnti_res(euk.arc.rare, "Rare", "Eukaryote", "Arctic")

full.bNTI.df = rbind(pro.ant.bNTI.df, pro.ant.abun.bNTI.df, pro.ant.int.bNTI.df, pro.ant.rare.bNTI.df, 
                     pro.arc.bNTI.df, pro.arc.abun.bNTI.df, pro.arc.int.bNTI.df, pro.arc.rare.bNTI.df,
                     euk.ant.bNTI.df, euk.ant.abun.bNTI.df, euk.ant.int.bNTI.df, euk.ant.rare.bNTI.df, 
                     euk.arc.bNTI.df,euk.arc.abun.bNTI.df, euk.arc.int.bNTI.df, euk.arc.rare.bNTI.df)

write.csv(full.bNTI.df, "../results/bNTI-results-table.csv")
