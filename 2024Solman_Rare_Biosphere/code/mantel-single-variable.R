# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Get data editing function
# 4. Define function for getting physical data
# 5. Get Bray-Curtis dissimilarity matrix for community data
# 6. Get get distance matrix for environmental factors
# 7. Get distance matrix for spatial factors
# 8. Mantel test for correlation between community composition and environmental factors
# 9. Mantel test for correlation between community composition and spatial factors.
# 10. Partial Mantel test for correlation between community composition and environmental factors, controlling for spatial factors
# 11. Partial Mantel test for correlation between community composition and spatial factors, controlling for environmental factors

# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(vegan) #for diversity indices 
library(dplyr) #for coalesce function
library(stringr) #for subsetting strings in lat long function
# #install.packages("corrplot")
# library(corrplot)
library(phyloseq)
source("00-solman-functions.R")

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) >= 1, TRUE)

#Antarctic prokaryotes
pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun.rds") 
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int.rds") 
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare.rds") 

#Antarctic eukaryotes
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun.rds") 
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int.rds") 
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare.rds") 

# Method

#Prior to analysis environmental variables were z-score transformed and variables with high multicollinearity removed. Correlation between community composition and environmental variables was investigated using mantel tests on Bray-Curtis dissimilarity matrices and Euclidean distance matrix of combined and single environmental variables and geographic coordinates (Standardized Euclidean Distance). Partial mantel tests were performed, controlling for the effect of geographic coordinates and environmental distances. Bray-Curtis dissimilarities and Euclidean distances were calculated using the vegdist function of the vegan package. As different sets of environmental variables were available for different samples two analyses were performed.** Non-hierarchical cluster analysis of Bray-Curtis dissimilarities was used to explore community groupings. Multiple regression was performed on Bray-Curtis dissimilarity matrices, Euclidean distance matrix of environmental variables and geographic distance matrix to test for a causal relationship between environmental and spatial factors. Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarities was performed to explore variation between community compositions among sites. Distance-based redundancy analysis (db-RDA) was performed to explore potentially non-linear relationships between environmental variables and community composition. db-RDA was performed by carrying out PCoA on the Bray-Curtis dissimilarity matrix. The eigenvectors obtained through PCoA were then used, along with environmental variables, to carry out redundancy analysis. Canonical correspondence analysis (CCA) was carried out to explore the relationship between community structure and explanatory variables. Variation partition analysis (VPA) was used to partition the explanatory power of environmental and physical factors. 


#Should I do it like this instead: Geographic distance matrices were generated using earth.dist function from the ‘fossil’ package (Vavrek, 2011). Correlation between geographic distance matrices and community dissimilarity matrices was tested for significance with Mantel tests using Spearman’s correlation with 9999 permutations. 

# Results

#Simple Mantel tests on single variables

#Compare Bray-Curtis dissimilarity matrix with Euclidean distance matrix of each variable of interest

#phylo1 = pro.ant.abun.nt #we'll use the non-transformed variables for simple mantel tests
# phylo = euk.ant.int
# subcom = "Intermediate"


simple_mantel <- function(phylo, subcom){
  
  #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
  #remove taxa or samples with zero counts
  ps = filter_taxa(phylo, function(x) sum(x) >0, TRUE)
  ps = prune_samples(sample_sums(ps)>0, ps)
  phylo.rel = transform_sample_counts(ps, function(x) x / sum(x) )
  
  #get community data
  comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
  
  #get environmental data
  vars = data.frame(sample_data(phylo.rel)) #get metadata
  keep = c("Latitude", "DistanceToSea", "Altitude", "WaterDepth", "WaterVolume", "TotalDepth", "SedimentDepth", "SedimentVolume", "MedianGrainSize", "Rock.Water", "Area", "Ice.lid", "Mass","Temp","EC","pH","pCO2", "DOC","DIC","C.N","Cl.age","DO","N","C","H","HCO3",  "NH4", "TOxN",            "NO2","NO3","TN","TIN","TON","DON","TDN","TDP","DOP","PO4","SiO2","Cl","SO4","Na","K","Mg","Ca","F","DRP")
  vars.keep = vars[,names(vars) %in% keep]
  
  #save output
  variname = vector()
  rho = vector()
  pval = vector()
  
  for (i in 1:ncol(vars.keep)){
    
    # Get Bray-Curtis dissimilarity matrix for community data
    commdist<-vegdist(comm, method="bray")
    
    vardist<-vegdist(vars.keep[,i], method="euclidean", na.rm = TRUE)
    
    set.seed(666)
    
    comm_env <- vegan::mantel(commdist, vardist, method="spearman", permutations=999, na.rm=TRUE)
    comm_env
    
    variname = c(variname, names(vars.keep)[[i]])
    rho = c(rho, round(comm_env$statistic, 2))
    pval = c(pval, round(comm_env$signif, 4))
    
  }
  
  res.df = data.frame(Variable = variname, Rho = rho, Pval = pval)
  names(res.df) = c("Variable", paste0(subcom, " MantelR"), paste0(subcom, " P val"))
  
  return(res.df)
  
}

pro.abun.sm = simple_mantel(pro.ant.abun, "Abundant")
pro.int.sm = simple_mantel(pro.ant.int, "Intermediate")
pro.rare.sm = simple_mantel(pro.ant.rare, "Rare")
euk.abun.sm = simple_mantel(euk.ant.abun, "Abundant")
euk.int.sm = simple_mantel(euk.ant.int, "Intermediate")
euk.rare.sm = simple_mantel(euk.ant.rare, "Rare")

#Results
pro.final = cbind(pro.abun.sm, pro.int.sm[, 2:3], pro.rare.sm[,2:3])
euk.final = cbind(euk.abun.sm, euk.int.sm[,2:3], euk.rare.sm[,2:3])

export_res = rbind(pro.final, euk.final)
export_res = cbind(Group=c(rep("Prokaryote", nrow(pro.final)), rep("Eukaryote", nrow(euk.final))), export_res)

write.csv(export_res, "../results/mantel-single-variable-results.csv")
