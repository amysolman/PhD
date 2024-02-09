# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. z-transform variables
# 4. Remove variables with high collinearity
# 5. Export data

# Method

#*Prior to analysis environmental variables were z-score transformed and variables with high multicollinearity removed.** Correlation between community composition and environmental variables was investigated using mantel tests on Bray-Curtis dissimilarity matrices and Euclidean distance matrix of environmental variables (Standardized Euclidean Distance). Partial mantel tests were performed, controlling for the effect of a geographic distance matrix. Non-hierarchical cluster analysis of Bray-Curtis dissimilarities was used to explore community groupings. Multiple regression was performed on Bray-Curtis dissimilarity matrices, Euclidean distance matrix of environmental variables and geographic distance matrix to test for a causal relationship between environmental and spatial factors. Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarities was performed to explore variation between community compositions among sites. Distance-based redundancy analysis (db-RDA) was performed to explore potentially non-linear relationships between environmental variables and community composition. db-RDA was performed by carrying out PCoA on the Bray-Curtis dissimilarity matrix. The eigenvectors obtained through PCoA were then used, along with environmental variables, to carry out redundancy analysis. Canonical correspondence analysis (CCA) was carried out to explore the relationship between community structure and explanatory variables. Variation partition analysis (VPA) was used to partition the explanatory power of environmental and physical factors. 

rm(list=ls())
graphics.off()

library(phyloseq)
library(dplyr) #coalesce function
#library(corrplot) #corrplot function

# 2. Import data
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


var_trans <- function(phylo){
  
  #get the dataframe
  df = data.frame(sample_data(phylo))
  
  #columns to be z-transformed
  keep = c("Latitude", "DistanceToSea", "Altitude", "WaterDepth", "WaterVolume", "TotalDepth", "SedimentDepth", "SedimentVolume", "MedianGrainSize", "Rock.Water", "Area", "Ice.lid", "Mass","Temp","EC","pH","pCO2", "DOC","DIC","C.N","Cl.age","DO","N","C","H","HCO3",  "NH4", "TOxN",            "NO2","NO3","TN","TIN","TON","DON","TDN","TDP","DOP","PO4","SiO2","Cl","SO4","Na","K","Mg","Ca","F","DRP")
  x = df[,names(df) %in% keep]
  
  names(x)[1] = "Latitude.Z" 
  
  #make sure data is numeric
  x.num = as.data.frame(lapply(x,as.numeric))
  names(x.num) = names(x)
  
  #z-score transform the data
  df.z = data.frame(scale(x.num))
  names(df.z) = names(x)
  
  #add additional info for metadata
  final.to.export = cbind(df$SampleID, df$Glacier, df$Region, df$Pole, df$Latitude, df$Longitude, df.z)
  
  colnames(final.to.export) <- c("SampleID", "Glacier", "Region", "Pole", "Latitude", "Longitude", colnames(df.z))
  
  rownames(final.to.export) = final.to.export$SampleID
  
  new.phylo = phylo
  
  sample_data(new.phylo) = final.to.export
  
  return(new.phylo)
  
  
}

pro.res = var_trans(pro)
euk.res = var_trans(euk)
pro.abun.res = var_trans(pro.ant.abun)
pro.int.res = var_trans(pro.ant.int)
pro.rare.res = var_trans(pro.ant.rare)
euk.abun.res = var_trans(euk.ant.abun)
euk.int.res = var_trans(euk.ant.int)
euk.rare.res = var_trans(euk.ant.rare)

#export new phyloseq objects
saveRDS(pro.res, "../results/16S-phylo-object-var-trans.rds")
saveRDS(euk.res, "../results/18S-phylo-object-var-trans.rds")
saveRDS(pro.abun.res, "../results/16S-phylo-object-ant-abun-var-trans.rds")
saveRDS(pro.int.res, "../results/16S-phylo-object-ant-int-var-trans.rds")
saveRDS(pro.rare.res, "../results/16S-phylo-object-ant-rare-var-trans.rds")
saveRDS(euk.abun.res, "../results/18S-phylo-object-ant-abun-var-trans.rds")
saveRDS(euk.int.res, "../results/18S-phylo-object-ant-int-var-trans.rds")
saveRDS(euk.rare.res, "../results/18S-phylo-object-ant-rare-var-trans.rds")
