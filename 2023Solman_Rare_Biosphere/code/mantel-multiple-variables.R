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
pro <- readRDS("../results/16S-phylo-object-rarefied-var-trans.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied-var-trans.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")


#And data with un-transformed variables
# 2. Import data
# pro.nt <- readRDS("../results/16S-phylo-object-rarefied.rds") 
# euk.nt <- readRDS("../results/18S-phylo-object-rarefied.rds")
# 
# pro.ant.nt <- subset_samples(pro.nt, Pole=="Antarctic")
# pro.ant.nt = filter_taxa(pro.ant.nt, function(x) sum(x) >= 1, TRUE)
# euk.ant.nt <- subset_samples(euk.nt, Pole=="Antarctic")
# euk.ant.nt = filter_taxa(euk.ant.nt, function(x) sum(x) >= 1, TRUE)
# 
# pro.ant.abun.nt <- readRDS("../results/16S-phylo-object-rarefied-ant-abun.rds")
# pro.ant.int.nt <- readRDS("../results/16S-phylo-object-rarefied-ant-int.rds")
# pro.ant.rare.nt <- readRDS("../results/16S-phylo-object-rarefied-ant-rare.rds")
# euk.ant.abun.nt <- readRDS("../results/18S-phylo-object-rarefied-ant-abun.rds")
# euk.ant.int.nt <- readRDS("../results/18S-phylo-object-rarefied-ant-int.rds")
# euk.ant.rare.nt <- readRDS("../results/18S-phylo-object-rarefied-ant-rare.rds")


# Method

# Prior to analysis environmental variables were z-score transformed and variables with high multicollinearity removed. Correlation between community composition and environmental variables was investigated using mantel tests on Bray-Curtis dissimilarity matrices and Euclidean distance matrix of combined and single environmental variables and geographic coordinates (Standardized Euclidean Distance). Partial mantel tests were performed, controlling for the effect of geographic coordinates and environmental distances. Bray-Curtis dissimilarities and Euclidean distances were calculated using the vegdist function of the vegan package. As different sets of environmental variables were available for different samples two analyses were performed.** Non-hierarchical cluster analysis of Bray-Curtis dissimilarities was used to explore community groupings. Multiple regression was performed on Bray-Curtis dissimilarity matrices, Euclidean distance matrix of environmental variables and geographic distance matrix to test for a causal relationship between environmental and spatial factors. Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarities was performed to explore variation between community compositions among sites. Distance-based redundancy analysis (db-RDA) was performed to explore potentially non-linear relationships between environmental variables and community composition. db-RDA was performed by carrying out PCoA on the Bray-Curtis dissimilarity matrix. The eigenvectors obtained through PCoA were then used, along with environmental variables, to carry out redundancy analysis. Canonical correspondence analysis (CCA) was carried out to explore the relationship between community structure and explanatory variables. Variation partition analysis (VPA) was used to partition the explanatory power of environmental and physical factors. 
# 
# 
# Should I do it like this instead: Geographic distance matrices were generated using earth.dist function from the ‘fossil’ package (Vavrek, 2011). Correlation between geographic distance matrices and community dissimilarity matrices was tested for significance with Mantel tests using Spearman’s correlation with 9999 permutations. 
# 
# #Multicollinearity
# When testing the effect of multiple independent variables on a dependent variable we need to make sure that the independent variables are not highly correlated. If they are highly correlated the regression model can't give us accurate information about what/how the independent variables are affecting the dependent variable. I.e. the parameters of the model become indeterminate and the standard errors of the estimates become infinitely large.
# 
# We can detect multicollinearity in a model by measuring the variance inflation factor (VIF) of each variable - this is the correlation and strength of correlation between independent variables in a regression model. We want to remove variables with VIF above 20 (ideally 5 or 10 but this is not always possible) as this indicates very high multicollinearity.
# 
# When carrying out Mantel and partial mantel tests of individual variables it is appropriate to use only those variables with VIF < 20 indicated by CCA and RDA (https://www.sciencedirect.com/science/article/pii/S0304389420323621). These are the variables we include in VPA - and these Mantel and Partial Mantel tests are here to support or deny the VPA results so we just include those variables used in VPA (https://academic.oup.com/femsec/article/96/6/fiaa071/5822764).

# Results

# phylo = pro.ant
# which_df = 1

mantel_matrix_func <- function(phylo, which_df){
  
  #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
  phylo.rel = transform_sample_counts(phylo, function(x) x / sum(x) )
  
  #get community data
  comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
  
  #get environmental data
  # vars = data.frame(sample_data(phylo.rel)) #get metadata
  # keep = c("DistanceToSea", "Altitude", "WaterDepth", "WaterVolume", "TotalDepth", "SedimentVolume", "Rock.Water", "Area", "Ice.lid", "Mass","Temp","EC","pH", "DOC","DIC","C.N","Cl.age","DO", "H","HCO3",     "NO3", "DON","TDN","TDP","DOP", "Cl","SO4","Na","K","Mg","Ca","F","DRP")
  # chem = vars[,names(vars) %in% keep]
  # # twochem = meta_df(phylo.rel)
  # # chem = twochem[[which_df]]
  
    #env data
  vars = data.frame(sample_data(phylo.rel))
 if (which_df == 1){
    keep = c("Ca", "Mg", "K", "Na", "SO4", "Cl", "NO3", "EC", "pH", "Area", "TotalDepth")
  } else if (which_df == 2){
    keep = c("DistanceToSea", "Altitude", "HCO3")
  }
  var.trim <- vars[,names(vars) %in% keep]
  chem <- var.trim[complete.cases(var.trim),]
  
  #match community data samples
  comm.trim <- comm[rownames(comm) %in% rownames(chem),]
  
  #check for vif
  vif.out = data.frame(vif.cca(cca(comm.trim ~., chem)))
  #vif.cca(capscale(comm.trim ~ ., chem, dist="bray", add=TRUE)) #similar results to cca
  #keep only variables with VIF of 20 or less
  chem.keep = rownames(vif.out)[vif.out$vif.cca.cca.comm.trim......chem.. <= 20]
  chem2 = data.frame(chem[,names(chem) %in% chem.keep])
  #chem.out <- chem2[complete.cases(chem2),]
  names(chem2) = chem.keep
  rownames(chem2) = rownames(chem)
  chem = chem2

  #get distance data
  geo <- get_my_long_lat(phylo.rel)
  
  #make sure all the same samples are used 
  comm <- comm[rownames(comm) %in% rownames(chem),]
  geo <- geo[rownames(geo) %in% rownames(chem),]
  
  # 5. Get Bray-Curtis dissimilarity matrix for community data
  commdist<-vegdist(comm, method="bray")
  
  #6. Get get distance matrix for environmental factors
  chemdist<-vegdist(chem, method="euclidean")
  
  # 7. Get distance matrix for spatial factors
  geodist<-vegdist(geo, method="euclidean")
  
  # 8. Mantel test for correlation between community composition and environmental factors
  comm_env <- vegan::mantel(commdist, chemdist, method="spearman", permutations=999)
  comm_env
  
  # 9. Mantel test for correlation between community composition and spatial factors.
  comm_geo <- vegan::mantel(commdist, geodist, method="spearman", permutations=999)
  comm_geo
  
  # 10. Partial Mantel test for correlation between community composition and environmental factors, controlling for spatial factors
  comm_env_geo <- mantel.partial(commdist, chemdist, geodist, method="spearman", permutations=999)
  comm_env_geo
  
  # 11. Partial Mantel test for correlation between community composition and spatial factors, controlling for environmental factors
  comm_geo_env <- mantel.partial(commdist, geodist, chemdist, method="spearman", permutations=999)
  comm_geo_env
  
  res.list = list(rownames(chem), names(chem), which_df, comm_env, comm_geo, comm_env_geo, comm_geo_env)
  
  return(res.list)
  
  }

# mantel_matrix_func2 <- function(phylo, which_df){
#   
#   #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
#   phylo.rel = transform_sample_counts(phylo, function(x) x / sum(x) )
#   
#   #get community data
#   comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
#   
#   #Get metadata of samples
#   samp_data <- data.frame(sample_data(phylo.rel))
#   samp_data = samp_data[,10:(ncol(samp_data))]
#   chem  <- samp_data [rowSums(is.na(samp_data )) != ncol(samp_data), ]
# 
#   #get distance data
#   geo <- get_my_long_lat(phylo.rel)
#   
#   #make sure all the same samples are used 
#   comm <- comm[rownames(comm) %in% rownames(chem),]
#   geo <- geo[rownames(geo) %in% rownames(chem),]
#   
#   # 5. Get Bray-Curtis dissimilarity matrix for community data
#   commdist<-vegdist(comm, method="bray")
#   
#   # 6. Get get distance matrix for environmental factors
#   chemdist<-vegdist(chem, method="euclidean", na.rm = TRUE)
#   
#   # 7. Get distance matrix for spatial factors
#   geodist<-vegdist(geo, method="euclidean")
#   
#   # 8. Mantel test for correlation between community composition and environmental factors
#   comm_env <- mantel(commdist, chemdist, method="spearman", permutations=999, na.rm=TRUE)
#   comm_env
#   
#   # 9. Mantel test for correlation between community composition and spatial factors.
#   comm_geo <- mantel(commdist, geodist, method="spearman", permutations=999, na.rm=TRUE)
#   comm_geo
#   
#   # 10. Partial Mantel test for correlation between community composition and environmental factors, controlling for spatial factors
#   comm_env_geo <- mantel.partial(commdist, chemdist, geodist, method="spearman", permutations=999, na.rm=TRUE)
#   comm_env_geo
#   
#   # 11. Partial Mantel test for correlation between community composition and spatial factors, controlling for environmental factors
#   comm_geo_env <- mantel.partial(commdist, geodist, chemdist, method="spearman", permutations=999, na.rm=TRUE)
#   comm_geo_env
#   
#   res.list = list(rownames(chem), names(chem), comm_env, comm_geo, comm_env_geo, comm_geo_env)
#   
#   return(res.list)
#   
# }


pro1 = mantel_matrix_func(pro.ant, 1)
pro2 = mantel_matrix_func(pro.ant, 2)

euk1 = mantel_matrix_func(euk.ant, 1)
euk2 = mantel_matrix_func(euk.ant, 2)

pro.abun1 = mantel_matrix_func(pro.ant.abun, 1)
pro.abun2 = mantel_matrix_func(pro.ant.abun, 2)

pro.int1 = mantel_matrix_func(pro.ant.int, 1)
pro.int2 = mantel_matrix_func(pro.ant.int, 2)

pro.rare1 = mantel_matrix_func(pro.ant.rare, 1)
pro.rare2 = mantel_matrix_func(pro.ant.rare, 2)

euk.abun1 = mantel_matrix_func(euk.ant.abun, 1)
euk.abun2 = mantel_matrix_func(euk.ant.abun, 2)

euk.int1 = mantel_matrix_func(euk.ant.int, 1)
euk.int2 = mantel_matrix_func(euk.ant.int, 2)

euk.rare1 = mantel_matrix_func(euk.ant.rare, 1)
euk.rare2 = mantel_matrix_func(euk.ant.rare, 2)

#Results table
Gene=c(rep("16S", 8), rep("18S", 8))
Com = c(rep("Full", 2), rep("Abun",2), rep("Int",2), rep("Rare", 2))
Com = rep(Com,2)
Analysis=c(rep(c("1", "2"), 8))
EnvP=c(round(pro1[[4]]$signif, 4), round(pro2[[4]]$signif,4),
       round(pro.abun1[[4]]$signif, 4), round(pro.abun2[[4]]$signif, 4),
       round(pro.int1[[4]]$signif, 4), round(pro.int2[[4]]$signif, 4),
       round(pro.rare1[[4]]$signif, 4), round(pro.rare2[[4]]$signif, 4),
       round(euk1[[4]]$signif, 4), round(euk2[[4]]$signif, 4),
       round(euk.abun1[[4]]$signif, 4), round(euk.abun2[[4]]$signif, 4),
       round(euk.int1[[4]]$signif, 4), round(euk.int2[[4]]$signif, 4),
       round(euk.rare1[[4]]$signif, 4), round(euk.rare2[[4]]$signif, 4))
EnvR=c(round(pro1[[4]]$statistic, 2), round(pro2[[4]]$statistic, 2),
       round(pro.abun1[[4]]$statistic, 2), round(pro.abun2[[4]]$statistic, 2),
       round(pro.int1[[4]]$statistic, 2), round(pro.int2[[4]]$statistic, 2),
       round(pro.rare1[[4]]$statistic, 2), round(pro.rare2[[4]]$statistic, 2),
       round(euk1[[4]]$statistic, 2), round(euk2[[4]]$statistic, 2),
       round(euk.abun1[[4]]$statistic, 2), round(euk.abun2[[4]]$statistic, 2),
       round(euk.int1[[4]]$statistic, 2), round(euk.int2[[4]]$statistic, 2),
       round(euk.rare1[[4]]$statistic, 2), round(euk.rare2[[4]]$statistic, 2))
GeoP=c(round(pro1[[5]]$signif, 4), round(pro2[[5]]$signif, 4),
       round(pro.abun1[[5]]$signif, 4), round(pro.abun2[[5]]$signif, 4),
       round(pro.int1[[5]]$signif, 4), round(pro.int2[[5]]$signif, 4),
       round(pro.rare1[[5]]$signif, 4), round(pro.rare2[[5]]$signif, 4),
       round(euk1[[5]]$signif, 4), round(euk2[[5]]$signif, 4),
       round(euk.abun1[[5]]$signif, 4), round(euk.abun2[[5]]$signif, 4),
       round(euk.int1[[5]]$signif, 4), round(euk.int2[[5]]$signif, 4),
       round(euk.rare1[[5]]$signif, 4), round(euk.rare2[[5]]$signif, 4))
GeoR=c(round(pro1[[5]]$statistic, 2), round(pro2[[5]]$statistic, 2),
       round(pro.abun1[[5]]$statistic, 2), round(pro.abun2[[5]]$statistic, 2),
       round(pro.int1[[5]]$statistic, 2), round(pro.int2[[5]]$statistic, 2),
       round(pro.rare1[[5]]$statistic, 2), round(pro.rare2[[5]]$statistic, 2),
       round(euk1[[5]]$statistic, 2), round(euk2[[5]]$statistic, 2),
       round(euk.abun1[[5]]$statistic, 2), round(euk.abun2[[5]]$statistic, 2),
       round(euk.int1[[5]]$statistic, 2), round(euk.int2[[5]]$statistic, 2),
       round(euk.rare1[[5]]$statistic, 2), round(euk.rare2[[5]]$statistic, 2))
Env_GeoP=c(round(pro1[[6]]$signif, 4), round(pro2[[6]]$signif, 4),
       round(pro.abun1[[6]]$signif, 4), round(pro.abun2[[6]]$signif, 4),
       round(pro.int1[[6]]$signif, 4), round(pro.int2[[6]]$signif, 4),
       round(pro.rare1[[6]]$signif, 4), round(pro.rare2[[6]]$signif, 4),
       round(euk1[[6]]$signif, 4), round(euk2[[6]]$signif, 4),
       round(euk.abun1[[6]]$signif, 4), round(euk.abun2[[6]]$signif, 4),
       round(euk.int1[[6]]$signif, 4), round(euk.int2[[6]]$signif, 4),
       round(euk.rare1[[6]]$signif, 4), round(euk.rare2[[6]]$signif, 4))
Env_GeoR=c(round(pro1[[6]]$statistic, 2), round(pro2[[6]]$statistic, 2),
       round(pro.abun1[[6]]$statistic, 2), round(pro.abun2[[6]]$statistic, 2),
       round(pro.int1[[6]]$statistic, 2), round(pro.int2[[6]]$statistic, 2),
       round(pro.rare1[[6]]$statistic, 2), round(pro.rare2[[6]]$statistic, 2),
       round(euk1[[6]]$statistic, 2), round(euk2[[6]]$statistic, 2),
       round(euk.abun1[[6]]$statistic, 2), round(euk.abun2[[6]]$statistic, 2),
       round(euk.int1[[6]]$statistic, 2), round(euk.int2[[6]]$statistic, 2),
       round(euk.rare1[[6]]$statistic, 2), round(euk.rare2[[6]]$statistic, 2))
Geo_EnvP=c(round(pro1[[7]]$signif, 4), round(pro2[[7]]$signif, 4),
       round(pro.abun1[[7]]$signif, 4), round(pro.abun2[[7]]$signif, 4),
       round(pro.int1[[7]]$signif, 4), round(pro.int2[[7]]$signif, 4),
       round(pro.rare1[[7]]$signif, 4), round(pro.rare2[[7]]$signif, 4),
       round(euk1[[7]]$signif, 4), round(euk2[[7]]$signif, 4),
       round(euk.abun1[[7]]$signif, 4), round(euk.abun2[[7]]$signif, 4),
       round(euk.int1[[7]]$signif, 4), round(euk.int2[[7]]$signif, 4),
       round(euk.rare1[[7]]$signif, 4), round(euk.rare2[[7]]$signif, 4))
Geo_EnvR=c(round(pro1[[7]]$statistic, 2), round(pro2[[7]]$statistic, 2),
       round(pro.abun1[[7]]$statistic, 2), round(pro.abun2[[7]]$statistic, 2),
       round(pro.int1[[7]]$statistic, 2), round(pro.int2[[7]]$statistic, 2),
       round(pro.rare1[[7]]$statistic, 2), round(pro.rare2[[7]]$statistic, 2),
       round(euk1[[7]]$statistic, 2), round(euk2[[7]]$statistic, 2),
       round(euk.abun1[[7]]$statistic, 2), round(euk.abun2[[7]]$statistic, 2),
       round(euk.int1[[7]]$statistic, 2), round(euk.int2[[7]]$statistic, 2),
       round(euk.rare1[[7]]$statistic, 2), round(euk.rare2[[7]]$statistic, 2))


res.df = cbind(Gene, Com, Analysis,  EnvR, EnvP, GeoR, GeoP,  Env_GeoR, Env_GeoP, Geo_EnvR, Geo_EnvP)
res.df = as.data.frame(res.df)
names(res.df) = c("Gene", "Comm", "Num", "Env rho", "Env P", "Geo rho", "Geo P", "Env - Geo rho",  "Env - Geo P", "Geo - Env rho", "Geo - Env P")

write.csv(res.df, "../results/Mantel-Partial-Mantel-results.csv")

# out = pro1[[2]][1]
# for (i in 2:length(pro1[[2]])){
#   out = paste(out, pro1[[2]][i], sep=", ")
# }
# 
# out2 = pro2[[2]][1]
# for (i in 2:length(pro2[[2]])){
#   out2 = paste(out2, pro2[[2]][i], sep=", ")
# }

# pro = mantel_matrix_func2(pro.ant)
# euk = mantel_matrix_func2(euk.ant)
# pro.abun = mantel_matrix_func2(pro.ant.abun)
# pro.int = mantel_matrix_func2(pro.ant.int)
# pro.rare = mantel_matrix_func2(pro.ant.rare)
# euk.abun = mantel_matrix_func2(euk.ant.abun)
# euk.int = mantel_matrix_func2(euk.ant.int)
# euk.rare = mantel_matrix_func2(euk.ant.rare)
# 
# #Results table
# Gene=c(rep("16S", 4), rep("18S", 4))
# Com = c("Full","Abun","Int","Rare")
# Com = rep(Com,2)
# #Analysis=c(rep(c("1", "2"), 8))
# EnvP=c(round(pro[[3]]$signif, 4),
#        round(pro.abun[[3]]$signif, 4),
#        round(pro.int[[3]]$signif, 4),
#        round(pro.rare[[3]]$signif, 4),
#        round(euk[[3]]$signif, 4),
#        round(euk.abun[[3]]$signif, 4),
#        round(euk.int[[3]]$signif, 4),
#        round(euk.rare[[3]]$signif, 4))
# 
# EnvR=c(round(pro1[[4]]$statistic, 2), round(pro2[[4]]$statistic, 2),
#        round(pro.abun1[[4]]$statistic, 2), round(pro.abun2[[4]]$statistic, 2),
#        round(pro.int1[[4]]$statistic, 2), round(pro.int2[[4]]$statistic, 2),
#        round(pro.rare1[[4]]$statistic, 2), round(pro.rare2[[4]]$statistic, 2),
#        round(euk1[[4]]$statistic, 2), round(euk2[[4]]$statistic, 2),
#        round(euk.abun1[[4]]$statistic, 2), round(euk.abun2[[4]]$statistic, 2),
#        round(euk.int1[[4]]$statistic, 2), round(euk.int2[[4]]$statistic, 2),
#        round(euk.rare1[[4]]$statistic, 2), round(euk.rare2[[4]]$statistic, 2))
# GeoP=c(mod_my_p_val(pro1[[5]]$signif), mod_my_p_val(pro2[[5]]$signif),
#        mod_my_p_val(pro.abun1[[5]]$signif), mod_my_p_val(pro.abun2[[5]]$signif),
#        mod_my_p_val(pro.int1[[5]]$signif), mod_my_p_val(pro.int2[[5]]$signif),
#        mod_my_p_val(pro.rare1[[5]]$signif), mod_my_p_val(pro.rare2[[5]]$signif),
#        mod_my_p_val(euk1[[5]]$signif), mod_my_p_val(euk2[[5]]$signif),
#        mod_my_p_val(euk.abun1[[5]]$signif), mod_my_p_val(euk.abun2[[5]]$signif),
#        mod_my_p_val(euk.int1[[5]]$signif), mod_my_p_val(euk.int2[[5]]$signif),
#        mod_my_p_val(euk.rare1[[5]]$signif), mod_my_p_val(euk.rare2[[5]]$signif))
# GeoR=c(round(pro1[[5]]$statistic, 2), round(pro2[[5]]$statistic, 2),
#        round(pro.abun1[[5]]$statistic, 2), round(pro.abun2[[5]]$statistic, 2),
#        round(pro.int1[[5]]$statistic, 2), round(pro.int2[[5]]$statistic, 2),
#        round(pro.rare1[[5]]$statistic, 2), round(pro.rare2[[5]]$statistic, 2),
#        round(euk1[[5]]$statistic, 2), round(euk2[[5]]$statistic, 2),
#        round(euk.abun1[[5]]$statistic, 2), round(euk.abun2[[5]]$statistic, 2),
#        round(euk.int1[[5]]$statistic, 2), round(euk.int2[[5]]$statistic, 2),
#        round(euk.rare1[[5]]$statistic, 2), round(euk.rare2[[5]]$statistic, 2))
# Env_GeoP=c(mod_my_p_val(pro1[[6]]$signif), mod_my_p_val(pro2[[6]]$signif),
#        mod_my_p_val(pro.abun1[[6]]$signif), mod_my_p_val(pro.abun2[[6]]$signif),
#        mod_my_p_val(pro.int1[[6]]$signif), mod_my_p_val(pro.int2[[6]]$signif),
#        mod_my_p_val(pro.rare1[[6]]$signif), mod_my_p_val(pro.rare2[[6]]$signif),
#        mod_my_p_val(euk1[[6]]$signif), mod_my_p_val(euk2[[6]]$signif),
#        mod_my_p_val(euk.abun1[[6]]$signif), mod_my_p_val(euk.abun2[[6]]$signif),
#        mod_my_p_val(euk.int1[[6]]$signif), mod_my_p_val(euk.int2[[6]]$signif),
#        mod_my_p_val(euk.rare1[[6]]$signif), mod_my_p_val(euk.rare2[[6]]$signif))
# Env_GeoR=c(round(pro1[[6]]$statistic, 2), round(pro2[[6]]$statistic, 2),
#        round(pro.abun1[[6]]$statistic, 2), round(pro.abun2[[6]]$statistic, 2),
#        round(pro.int1[[6]]$statistic, 2), round(pro.int2[[6]]$statistic, 2),
#        round(pro.rare1[[6]]$statistic, 2), round(pro.rare2[[6]]$statistic, 2),
#        round(euk1[[6]]$statistic, 2), round(euk2[[6]]$statistic, 2),
#        round(euk.abun1[[6]]$statistic, 2), round(euk.abun2[[6]]$statistic, 2),
#        round(euk.int1[[6]]$statistic, 2), round(euk.int2[[6]]$statistic, 2),
#        round(euk.rare1[[6]]$statistic, 2), round(euk.rare2[[6]]$statistic, 2))
# Geo_EnvP=c(mod_my_p_val(pro1[[7]]$signif), mod_my_p_val(pro2[[7]]$signif),
#        mod_my_p_val(pro.abun1[[7]]$signif), mod_my_p_val(pro.abun2[[7]]$signif),
#        mod_my_p_val(pro.int1[[7]]$signif), mod_my_p_val(pro.int2[[7]]$signif),
#        mod_my_p_val(pro.rare1[[7]]$signif), mod_my_p_val(pro.rare2[[7]]$signif),
#        mod_my_p_val(euk1[[7]]$signif), mod_my_p_val(euk2[[7]]$signif),
#        mod_my_p_val(euk.abun1[[7]]$signif), mod_my_p_val(euk.abun2[[7]]$signif),
#        mod_my_p_val(euk.int1[[7]]$signif), mod_my_p_val(euk.int2[[7]]$signif),
#        mod_my_p_val(euk.rare1[[7]]$signif), mod_my_p_val(euk.rare2[[7]]$signif))
# Geo_EnvR=c(round(pro1[[7]]$statistic, 2), round(pro2[[7]]$statistic, 2),
#        round(pro.abun1[[7]]$statistic, 2), round(pro.abun2[[7]]$statistic, 2),
#        round(pro.int1[[7]]$statistic, 2), round(pro.int2[[7]]$statistic, 2),
#        round(pro.rare1[[7]]$statistic, 2), round(pro.rare2[[7]]$statistic, 2),
#        round(euk1[[7]]$statistic, 2), round(euk2[[7]]$statistic, 2),
#        round(euk.abun1[[7]]$statistic, 2), round(euk.abun2[[7]]$statistic, 2),
#        round(euk.int1[[7]]$statistic, 2), round(euk.int2[[7]]$statistic, 2),
#        round(euk.rare1[[7]]$statistic, 2), round(euk.rare2[[7]]$statistic, 2))
# 
# res.df = cbind(Gene, Com, Analysis, EnvP, EnvR, GeoP, GeoR, Env_GeoP, Env_GeoR, Geo_EnvP, Geo_EnvR)
# res.df = as.data.frame(res.df)
# names(res.df) = c("Gene", "Comm", "Num", "Env P", "Env R2", "Geo P", "Geo R2", "Env - Geo P", "Env - Geo R2", "Geo - Env P", "Geo - Env R2")
# 
# out = pro1[[2]][1]
# for (i in 2:length(pro1[[2]])){
#   out = paste(out, pro1[[2]][i], sep=", ")
# }
# 
# out2 = pro2[[2]][1]
# for (i in 2:length(pro2[[2]])){
#   out2 = paste(out2, pro2[[2]][i], sep=", ")
# }
