# **Script breakdown**

# 1. Plot our sample coordinates
# 2. Calculate PCNMs from a Euclidean distance matrix of sample coordinates and extract scores associated with these new variables.
# 3. Plot the first 6 PCNM axis against geographic position
# 4. Choose which variables are chemical and which are physical
# 5. Partition variation among three predictor tables and run varpart() 
# 6. Plot model
# 7. Use ordistep to select significant PCNM axes
# 8. Use ordistep to select significant chemical variables
# 9. Use ordistep to select significant physical variables
# 10. Partition variation among three predictor tables and run varpart() 
# 11. Plot model
# 12. Test significance of partition 1
# 13. Test significance of partition 2
# 14. Test significance of partition 3
# 15. Test significance of partition 1 controlling for partition 2 and 3
# 16. Test significance of partition 2 controlling for partition 1 and 3 
# 17. Test significance of partition 3 controlling for partition 1 and 2 
# 
# Report: Terms in first model, plot of first model,
#         Terms in final model, plot of final model,
#         significance of partition 1, 2 and 3,
#         significance of each partition, controlling for other partitions.


# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(vegan) #for diversity indices 
# library(dplyr) #for coalesce function
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

# 
# 
# # Method
# 
# Prior to analysis environmental variables were z-score transformed and variables with high multicollinearity removed. Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarities was performed to explore variation between community compositions among sites. Correlation between community composition and environmental variables was investigated using mantel tests on Bray-Curtis dissimilarity matrices and Euclidean distance matrix of environmental variables and geographic coordinates (Standardized Euclidean Distance). Partial mantel tests were performed, controlling for the effect of geographic coordinates and environmental distances. Bray-Curtis dissimilarities and Euclidean distances were calculated using the vegdist function of the vegan package. As different sets of environmental variables were available for different samples two analyses were performed. Canonical correspondence analysis (CCA) was carried out to explore the relationship between community structure and explanatory variables. Distance-based redundancy analysis (db-RDA) was performed to explore linear relationships between environmental/ geographic variables and community composition. db-RDA was performed by carrying out PCoA on the Bray-Curtis dissimilarity matrix. The eigenvectors obtained through PCoA were then used, along with environmental/geographic variables, to carry out redundancy analysis. Forward variable selection was carried out to reduce models to their significant terms. db-RDA was carried out using the capscale function in the vegan package. **Variation partition analysis (VPA) was used to partition the explanatory power of environmental, physical and spatial factors on Hellinger-transformed community matrices.**
#   
#   ## Variation Paritioning Analysis Explained (VPA) 
#   
#   VPA attempts to partition the explanatory power of difference explanatory matrices or groups (e.g. environmental/chemical/physical/geographic) in relation to the same response matrix (community dissimilarity matrix). VPA can be used with RDA or CCA. When using RDA multiple partial RDAs will be run to determine the partial, linear effect of each explanatory matrix on the response matrix. With CCA the total inertia of the response matrix is partitioned.
# 
# The shared partition is not an interaction term and cannot be assigned significance. It illustrates that the explanatory matrices are redundant in this partition and highlight multicollinearity in the model. Bigger shared partition = more collinearity! Evidence of the redundancy (collinearity) of the explanatory matrices can also be tested using mantel test. If they are collinear they are likely to share explanatory power and have a large partition B. These redundant variables should be removed from further investigation.
# 
# Using VPA with CCA is not universally accepted as the inertia partitions from CCA are not truly comparable - partial CCA may be a more appropriate method for this.
# 
# We can incorporating information on the physical locations of our samples in a slightly different way. The patterns linked to spatial distribution can be associated with unmeasured (and spatially autocorrelated) environmental variation and species dispersal among each location. We use the function pcnm() in the vegan package to generate a dataframe of variables representing different spatial scales. 

# Results


# phylo = pro.ant.abun
# which_df = 1


simple_VPA <- function(phylo, which_df){
  
  set.seed(666) # so we get the same results for each ANOVA-like permutation 
  
  #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
  phylo.rel = transform_sample_counts(phylo, function(x) x / sum(x) )
  
  #get community data
  comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
  
  #Hellinger transform the data to correct for double zero problem
  #comm.hel <- decostand(comm, method = "hellinger") #converts communities to relative abundances and square roots them
  
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
  vif.out = data.frame(vif.cca(rda(comm.trim ~., chem)))
  #vif.cca(capscale(comm.trim ~ ., chem, dist="bray", add=TRUE)) #similar results to cca
  #keep only variables with VIF of 20 or less
  chem.keep = rownames(vif.out)[vif.out$vif.cca.rda.comm.trim......chem.. <= 20]
  chem2 = data.frame(chem[,names(chem) %in% chem.keep])
  #chem.out <- chem2[complete.cases(chem2),]
  names(chem2) = chem.keep
  rownames(chem2) = rownames(chem)
  chem = chem2
  
  #get distance data
  geo <- get_my_long_lat(phylo)
  
  #make sure all the same samples are used 
  comm.trim <- comm[rownames(comm) %in% rownames(chem),]
  geo.trim <- geo[rownames(geo) %in% rownames(chem),]
  
  #Calculate PCNMs from a Euclidean distance matrix of sample coordinates and extract scores associated with these new variables.
  geo.pcnm <- as.data.frame(scores(pcnm(dist(geo.trim))))
  
  # Use ordistep to select significant PCNM axes
  fwd.sel.geo <- ordistep(rda(comm.trim ~ 1, geo.pcnm, scale=T),
                          scope = formula(rda(comm.trim ~ ., geo.pcnm, scale=T)))
  
  axis.keep = rownames(fwd.sel.geo$anova)
  axis.keep2 = sub('..', '', axis.keep)
  geo.pcnm.sel = geo.pcnm[,names(geo.pcnm) %in% axis.keep2]
  
  if (length(axis.keep2) > 0){
    
    comm.var <- varpart(comm.trim, ~ ., geo.pcnm.sel, data=chem)
    
  } else {
    
    comm.var = "No significant spatial axis"
  }
  
  
  
  return(comm.var)
}



out1 = simple_VPA(pro.ant.abun, 1)
out2 = simple_VPA(pro.ant.abun, 2)
out3 = simple_VPA(pro.ant.int, 1)
out4 = simple_VPA(pro.ant.int, 2)
out5 = simple_VPA(pro.ant.rare, 1)
out6 = simple_VPA(pro.ant.rare, 2)

#For adding label to base plot
line = -2
cex = 2
adj  = 0.025

pdf("../results/VPA-results-prokaryote.pdf", width=10, height=11)
par(mfrow=c(3,2),
    mar=c(1,1, 1, 1))
plot(out1, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Abundant 1",cex.main=cex,col="black",font=2,line=line)
plot(out2, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Abundant 2",cex.main=cex,col="black",font=2,line=line)
plot(out3, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Intermediate 1",cex.main=cex,col="black",font=2,line=line)
plot(out4, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Intermediate 2",cex.main=cex,col="black",font=2,line=line)
plot(out5, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Rare 1",cex.main=cex,col="black",font=2,line=line)
plot(out6, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Rare 2",cex.main=cex,col="black",font=2,line=line)
dev.off()


out1 = simple_VPA(euk.ant.abun, 1)
out2 = simple_VPA(euk.ant.abun, 2)
out3 = simple_VPA(euk.ant.int, 1)
out4 = simple_VPA(euk.ant.int, 2)
out5 = simple_VPA(euk.ant.rare, 1)
out6 = simple_VPA(euk.ant.rare, 2)

#For adding label to base plot
line = -2
cex = 2
adj  = 0.025

pdf("../results/VPA-results-eukaryote.pdf", width=10, height=11)
par(mfrow=c(3,2),
    mar=c(1,1, 1, 1))
plot(out1, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Abundant 1",cex.main=cex,col="black",font=2,line=line)
plot(out2, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Abundant 2",cex.main=cex,col="black",font=2,line=line)
plot(out3, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Intermediate 1",cex.main=cex,col="black",font=2,line=line)
plot(out4, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Intermediate 2",cex.main=cex,col="black",font=2,line=line)
plot(out5, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Rare 1",cex.main=cex,col="black",font=2,line=line)
plot(out6, bg=1:2, Xnames=c('Environmental', 'Spatial'), id.size=0.75)
title(outer=outer,adj=adj,main="Rare 2",cex.main=cex,col="black",font=2,line=line)
dev.off()