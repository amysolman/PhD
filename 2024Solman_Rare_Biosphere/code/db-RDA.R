# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Get environmental data and reduce count data to match samples
# 4. Perform db-RDA
# 5. Perform forward selection on variables
# 6. Report term % Explained, pseudo-F and P values for all terms tested
# 7. Create plot of final model


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
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(dplyr)

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied-var-trans.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied-var-trans.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
# pro.arc <- subset_samples(pro, Pole=="Arctic")
# pro.arc = filter_taxa(pro.arc, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)
# euk.arc <- subset_samples(euk, Pole=="Arctic")
# euk.arc = filter_taxa(euk.arc, function(x) sum(x) >= 1, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")

# Method

# Prior to analysis environmental variables were z-score transformed and variables with high multicollinearity removed. Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarities was performed to explore variation between community compositions among sites. Correlation between community composition and environmental variables was investigated using mantel tests on Bray-Curtis dissimilarity matrices and Euclidean distance matrix of environmental variables and geographic coordinates (Standardized Euclidean Distance). Partial mantel tests were performed, controlling for the effect of geographic coordinates and environmental distances. Bray-Curtis dissimilarities and Euclidean distances were calculated using the vegdist function of the vegan package. As different sets of environmental variables were available for different samples two analyses were performed. Canonical correspondence analysis (CCA) was carried out to explore the relationship between community structure and explanatory variables. **Distance-based redundancy analysis (db-RDA) was performed to explore linear relationships between environmental/ geographic variables and community composition. db-RDA was performed by carrying out PCoA on the Bray-Curtis dissimilarity matrix. The eigenvectors obtained through PCoA were then used, along with environmental/geographic variables, to carry out redundancy analysis. Forward variable selection was carried out to reduce models to their significant terms. db-RDA was carried out using the capscale function in the vegan package.** Variation partition analysis (VPA) was used to partition the explanatory power of environmental and physical factors. 
# 
# ## Distance-based Redundancy Analysis Explained
# 
# Distance-based redundancy analysis (db-RDA) is a method to conduct redundancy analysis (RDA) using dissimilarity metrics that are appropriate for microbiome count data. A dissimilarity matrix appropriate to the response (count) data is calculated and used as input in Principal Coordinates Analysis (PCoA). The most informative PCoA axis for each sample are then used as the response variable in place of count data dissimilarity matrices. These axis represent the community dissimilarities in **Euclidean space**. PCA and RDA are based on Euclidean distance, which is inappropriate for microbiome data as they are senstivite to the double-zero problem. Essentially, db-RDA allows ecologically relevant dissimilarity metrics (like Bray-Curtis) to be used and converted to Euclidean space, rather than just using Euclidean space directly. Partial RDA can also be performed to control for the effect of another set of explanatory variables (e.g. geographic coordinates).
# 
# db-RDA assumes that there is a linear relationship between explanatory and response variables. Usually in ecology responses to predicotr variables are unimodal instead of linear. When investigating unimodal relationships we should use CCA. In the case of unimodal responses to variables, distance-based CCA is more appropriate.
# 
# ### Choosing the Right Dissimilarity Metric
# 
# When choosing the right dissimilarity metric for you analysis you can:
#   A) Choose the one that is best for addressing the questions you are interested in (e.g. Phylogenetic turnover (UniFrac) or structural (Bray-Curtis/Sorenson/Jaccard) turnover? Abundance (Bray-Curtis/Weighted Unifrac) or presence/absence (Sorenson/Unweighted Unifrac)?)
# B) By looking at the rank correlations between each dissimilarity indices and gradient separation (higher values are better!). A good dissiilarity metric should have a high rank-order similarity with gradient separation. Metrics automatically included in the rankindex() function are:
#   Euclidean
# Mahhattan
# Gower
# Bray-Curtis
# Kulczynski
# 
# ### A note on ANOVA output
# 
# ####What is ANOVA?
# 
# When we run ANOVA on our model (CCA/RDA etc) we are performing an ANOVA-like permutation test to assess the significance of the constraints (explanatory variables). We can use this to test the significance of the model overall, or we can use it to assess the significance of the individual terms in the model. A permutation test is a non-parametric test whereby the REAL relationship (regression) betweeen our explanatory variables and our response variable is compared to the relationship in the values of the explanatory variable were shuffled randomly. This random shuffle gives us our test statistic. We do this many times (permutations) until we get our test statistic distribution. We can then see where our REAL test statistic falls within this distribution. This gives us the MAGICAL p-value! The probability of observing this relationship (regression) by chance or under out null-hypothesis scenario. 
# 
# Analysis of Variance (ANOVA) is a statistical method that separates observed variance (in the responses variable) into different groups. Essentially it consideres two sources of variance, between-group variance and within-group variance. The model calculates the mean of each group and compares this to the overall mean of the data. The within-group variance is the variation of each datapoint (observation) from that groups mean. Between-group and within-group variance is quantified by calculating a sum of squares (this is the distance of each point to the mean). The ratio of between group SS and within group SS gives us the F-statistic. The F-statistic is combined with degrees of freedom to calculate the p-statistic. So LARGE F-statistics suggest there is greater between group variation than within group variation. This then suggests that there are SIGNIFICANT differences between groups = a small p-value. So when we do an ANOVA-like permutation test we are comparing the TRUE data to RANDOMLY assembled data: we calculate variation (sum of squares) within REAL and RANDOM (PERMUTATION) results and between REAL and RANDOM results. We find the ratio of SS BETWEEN and WITHIN results (high ratio = more difference between than within) and use this with degrees of freedom (number of explanatory variables) to calculate the likelihood that these differences between REAL and RANDOM results would be observed by chance.  



# Results

# phylo = pro.ant.abun
# which_df = 1

dbRDA_func <- function(phylo, which_df){
  
  set.seed(666) # so we get the same results for each ANOVA-like permutation 
  
  #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
  phylo.rel = transform_sample_counts(phylo, function(x) x / sum(x) )
  
  #get community data
  comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
  
  # 3. Get environmental data and reduce count data to match samples
  
  #get environmental data
  vars = data.frame(sample_data(phylo.rel))
  if (which_df == 1){
    keep = c("Ca", "Mg", "K", "Na", "SO4", "Cl", "NO3", "EC", "pH", "Area", "TotalDepth")
  } else if (which_df == 2){
    keep = c("DistanceToSea", "Altitude", "HCO3")
  }
  var.trim <- vars[,names(vars) %in% keep]
  chem <- var.trim[complete.cases(var.trim),]
  
  #make sure all the same samples are used 
  comm <- comm[rownames(comm) %in% rownames(chem),]
  
  # 4. Perform db-RDA
  
  #perform db-RDA
  dbRDA_env = capscale(comm ~ ., chem, dist="bray", add=TRUE)
  
  #check for multicollinearity
  vif.out = data.frame(vif.cca(dbRDA_env))
  
  #keep only variables with VIF of 20 or less
  chem.keep = rownames(vif.out)[vif.out$vif.cca.dbRDA_env. <= 20]
  
  #re-run model
  chem = data.frame(chem[,names(chem) %in% chem.keep])
  names(chem) = chem.keep
  dbRDA_env = capscale(comm ~ ., chem, dist="bray", add=TRUE)
  
  #simplify the model by selecting only significant variables
  fwd.sel = ordistep(capscale(comm ~ 1, chem, dist="bray", add=TRUE), # lower model limit (simple!)
                     scope = formula(dbRDA_env))
  
  #Look at the new model with forward selected variables
  fwd.sel$call
  
  # Write our new model
  my.formula = fwd.sel$terms
  dbRDA_env_final = capscale(noquote(my.formula), chem, dist="bray", add=TRUE)
  
  return(list(dbRDA_env_final, chem.keep))
  
}



# phylo = euk.ant.rare
# subcom = "Rare"

get_sig_terms <- function(phylo, subcom){
  
  set.seed(666)
  
  out1 <- dbRDA_func(phylo, 1)
  out2 <- dbRDA_func(phylo, 2)
  
  #%Explained of abundant prokaryotes
  sig1 = anova(out1[[1]], by="terms")
  sig2 = anova(out2[[1]], by="terms")
  
  t1 = rownames(sig1)[1:nrow(sig1)-1]
  t2 = rownames(sig2)[1:nrow(sig2)-1]
  
  #save my values
  perc1 = vector()
  p_F1 = vector()
  p_val1 = vector()
  
  for (i in 1:length(t1)){
    #percentage
    perc1 = c(perc1, round(sig1$SumOfSqs[[i]]/sum(sig1$SumOfSqs), 4))
    #pseudo F
    p_F1 = c(p_F1, sig1$F[[i]])
    #p
    p_val1 = c(p_val1,sig1$`Pr(>F)`[[i]])
  }
  
  #save my values
  perc2 = vector()
  p_F2 = vector()
  p_val2 = vector()
  
  for (j in 1:length(t2)){
    #percentage
    perc2 = c(perc2, round(sig2$SumOfSqs[[j]]/sum(sig2$SumOfSqs), 4))
    #pseudo F
    p_F2 = c(p_F2, sig2$F[[j]])
    #p
    p_val2 = c(p_val2,sig2$`Pr(>F)`[[j]])
  }
  
  df = data.frame(Analysis = c(rep(1, length(t1)), rep(2, length(t2))), Variables = c(t1, t2), "% Explained" = c(perc1*100, perc2*100), "pseudo-F" = c(round(p_F1, 2), round(p_F2, 2)), "p" = c(round(p_val1, 4), round(p_val2, 4)))
  names(df) = c("Analysis", "Variables", paste(subcom, "% Explained"), paste(subcom, "pseudo-F"), paste(subcom, "p"))
  
  #add all variables that were tested (i.e. those with VIF <= 20)
  vars.tested = c(out1[[2]], out2[[2]]) #all the variables tested
  x1 = setdiff(out1[[2]], df$Variables) #variables tested that were removed from the model selection process 
  x2 = setdiff(out2[[2]], df$Variables) 
  df2 = data.frame(Analysis=c(rep(1, length(x1)), rep(2, length(x2))), Variables = c(x1, x2), "% Explained" = 0, "pseudo-F" = 0, "p"=NA)
  names(df2) = c("Analysis", "Variables", paste(subcom, "% Explained"), paste(subcom, "pseudo-F"), paste(subcom, "p"))
  final.df = rbind(df, df2)
  
  return(final.df)
  
}


get_multi_sig_terms <- function(phylo1, phylo2, phylo3){
  
  set.seed(666)
  
  ps.out1 = get_sig_terms(phylo1, "Abun")
  ps.out2 = get_sig_terms(phylo2, "Int")
  ps.out3 = get_sig_terms(phylo3, "Rare")
  
  x = full_join(ps.out1, ps.out2,
                by=NULL)
  
  y = full_join(x, ps.out3,
                by=NULL)
  
  #remove any rows that say Model
  y2 = y[y$Variables != "Model",]
  
  #add other variables that were tested by removed from analysis
  #list of all the variables tested
  vars_test1 = c("Ca", "Mg", "K", "Na", "SO4", "Cl", "NO3", "EC", "pH", "Area", "TotalDepth")
  vars_test2 = c("DistanceToSea", "Altitude", "HCO3")
  
  vars_test = c(vars_test1, vars_test2)
  
  z = data.frame(matrix(ncol = ncol(y2), nrow=length(setdiff(vars_test, y2$Variables))) )
  names(z) = names(y2)
  z$Analysis = c(rep(1, length(setdiff(vars_test1, y2$Variables))), rep(2, length(setdiff(vars_test2, y2$Variables))))
  z$Variables = setdiff(vars_test, y2$Variables)
  
  final.df = rbind(y2, z)
  
  return(final.df)
  
}


#Get results
pro.res = get_multi_sig_terms(pro.ant.abun, pro.ant.int, pro.ant.rare)
euk.res = get_multi_sig_terms(euk.ant.abun, euk.ant.int, euk.ant.rare)

export_res = rbind(pro.res, euk.res)
export_res = cbind(Group=c(rep("Prokaryote", nrow(pro.res)), rep("Eukaryote", nrow(euk.res))), export_res)

write.csv(export_res, "../results/dbRDA-results.csv")

