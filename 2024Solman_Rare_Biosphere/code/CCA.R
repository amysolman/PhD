# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Get environmental data and reduce count data to match samples
# 4. Perform CCA
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
library(car)

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

## Method

#A Canonical Correspondence Analysis (CCA) was implemented to assess the impact of environmental and physical variables on Antarctic subcommunity structures. Two sets of environmental variables were used as not all environmental data was available for all variables. A step-wise model was used to select the environmental variables that best explain community structure.

## CCA Explained

#CCA is a form of constrained ordination that associates two or more quantitative data sets (where simple/unconstrained ordination analyses a single data matrix). CCA is the constrained form of correspondence analysis (CA). CA is a multivariate ordination analysis which takes variables from a single quantitative data set and finds new variables or axis (i.e. orthogonal components) that explain the structure of the data set. Where Principal Component Analysis (PCA) assumes that the data are linearly related to each other and to environmental gradients (more species at higher environmental variable levels). CA expects data sets to have a unimodal distribution which is what is usually observed in species data (e.g. more species at intermediate environmental variable level, less species at high and low levels). CCA is a CA that looks for associations among two sets of variables (two matrices of data). One matrix will have species observations, while the other will have environmental variables. The CCA is appropriate if the relationship between observed species and envirnmental gradient is unimodal. If the relationship is linear a Canonical Correlation Analysis (CANCOR) or Redundancy Analysis (RDA) should be used. Confusingly CANCOR is also known as CCA. RDA is a multivariate analogue of simple linear regression (i.e. what is the linear relationship between community structure and an environmental variable).

#Partial CCA controls for the effect of a third (conditioning) matrix. The effect of the conditioning matrix is removed from the observation (species count) matrix. The resulting reduced matrix is then used in CCA with the variable matrix of interest (e.g. impacts of environmental data when controlling for spatial effects).

## How do I know if I should use CA/CCA or PCA/RDA?

#PCA/RDA should be used with variables that have a roughly normal distribution and with few zeros (i.e. they have low turnover - similar communities - or a short gradient). CA/CCA should be used for variables with longer gradients (i.e. variables are present in some samples and not others). We can estimate gradient length using the decorana function. When Axis lengths is greater than 3 CA/CCA is appropriate. When the Axis length < 3 PCA/RDA should be used. 

# Results


phylo = euk.ant.abun
which_df = 1

CCA_func <- function(phylo, which_df){
  
  set.seed(666)
  
  #normalise the data by transforming into relative abundance because our subcommunity sizes will be different
  phylo.rel = transform_sample_counts(phylo, function(x) x / sum(x) )
  
  #get community data
  comm <- data.frame(t(otu_table(phylo.rel)), check.names=FALSE)
  
  # 3. Get environmental data and reduce count data to match samples
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
  
  # 4. Perform CCA
  
  #perform CCA
  ccamodel <- cca(comm ~., chem)
  
  #check for multicollinearity
  vif.out = data.frame(vif.cca(ccamodel))
  
  #keep only variables with VIF of 20 or less
  chem.keep = rownames(vif.out)[vif.out$vif.cca.ccamodel. <= 20]
  
  #re-run model
  chem = data.frame(chem[,names(chem) %in% chem.keep])
  names(chem) = chem.keep
  ccamodel <- cca(comm ~., chem)
  
  ccamodel
  
  #simplify the model by selecting only significant variables
  fwd.sel = ordistep(cca(comm ~ 1, chem), # lower model limit (simple!)
                     scope = formula(ccamodel))
  
  #Look at the new model with forward selected variables
  fwd.sel$call
  
  # Write our new model
  my.formula = fwd.sel$terms
  cca_final = cca(noquote(my.formula), chem)
  
  # plot(fwd.sel)
  # fwd.sel
  # summary(fwd.sel)
  
  return(list(cca_final, chem.keep))
  
}

# phylo = euk.ant.rare
# subcom = "Rare"

#phylo = phylo1

get_sig_terms <- function(phylo, subcom){
  
  set.seed(666)
  
  out1 <- CCA_func(phylo, 1)
  out2 <- CCA_func(phylo, 2)
  
  #%Explained of prokaryotes
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
    perc1 = c(perc1, round(sig1$ChiSquare[[i]]/sum(sig1$ChiSquare), 4))
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
    perc2 = c(perc2, round(sig2$ChiSquare[[j]]/sum(sig2$ChiSquare), 4))
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
  
  # df = data.frame(Variables = t1, "% Explained" = perc1*100, "pseudo-F" = c(round(p_F1, 2)), "p" = c(round(p_val1, 4)))
  # names(df) = c("Variables", paste(subcom, "% Explained"), paste(subcom, "pseudo-F"), paste(subcom, "p"))
  
  return(final.df)
  
}


# phylo1 = pro.ant.abun
# phylo2 = pro.ant.int
# phylo3 = pro.ant.rare

get_multi_sig_terms <- function(phylo1, phylo2, phylo3){
  
  set.seed(666)
  
  ps.out1 = get_sig_terms(phylo1, "Abun")
  ps.out2 = get_sig_terms(phylo2, "Int")
  ps.out3 = get_sig_terms(phylo3, "Rare")
  
  #do we have the same variables being tested for each group
  setdiff(ps.out1$Variables, ps.out2$Variables)
  
  x = full_join(ps.out1, ps.out2,
                by=NULL)
  
  y = full_join(x, ps.out3,
                by=NULL)
  
  #remove any rows that say Model
  y2 = y[y$Variables != "Model",]
  
  #add other variables that were tested but removed from analysis
  #list of all the variables tested
  # variables = meta_df(phylo1)
  # vars_test = c(names(variables[[1]]), names(variables[[2]]))
  # vars_test = c("Ca", "Mg", "K", "Na", "SO4", "Cl", "NO3", "EC", "pH", "Area", "TotalDepth", "DistanceToSea", "Altitude", "HCO3")
  
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

write.csv(export_res, "../results/CCA-results.csv")
