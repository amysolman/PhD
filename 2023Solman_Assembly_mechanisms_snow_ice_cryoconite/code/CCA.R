

#Method

# Before CCA, db-RDA, VPA and Mantel analyses, all environmental variables were z-score transformed to improve normality and homoscedasticity. 
# 
# Canonical Correspondence Analysis (CCA) was implemented to assess the impact of environmental and physical variables on community structures. CCA is a form of constrained ordination that associated two or more quantitative datasets, in this case community profiles and environmental data. CCA is the constrained form of correspondence analysis (CA). CA expects data to have unimodal distributions. As such, CCA is appropriate if the relationship between species observations and environmental gradient is unimodal (increases at moderate values i.e. neutral pH). If the relationship is linear than redundancy analysis (a linear method) is appropriate. RDA should be used for data with a roughly normal distribution and few zeros (i.e. a community with low turnover/ high similarity - a short gradient). CCA should be used for variables with longer gradients (i.e. high turnover, many zeros) which is common for microbiome data. Gradient length can be estimated using the decorana function. When Axis lengths is greater than 3 CA/CCA is appropriate. When the Axis length < 3 PCA/RDA should be used.
# 
# Method: Proportionally transform community data. Perform CCA with community data and z-transformed environmental data using cca function from vegan package. Multicollinearity (linear dependancies) of variables (constraints) was checked within the model object using vif.cca function of vegan package. Constraints with VIF of 20 or less are retained and the model was re-run using this reduced list of variables. The resulting model was then simplified using forward stepwise model selection using the ordistep function of the vegan package. An ANOVA-like non-parametric permutation test was run on our final model to assess the significance of the individual terms compared to randomly assembled variable values to get our test statistic distribution. This gives us our F-value (ratio of between group sum of squares and within group sum of squares - high F statistic is good!) and p-value or the chance of observing these results by chance. The forward selected variables (significant terms) were reported, along with the percentage variation explained by each term, the pseudo-F value and P-value.
# 
# Distance-based Redundancy Analysis (db-RDA). db-RDA is a method for conducting redundancy analysis (an extension of simple linear regression - constrained ordination that assumes linear relationships), that uses dissimilarity matrices that are appropriate for microbiome data. The reason for this is that RDA is based on Euclidean space, which is inappropriate for microbiome data. The reason Euclidean space is inappropriate is because of the double-zero problem. That is, where two sites are missing the same organisms they are considered as similar as if they shared that feature/ASV in common. This is clearly biologically untrue. The absence of a bacterial taxa does not make two sites more similar. In db-RDA we calculate Bray-Curtis dissimilarity matrix for our communities. We then use this as input in Principal Coordinates Analysis (PCoA). The most informative PCoA axis are used as the response variable for RDA. These axis represent conversion of ecologically relevant beta diversity dissimilarity matrices to Euclidean space. 
# 
# Method: Proportionally transform community data. Perform db-RDA with bray-curtis dissimilarities of community data and z-transformed environmental data using capscale function from vegan package. Multicollinearity (linear dependancies) of variables (constraints) was checked within the model object using vif.cca function of vegan package. Constraints with VIF of 20 or less are retained and the model was re-run using this reduced list of variables. The resulting model was then simplified using forward stepwise model selection using the ordiR2step (was ordistep, ordiR2step used instead as it is designed for use with rda and capscale functions) function of the vegan package. An ANOVA-like non-parametric permutation test was run on our final model to assess the significance of the individual terms compared to randomly assembled variable values to get our test statistic distribution. This gives us our F-value (ratio of between group sum of squares and within group sum of squares - high F statistic is good!) and p-value or the chance of observing these results by chance. The forward selected variables (significant terms) were reported, along with the percentage variation explained by each term, the pseudo-F value and P-value.
# 
# Variation Partitioning Analysis (VPA). VPA partitioned the explanatory power of environmental and geographic matrices in relation to a community dissimilarity matrix. Multiple partial RDAs are run to determine the linear effect of each explanatory matrix. The shared partition of explanatory power highlighted by VPA is not an interaction term, but highlights multicollinearity in the model. 
# 
# Method: Proportionally transform community data. db-RDA was performed using community data and environmental variables. Multicollinearity (linear dependancies) of variables (constraints) was checked within the model object using vif.cca function of vegan package. Constraints with VIF of 20 or less are retained. Latitude and longitude coordinates were retained. A Euclidean distance matrix using geographic coordinates was generated using the dist function of the stats package. Principal coordinates of neighbourhood matrix were computed by the principal coordinate analysis of the distance matrix using the pcnm function of the vegan package. This method is used to transform spatial distances to rectangular data that is suitable for constrained ordination. The scores of the new PCNM variables was extracted using the scores function in vegan. Forward selection was carried out to select significant geographic PCNM axis using the ordiR2step (was ordistep, ordiR2step used instead as it is designed for use with rda and capscale functions) function of the vegan package. The variation of our community matrix was then partitioned by our reduced environmental variable matrix and reduced geographic PCNM matrix using the varpart function of the vegan package. Visualisation of VPA analysis was carried out using base R.
# 
# Mantel and Partial Mantel Tests. Proportionally transform community data. RDA (previously was CCA) was performed using community data and environmental variables. Multicollinearity (linear dependancies) of variables (constraints) was checked within the model object using vif.cca function of vegan package. Constraints with VIF of 20 or less were retained in order to only use that incorporated into VPA. Latitude and longitude coordinates were retained. A Bray-Curtis dissimilarity matrix was generated for the community data. A Euclidean distance matrix using geographic coordinates, and another using environmental data, was generated using the vegdist function of the vegan package. Mantel tests using spearman’s rank correlations was carried out between the community dissimilarity matrix and each of the environmental and geographic matrices with 999 permutations using the mantel function of the vegan package. Partial Mantel tests using spearman’s rank correlations were carried out between the community dissimilarity matrix and each of the environmental and geographic matrices while controlling for the other, with 999 permutations using the mantel function of the vegan package. Mantel statistics (rho) and p-values were reported.


# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

# source("00-solman-functions.R")
library(phyloseq)
# library(RColorBrewer)
# library(tidyr) #trans data between wide and long format
# library(ggplot2)
library(vegan) #for vegdist
# library(BiodiversityR) #for rank abundance curves
# library(dplyr) #for %>% function
# library(cowplot) #plot_grid
# library(reshape2) #for function melt
# library(funrar) #for make relative
# library(stringr) #mutate function, to split columns
# library(gridExtra) #for exporting as pdf
# library(scales) #for percentages
# library(microbiome) #for summarize_phyloseq
# library(picante) #for faith's pd and comdistnt functions
# library(ggpmisc) #for adding lm stats to the plot
# library(broom) #for adding stats to the plot
# library(fossil) #for earth.dist function
# library(Hmisc) #binconf function
# library(minpack.lm) #for non-linear model fitting
# library(tibble) #function rownames_to_column
# library(stats4) #mle function
# library(ecodist) #for distance() function
# library(parallel)
# library(svglite)
library(Hmisc) #binconf and rcorr functions
library(corrplot) #for plotting our spearmans correlations
library(ggvegan) #for cca plots
library(cowplot)
library(dplyr)
# install.packages("ggpp")
library(ggpp)
library(patchwork)
library(egg)
library(ggpubr)
library(tibble)

#prokaryotes
ps.pro <- readRDS("../results-proportional/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results-proportional/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results-proportional/18S-mm-only-ps-norm.rds") 


#Pairwise scatter plots of environmental variables to see which are highly correlated

#correlation test function

# ps = ps.pro
# hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID
# habitat = "Cryoconite"
# group = "prokaryotes"

# correlated_vars <- function(ps, hab, habitat, group){
#   
#   #subset data by habitat
#   ps.hab <- prune_samples(hab, ps)
#   
#   #get metadata
#   meta = data.frame(sample_data(ps))
#   meta = meta[,c(9:ncol(meta))]
#   
#   #remove columns with all NA values
#   #meta.trim <- meta[ , colSums(is.na(meta)) < nrow(meta)]
#   
#   #remove columns with all zero or NA values
#   meta.no0.na = Filter(function(x) !all(is.na(x)|x == 0), meta)
#   
#   #compute correlation matrix
#   res1 <- rcorr(as.matrix(meta.no0.na), type=c("spearman"))
#   
#   #format the correlation matrix
#   #flattenCorrMatrix function taken from http://sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
#   flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#     )
#   }
#   
#   #results dataframe
#   res.flat <- flattenCorrMatrix(res1$r, res1$P)
#   
#   #change NA values to 1 for P value and 0 for r value 
#   res1$P[is.na(res1$P)] <- 1
#   res1$r[is.na(res1$r)] <- 0
#   
#   corrplot(res1$r, type="upper", order="hclust", 
#          p.mat = res1$P, sig.level = 0.05, insig = "blank")
#   print(p)
#   
#   #plot
#   pdf(paste0("../results-proportional/", habitat, "-", group, "-variable-correlations.pdf"), width=13, height=13)
#   corrplot(res1$r, type="upper", order="hclust", 
#          p.mat = res1$P, sig.level = 0.05, insig = "blank")
#    dev.off()
#    
#   #which correlations are significant?
#   res.rm <- res.flat[which(res.flat$p < 0.05 & abs(res.flat$cor) > 0.9),]
#   
#   #vector of variables we would preferentially remove
#   var.to.rm <- c("TN", "Br", "Li", "Al", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn", "As", "Se", "Y",
#                "Rb", "Sr", "Mo", "Ag", "Cd", "Ba", "Lu", "Pb", "Zr", "Nb", "Ru", "Sn", "Cs", "Hf", "Ta",
#                "Re", "U", "Sc", "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
#                "Yb", "Sb")
#   #which of our unwanted variables are found in our dataframe?
#   x = intersect(var.to.rm, unique(c(res.rm$row, res.rm$column)))
#   # x.col = intersect(var.to.rm, unique(res.rm$column))
#   # 
#   # #remove those unwanted variables from the first column
#   # x2 = res.rm[!(res.rm$row %in% x.row),]
#   # x3 = x2[!(x2$column %in% x.col),]
#   # 
#   # #variables that need to be removed
#   # var.rm.me = unique(c(x.row,x.col))
#   
#   write.csv(x, paste0("../results/", habitat, "-", group, "-correlated-variables.csv"), row.names = FALSE)
#   
#   return(res.rm)
#   
# }


# pro.snow <- correlated_vars(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, habitat="Snow", group="prokaryotes")
# 
# pro.sp <- correlated_vars(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, habitat="SpringIce", group="prokaryotes")
# 
# pro.sum <- correlated_vars(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, habitat="SummerIce", group="prokaryotes")
# 
# pro.cry <- correlated_vars(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite", group="prokaryotes")
# 
# euk.snow <- correlated_vars(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, habitat="Snow", group="eukaryotes")
# 
# euk.sp <- correlated_vars(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, habitat="SpringIce", group="eukaryotes")
# 
# euk.sum <- correlated_vars(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, habitat="SummerIce", group="eukaryotes")
# 
# euk.cry <- correlated_vars(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite", group="eukaryotes")
# 
# mm.snow <- correlated_vars(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, habitat="Snow", group="micrometazoans")
# 
# mm.sp <- correlated_vars(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, habitat="SpringIce", group="micrometazoans")
# 
# mm.sum <- correlated_vars(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, habitat="SummerIce", group="micrometazoans")
# 
# mm.cry <- correlated_vars(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite", group="micrometazoans")


#CCA
#Methods for dealing with high VIF values
#https://quantpalaeo.wordpress.com/2014/04/14/variance-inflation-factors-and-ordination-model-selection/
  
  #Tutorial for CCA
#   https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html
# https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html


# ps = ps.pro
# hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID
# group = "Prokaryote"
# habitat = "Spring Ice"
# col = "#a4c64d"

# min(taxa_sums(ps.pro))

cca_func <- function(ps, hab, group, habitat, col){
  
  ##############################PREP THE DATA##############################
  
  #set the seed to get repeatable results were using randomized function
  set.seed(666)
  
  #subset the phyloseq object by habitat
  h <- prune_samples(hab, ps)
  
  #remove taxa with zero counts after subsetting the data
  s <- prune_taxa(taxa_sums(h) > 0, h)
  
  #extract the community profile
  comm = data.frame(t(otu_table(s)), check.names = FALSE)
  min(colSums(comm)) == 0 #should be FALSE
  
  #get our environmental data from the phyloseq object, minus the columns we're not interested in
  #df = data.frame(sample_data(s))
  vars = data.frame(sample_data(s))[,-c(1:8, 11:13)]
  
  #remove rare earth elements/lanthanides
  rm <- c("Ln", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
          "Dy", "Ho", "Er", "Tm", "Yb", "Lu")
  vars = vars[,! names(vars) %in% rm]
  
  #change the name of conductivity
  colnames(vars)[colnames(vars) == 'Conductivity_muS'] <- 'Conductivity'
  
  #remove columns with less than 50% values != 0 or NA
  c <- vars
  c[c > 0 | c < 0] <- 1 #give values greater than 0 or less than 0 the value of 1
  c[is.na(c)] <- 0 #make NA 0
  var.trim <- vars[,colSums(c) >= round(nrow(vars)*0.50)] #only keep those with at least 50% values != 0 or NA
  #only keep complete cases
  var.comp <- var.trim[complete.cases(var.trim),]
  # var.comp <- vars[complete.cases(vars),]
  
  #reduce the number of variables for micrometazoans because we have less data points
  if(group == "Micrometazoan"){
    var.comp = var.comp[,names(var.comp) %in% c("pH", "Conductivity", "TC", "Br", "Fe", "Na", "Mg", "K", "Ca")]
  }
  
  #remove columns with all zero values
  var.comp <- var.comp[,colSums(var.comp) > 0]
  
  #z-score transform the environmental data
  var.z = data.frame(scale(var.comp))
  
  #make sure our samples are the same for comm and env
  comm = comm[rownames(comm) %in% rownames(var.z),]
  
  ###############################RUN THE MODEL###############################
  
  ########################ALL VARIABLES WITHOUT REMOVING COLLINEAR VARIABLES###############
  #perform CCA
  ccamodel1 <- cca(comm ~., var.z) #model with all explanatory variables
  
  ccamodel1 #look at the proportion of variation constrained by the model
  
  summary(ccamodel1)
  
  RsquareAdj(ccamodel1) #total variance explained by the model
  
  #check the significance of the whole model - GLOBAL SIGNIFICANCE
  anova.cca(ccamodel1, permutations=999) 
  
  #check for significance of the axis
  anova.cca(ccamodel1, by="axis")
  
  #check the significance of the terms in the model
  anova.cca(ccamodel1, by="terms") 
  
  #######################REMOVING COLLINEAR VARIABLES##########################
  
  #remove variables with high multicollinearity
  vif.out = data.frame(vif.cca(ccamodel1)) #using z-score transformed data has no effect on vif
  vif.out[is.na(vif.out)] <- 0
  
  # if one variable has VIF > 10 then remove the variable with highest VIF and re-run the ordination
  repeat{
    if (max(vif.out) > 10){
      vif.out$variable = rownames(vif.out)
      names(vif.out) = c("VIF", "Variable")
      vars.keep = vif.out[which(vif.out$VIF < max(vif.out$VIF)),]$Variable
      var.z.trim = data.frame(var.z[,colnames(var.z) %in% vars.keep])
      
      #check for multicollinearity
      vif.out = data.frame(vif.cca(cca(comm ~., var.z.trim))) #using z-score transformed data has no effect on vif
      vif.out[is.na(vif.out)] <- 0
      
    }
    
    if (max(vif.out) < 10){
      break
    }
    
  }  
  
  ########################MOLDEL WITH COLLINEAR VARIABLES REMOVED###############
  
  #re-run the model with collinear variables removed
  ccamodel2 <- cca(comm ~., var.z.trim) #model with all explanatory variables
  ccamodel2 #look at the proportion of variation constrained by the model
  summary(ccamodel2)
  ano.res <- anova.cca(ccamodel2)#check the significance of the whole model
  pval <- ano.res$`Pr(>F)`[[1]]
  print(pval)
  anova.cca(ccamodel2) #check the significance of the whole model
  anova.cca(ccamodel2, by="axis") #check the significance of the first two axis
  anova.cca(ccamodel2, by="terms") #check the significance of the terms in the model
  
  #global model r2
  global_r2 <- RsquareAdj(ccamodel2)$adj.r.squared
  
  # # Forward selection of variables:
  # fwd.sel <- ordiR2step(cca(comm ~1, var.z.trim), # lower model limit (simple!)
  #              scope = formula(ccamodel2), # upper model limit (the "full" model)
  #              direction = "forward",
  #              R2scope = TRUE, # can't surpass the "full" model's R2
  #              pstep = 1000,
  #              trace = FALSE) # change to TRUE to see the selection process!
  # fwd.sel$call
  # 
  # #global model r2
  # global_r2 <- RsquareAdj(ccamodel2)$adj.r.squared
  # 
  # #look at the new model
  # fwd.sel #look at the proportion of variation constrained by the model
  # summary(fwd.sel)
  # anova.cca(fwd.sel) #check the significance of the whole model
  # anova.cca(fwd.sel, by="axis")
  # anova.cca(fwd.sel, by="terms") #check the significance of the terms in the model
  # 
  # #stepwise model selection using ordistep
  # ccamodel0 <- cca(comm ~1, var.z.trim) #model with intercept only
  # finalmodel <-  ordistep(ccamodel0, scope = formula(ccamodel2))
  # finalmodel
  # anova(finalmodel)
  # anova(finalmodel, by="terms")
  
  
  if (ncol(var.z.trim) > 0){
    
    sig = anova(cca(comm ~., var.z.trim), by="terms")
    
    #put these values into a dataframe
    df = data.frame(Variables = rownames(sig), "% Explained" = round(sig$ChiSquare/sum(sig$ChiSquare), 4)*100, "pseudo-F" = round(sig$F, 2), "p" = round(sig$`Pr(>F)`, 4))
    names(df) = c("Variables", "% Explained", "pseudo-F", "p")
    
    cca.final.df = df
    cca.final.df$Group = group
    cca.final.df$Habitat = habitat
    cca.final.df$Ordination = "CCA"
    
  } 
  
  if (any(cca.final.df$p < 0.05, na.rm=TRUE)){
    
    #get the model data
    fmod <- fortify(ccamodel2) 
    
    #get site data
    data = data.frame(sample_data(s))
    data$Location = paste(data$Location, "Foxfonna")
    
    #subset to include only those samples kept in the model
    data.sub = data[data$SampleID %in% subset(fmod, Score == "sites")$Label,]
    
    #get the environment vars CC1 and CC2 scores
    df_environ  <- as.data.frame(scores(ccamodel2, display = 'bp') )
    
    #subset env data to significant terms only
    df_environ_sig <- df_environ[rownames(df_environ) %in% cca.final.df[cca.final.df$p < 0.05,]$Variables,]
    
    #subset env data to non-significant terms only
    df_environ_no_sig <- df_environ[rownames(df_environ) %in% cca.final.df[cca.final.df$p >= 0.05,]$Variables,]
    
    #Get percentage of variance explained by first axis
    cca1_varex<-round(summary(ccamodel2)$cont$importance[2,1]*100,2) 
    
    #Get percentage of variance explained by second axis
    cca2_varex<-round(summary(ccamodel2)$cont$importance[2,2]*100,2) 
    
    #set a scaling factor
    scaling_factor <- 2
    
    if(any(names(fmod) == "CCA2")){
      
      #plot with ggplot
      p = ggplot(fmod, aes(x = CCA1, y = CCA2))+
        geom_point(data = cbind(subset(fmod, Score == "sites"), Location = data.sub$Location),
                   aes(fill = col, shape=Location, alpha=0.7), size = 4) +
        #scale_colour_brewer("Location", palette = "Set1") +
        scale_fill_manual(values=c(col))+
        guides(fill="none", alpha="none")+
        scale_shape_manual(values=c(21,24))+
        coord_fixed() +
        theme(legend.position = "top")+
        #add significant variables
        geom_segment(data=df_environ_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CCA2*scaling_factor), #Ending coordinate in CCA2 
                     color="firebrick1", #set color
                     linewidth=1,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        #add non-significant variables
        geom_segment(data=df_environ_no_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CCA2*scaling_factor), #Ending coordinate in CCA2 
                     #color="blue", #set color
                     linewidth=0.2,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        #Add environmental vars text significant
        geom_text(data=df_environ_sig, 
                  aes(x=CCA1*scaling_factor, 
                      y=CCA2*scaling_factor,
                      label=rownames(df_environ_sig),
                      hjust=0.5*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(CCA2))),#Add the text of each environmental var at the end of the arrow 
                  color="firebrick1",
                  size=4,
        )+
        #Add environmental vars text non-significant
        geom_text(data=df_environ_no_sig, 
                  aes(x=CCA1*scaling_factor, 
                      y=CCA2*scaling_factor,
                      label=rownames(df_environ_no_sig),
                      hjust=0.5*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(CCA2))),#Add the text of each environmental var at the end of the arrow 
                  #color="blue"
                  size=2,
        )+
        labs(x=paste0("CCA1 (",cca1_varex," %)"),
             y=paste0("CCA2 (",cca2_varex," %)"))+
        theme_bw()+
        xlim(-1.5, 3.5)+
        ylim(-2.5, 2.5)
      #ggtitle(habitat)
      # geom_text(x=min(fmod$CCA1), y=min(fmod$CCA2), label=paste("AdjR2", round(global_r2,2)))+
      #geom_text_npc(aes(npcx = "left", npcy = "top", label = paste("AdjR2", round(global_r2,2), "PVal", pval)))
      
      
      print(p)
      
    } else {
      
      
      #plot with ggplot
      p = ggplot(fmod, aes(x = CCA1, y = CA1))+
        geom_point(data = cbind(subset(fmod, Score == "sites"), Location = data.sub$Location),
                   aes(colour = Location), size = 3) +
        scale_colour_brewer("Location", palette = "Set1") +
        coord_fixed() +
        theme(legend.position = "top")+
        geom_segment(data=df_environ_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CA1*scaling_factor), #Ending coordinate in CA1 
                     #color="firebrick1", #set color
                     linewidth=1,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        geom_segment(data=df_environ_no_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CCA1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CA1*scaling_factor), #Ending coordinate in CA1 
                     #color="firebrick1", #set color
                     linewidth=0.2,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        #Add environmental vars text
        geom_text(data=df_environ, 
                  aes(x=CCA1*scaling_factor, 
                      y=CA1*scaling_factor,
                      label=rownames(df_environ),
                      hjust=0.5*(1-sign(CCA1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(CA1))),#Add the text of each environmental var at the end of the arrow 
                  #color="firebrick1"
                  size=4,
        )+
        labs(x=paste0("CCA1 (",cca1_varex," %)"),
             y=paste0("CA1"))+
        theme_bw()+
        xlim(min(fmod$CCA1),max(fmod$CCA1)+1)+
        ggtitle(habitat)#+
      #geom_text_npc(aes(npcx = "left", npcy = "top", label = paste("AdjR2", round(global_r2,2), "PVal", pval)))
      
      
      
    }
    
    #plot(p)
    
    
    return(list(cca.final.df, p, pval, global_r2))
    
  } else if (all(cca.final.df$p > 0.05, na.rm=TRUE)){
    
    return(cca.final.df)
    
  }
  
  
  
}

```

```{r, echo=FALSE}
#prokaryotes
# cca.pro.sn <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow", col="#fa9f99")
# cca.pro.sn[[1]] #table only
cca.pro.sp <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice", col="#a4c64d")
cca.pro.sp[[1]] #results
cca.pro.sp[[2]] #plot
cca.pro.sp[[3]] #pval
cca.pro.sp[[4]] #adjusted r
# cca.pro.sm <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice", col="#4dd2d6")
# cca.pro.sm[[1]] #table only
cca.pro.cr <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite", col="#d8a4ff")
# cca.pro.cr[[1]] #results
# cca.pro.cr[[2]] #plot
# cca.pro.cr[[3]] #pval
# cca.pro.cr[[4]] #adjusted r

#eukaryotes
# cca.euk.sn <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Microbial eukaryote", "Snow", col="#fa9f99")
# cca.euk.sn[[1]] #results
# cca.euk.sn[[2]] #plot
# cca.euk.sn[[3]] #pval
# cca.euk.sn[[4]] #adjusted r
# cca.euk.sp <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Microbial eukaryote", "Spring Ice", col="#a4c64d")
# cca.euk.sp[[1]] #results
# cca.euk.sp[[2]] #plot
# cca.euk.sp[[3]] #pval
# cca.euk.sp[[4]] #adjusted r
# cca.euk.sm <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Microbial eukaryote", "Summer Ice", col="#4dd2d6")
# cca.euk.sm[[1]] #table only
cca.euk.cr <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Microbial eukaryote", "Cryoconite", col="#d8a4ff")
# cca.euk.cr[[1]] #results
# cca.euk.cr[[2]] #plot
# cca.euk.cr[[3]] #pval
# cca.euk.cr[[4]] #adjusted r

#microfauna
# cca.mm.sn <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Microfauna", "Snow", col="#fa9f99")
# cca.mm.sn[[1]] #table only
# cca.mm.sp <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Microfauna", "Spring Ice", col="#a4c64d")
# cca.mm.sp[[1]] #table only
# cca.mm.sm <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Microfauna", "Summer Ice", col="#4dd2d6")
# cca.mm.sm[[1]] #results
# cca.mm.sm[[2]] #plot
# cca.mm.sm[[3]] #pval
# cca.mm.sm[[4]] #adjusted r
cca.mm.cr <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Microfauna", "Cryoconite", col="#d8a4ff")
# cca.mm.cr[[1]] #results
# cca.mm.cr[[2]] #plot
# cca.mm.cr[[3]] #pval
# cca.mm.cr[[4]] #adjusted r
```


```{r, echo=FALSE}
#table order of variables
geochem.ord = c("pH", "Conductivity", "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
                "Cl", "Br", "Fl", "Fe", "Na", "Mg", "K", "Ca", "Li", "Al", "Ti", "V", "Cr", "Mn", "Co",
                "Ni", "Cu", "Zn", "Y", "Rb", "Sr", "Mo", "Ag", "Ba", "Lu", "Pb", "Zr", "Cd", "As", "Se", "Nb", "Sn", "Re", "U", "La", "Residual")

#bind into one dataframe
names(cca.pro.sn)[2:4] = paste(names(cca.pro.sn)[2:4], "16S-Snow")
names(cca.pro.sp[[1]])[2:4] = paste(names(cca.pro.sp[[1]])[2:4], "16S-Spring")
names(cca.pro.sm)[2:4] = paste(names(cca.pro.sm)[2:4], "16S-Summer")
names(cca.pro.cr[[1]])[2:4] = paste(names(cca.pro.cr[[1]])[2:4], "16S-Cryo")

#use dplyr to join dfs
cca.pro.df = full_join(cca.pro.sn[,1:4], cca.pro.sp[[1]][,1:4], by = 'Variables')
cca.pro.df = full_join(cca.pro.df, cca.pro.sm[,1:4], by = 'Variables')
cca.pro.df = full_join(cca.pro.df, cca.pro.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.pro.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.pro.df$Variables)]
x = cca.pro.df[match(geo.keep, cca.pro.df$Variables),]

write.csv(x, "../results-proportional/cca-pro-df.csv")

#bind into one dataframe
names(cca.euk.sn[[1]])[2:4] = paste(names(cca.euk.sn[[1]])[2:4], "18S-Snow")
names(cca.euk.sp[[1]])[2:4] = paste(names(cca.euk.sp[[1]])[2:4], "18S-Spring")
names(cca.euk.sm)[2:4] = paste(names(cca.euk.sm)[2:4], "18S-Summer")
names(cca.euk.cr[[1]])[2:4] = paste(names(cca.euk.cr[[1]])[2:4], "18S-Cryo")

#use dplyr to join dfs
cca.euk.df = full_join(cca.euk.sn[[1]][,1:4], cca.euk.sp[[1]][,1:4], by = 'Variables')
cca.euk.df = full_join(cca.euk.df, cca.euk.sm[,1:4], by = 'Variables')
cca.euk.df = full_join(cca.euk.df, cca.euk.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.euk.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.euk.df$Variables)]
y = cca.euk.df[match(geo.keep, cca.euk.df$Variables),]

write.csv(y, "../results-proportional/cca-euk-df.csv")

#bind into one dataframe - micrometazoans
names(cca.mm.sn)[2:4] = paste(names(cca.mm.sn)[2:4], "MM-Snow")
names(cca.mm.sp)[2:4] = paste(names(cca.mm.sp)[2:4], "MM-Spring")
names(cca.mm.sm[[1]])[2:4] = paste(names(cca.mm.sm[[1]])[2:4], "MM-Summer")
names(cca.mm.cr[[1]])[2:4] = paste(names(cca.mm.cr[[1]])[2:4], "MM-Cryo")

#use dplyr to join dfs
cca.mm.df = full_join(cca.mm.sn[,1:4], cca.mm.sp[,1:4], by = 'Variables')
cca.mm.df = full_join(cca.mm.df, cca.mm.sm[[1]][,1:4], by = 'Variables')
cca.mm.df = full_join(cca.mm.df, cca.mm.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.mm.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.mm.df$Variables)]
z = cca.mm.df[match(geo.keep, cca.mm.df$Variables),]

write.csv(z, "../results-proportional/cca-mm-df.csv")
```

#########RDA

#Tutorial for RDA
https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html
https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html
```{r, include=FALSE}
# ps = ps.pro
# hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID
# group = "Prokaryotes"
# habitat = "Cryoconite"
# col = "#d8a4ff"

# min(taxa_sums(ps.pro))

dbRDA_func <- function(ps, hab, group, habitat, col){
  
  ##############################PREP THE DATA##############################
  
  #set the seed to get repeatable results were using randomized function
  set.seed(666)
  
  #subset the phyloseq object by habitat
  h <- prune_samples(hab, ps)
  
  #remove taxa with zero counts after subsetting the data
  s <- prune_taxa(taxa_sums(h) > 0, h)
  
  #extract the community profile
  comm = data.frame(t(otu_table(s)), check.names = FALSE)
  min(colSums(comm)) == 0 #should be FALSE
  
  #get our environmental data from the phyloseq object, minus the columns we're not interested in
  #df = data.frame(sample_data(s))
  vars = data.frame(sample_data(s))[,-c(1:8, 11:13)]
  
  #remove rare earth elements/lanthanides
  rm <- c("Ln", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
          "Dy", "Ho", "Er", "Tm", "Yb", "Lu")
  vars = vars[,! names(vars) %in% rm]
  
  #change the name of conductivity
  colnames(vars)[colnames(vars) == 'Conductivity_muS'] <- 'Conductivity'
  
  #remove columns with less than 50% values != 0 or NA
  c <- vars
  c[c > 0 | c < 0] <- 1 #give values greater than 0 or less than 0 the value of 1
  c[is.na(c)] <- 0 #make NA 0
  var.trim <- vars[,colSums(c) >= round(nrow(vars)*0.50)] #only keep those with at least 50% values != 0 or NA
  #only keep complete cases
  var.comp <- var.trim[complete.cases(var.trim),]
  # var.comp <- vars[complete.cases(vars),]
  
  #reduce the number of variables for micrometazoans because we have less data points
  if(group == "Micrometazoan"){
    var.comp = var.comp[,names(var.comp) %in% c("pH", "Conductivity", "TC", "Br", "Fe", "Na", "Mg", "K", "Ca")]
  }
  
  #remove columns with all zero values
  var.comp <- var.comp[,colSums(var.comp) > 0]
  
  #z-score transform the environmental data
  var.z = data.frame(scale(var.comp))
  
  #make sure our samples are the same for comm and env
  comm = comm[rownames(comm) %in% rownames(var.z),]
  
  ###############################RUN THE MODEL###############################
  
  ########################ALL VARIABLES WITHOUT REMOVING COLLINEAR VARIABLES###############
  
  #perform db-RDA
  rdamodel1 <- capscale(comm ~., var.z, dist="bray", add=TRUE) #model with all explanatory variables
  
  rdamodel1 #look at the proportion of variation constrained by the model
  
  RsquareAdj(rdamodel1) #total variance explained by the model
  
  #check the significance of the whole model - GLOBAL SIGNIFICANCE
  anova.cca(rdamodel1, permutations=999) 
  
  #check for significance of the axis
  anova.cca(rdamodel1, by="axis")
  
  #check the significance of the terms in the model
  anova.cca(rdamodel1, by="terms") 
  
  #######################REMOVING COLLINEAR VARIABLES##########################
  
  #remove variables with high multicollinearity
  vif.out = data.frame(vif.cca(rdamodel1)) #using z-score transformed data has no effect on vif
  vif.out[is.na(vif.out)] <- 0
  
  # if one variable has VIF > 10 then remove the variable with highest VIF and re-run the ordination
  repeat{
    if (max(vif.out) > 10){
      vif.out$variable = rownames(vif.out)
      names(vif.out) = c("VIF", "Variable")
      vars.keep = vif.out[which(vif.out$VIF < max(vif.out$VIF)),]$Variable
      var.z.trim = data.frame(var.z[,colnames(var.z) %in% vars.keep])
      
      #check for multicollinearity
      vif.out = data.frame(vif.cca(capscale(comm ~., var.z.trim))) #using z-score transformed data has no effect on vif
      vif.out[is.na(vif.out)] <- 0
      
    }
    
    if (max(vif.out) < 10){
      break
    }
    
  }  
  
  ########################MOLDEL WITH COLLINEAR VARIABLES REMOVED###############
  
  #re-run the model with collinear variables removed
  rdamodel2 <- capscale(comm ~., var.z.trim, dist="bray", add=TRUE) #model with all explanatory variables
  rdamodel2 #look at the proportion of variation constrained by the model
  summary(rdamodel2)
  ano.res <- anova.cca(rdamodel2)#check the significance of the whole model
  pval <- ano.res$`Pr(>F)`[[1]]
  anova.cca(rdamodel2, by="axis")
  anova.cca(rdamodel2, by="terms") #check the significance of the terms in the model
  
  #global model r2
  global_r2 <- RsquareAdj(rdamodel2)$adj.r.squared
  
  if (ncol(var.z.trim) > 0){
    
    sig = anova(capscale(comm ~., var.z.trim, dist="bray", add=TRUE), by="terms")
    
    #put these values into a dataframe
    df = data.frame(Variables = rownames(sig), "% Explained" = round(sig$SumOfSqs/sum(sig$SumOfSqs), 4)*100, "pseudo-F" = round(sig$F, 2), "p" = round(sig$`Pr(>F)`, 4))
    names(df) = c("Variables", "% Explained", "pseudo-F", "p")
    
    rda.final.df = df
    rda.final.df$Group = group
    rda.final.df$Habitat = habitat
    rda.final.df$Ordination = "RDA"
    
  } 
  
  if (any(rda.final.df$p < 0.05, na.rm=TRUE)){
    
    #get the model data
    fmod <- fortify(rdamodel2) 
    
    #get site data
    data = data.frame(sample_data(s))
    data$Location = paste(data$Location, "Foxfonna")
    
    #subset to include only those samples kept in the model
    data.sub = data[data$SampleID %in% subset(fmod, Score == "sites")$Label,]
    
    #get the environment vars CAP1 and CAP2 scores
    df_environ  <- as.data.frame(scores(rdamodel2, display = 'bp') )
    
    #subset env data to significant terms only
    df_environ_sig <- df_environ[rownames(df_environ) %in% rda.final.df[rda.final.df$p < 0.05,]$Variables,]
    
    #subset env data to non-significant terms only
    df_environ_no_sig <- df_environ[rownames(df_environ) %in% rda.final.df[rda.final.df$p >= 0.05,]$Variables,]
    
    #Get percentage of variance explained by first axis
    rda1_varex<-round(summary(rdamodel2)$cont$importance[2,1]*100,2) 
    
    #Get percentage of variance explained by second axis
    rda2_varex<-round(summary(rdamodel2)$cont$importance[2,2]*100,2) 
    
    #set a scaling factor
    scaling_factor <- 2
    
    if(any(names(fmod) == "CAP2")){
      
      #plot with ggplot
      p = ggplot(fmod, aes(x = CAP1, y = CAP2))+
        geom_point(data = cbind(subset(fmod, Score == "sites"), Location = data.sub$Location),
                   aes(fill = col, shape=Location, alpha=0.7), size = 4) +
        scale_fill_manual(values=c(col))+
        guides(fill="none", alpha="none")+
        scale_shape_manual(values=c(21,24))+
        coord_fixed() +
        theme(legend.position = "top")+
        geom_segment(data=df_environ_sig, #add significant arrows
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CAP1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CAP2*scaling_factor), #Ending coordinate in CCA2 
                     color="firebrick1", #set color
                     linewidth=1,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        geom_segment(data=df_environ_no_sig, #add non-significant arrows
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CAP1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=CAP2*scaling_factor), #Ending coordinate in CCA2 
                     #color="firebrick1", #set color
                     linewidth=0.2,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        #Add significant environmental vars text
        geom_text(data=df_environ_sig,
                  aes(x=CAP1*scaling_factor,
                      y=CAP2*scaling_factor,
                      label=rownames(df_environ_sig),
                      hjust=0.5*(1-sign(CAP1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(CAP2))),#Add the text of each environmental var at the end of the arrow
                  color="firebrick1",
                  size=4,
        )+
        #Add non-significant environmental vars text
        geom_text(data=df_environ_no_sig,
                  aes(x=CAP1*scaling_factor,
                      y=CAP2*scaling_factor,
                      label=rownames(df_environ_no_sig),
                      hjust=0.5*(1-sign(CAP1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(CAP2))),#Add the text of each environmental var at the end of the arrow
                  #color="firebrick1",
                  size=2,
        )+
        labs(x=paste0("RDA1 (",rda1_varex," %)"),
             y=paste0("RDA2 (",rda2_varex," %)"))+
        theme_bw()+
        xlim(-1.5, 3.5)+
        ylim(-2.5, 2.5)
      #geom_text_npc(aes(npcx = "left", npcy = "top", label = paste("AdjR2", round(global_r2,2), "PVal", pval)))
      
      
      print(p)
      
    } else {
      
      
      #plot with ggplot
      p = ggplot(fmod, aes(x = CAP1, y = MDS1))+
        geom_point(data = cbind(subset(fmod, Score == "sites"), Location = data.sub$Location),
                   aes(colour = Location), size = 3) +
        scale_colour_brewer("Location", palette = "Set1") +
        coord_fixed() +
        theme(legend.position = "top")+
        geom_segment(data=df_environ_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CAP1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=MDS1*scaling_factor), #Ending coordinate in CA1 
                     #color="firebrick1", #set color
                     linewidth=1,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        geom_segment(data=df_environ_no_sig, 
                     aes(x=0, #Starting coordinate in CCA1 = 0 
                         xend=CAP1*scaling_factor,#Ending coordinate in CCA1  
                         y=0, #Start in CCA2 = 0
                         yend=MDS1*scaling_factor), #Ending coordinate in CA1 
                     #color="firebrick1", #set color
                     linewidth=0.2,
                     arrow=arrow(length=unit(0.01,"npc"))#Set the size of the lines that form the tip of the arrow
        )+
        #Add environmental vars text
        geom_text(data=df_environ, 
                  aes(x=CAP1*scaling_factor, 
                      y=MDS1*scaling_factor,
                      label=rownames(df_environ),
                      hjust=0.5*(1-sign(CAP1)),#Add the text of each environmental var at the end of the arrow
                      vjust=0.5*(1-sign(MDS1))),#Add the text of each environmental var at the end of the arrow 
                  #color="firebrick1"
                  size=4,
        )+
        labs(x=paste0("RDA1 (",rda1_varex," %)"),
             y=paste0("MDS1"))+
        theme_bw()+
        xlim(min(fmod$CAP1),max(fmod$CAP1)+1)+
        ggtitle(habitat)#+
      #geom_text_npc(aes(npcx = "left", npcy = "top", label = paste("AdjR2", round(global_r2,2), "PVal", pval)))
      
      print(p)
      
      
    }
    
    
    return(list(rda.final.df, p, pval, global_r2))
    
  } else if (all(rda.final.df$p > 0.05, na.rm=TRUE)){
    
    return(rda.final.df)
    
  }
  
  
  
}

```

```{r, echo=FALSE}
#prokaryotes
rda.pro.sn <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow", col="#fa9f99")
rda.pro.sn[[1]] #results
rda.pro.sn[[2]] #plot
rda.pro.sn[[3]] #pval
rda.pro.sn[[4]] #adjusted r
rda.pro.sp <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice", col="#a4c64d")
rda.pro.sp[[1]] #results
rda.pro.sp[[2]] #plot
rda.pro.sp[[3]] #pval
rda.pro.sp[[4]] #adjusted r
rda.pro.sm <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice", col="#4dd2d6")
rda.pro.sm[[1]] #results
rda.pro.sm[[2]] #plot
rda.pro.sm[[3]] #pval
rda.pro.sm[[4]] #adjusted r
rda.pro.cr <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite", col="#d8a4ff")
rda.pro.cr[[1]] #results
rda.pro.cr[[2]] #plot
rda.pro.cr[[3]] #pval
rda.pro.cr[[4]] #adjusted r

#eukaryotes
rda.euk.sn <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Eukaryote", "Snow", col="#fa9f99")
rda.euk.sn[[1]] #results
rda.euk.sn[[2]] #plot
rda.euk.sn[[3]] #pval
rda.euk.sn[[4]] #adjusted r
rda.euk.sp <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Eukaryote", "Spring Ice", col="#a4c64d")
rda.euk.sp[[1]] #results
rda.euk.sp[[2]] #plot
rda.euk.sp[[3]] #pval
rda.euk.sp[[4]] #adjusted r
rda.euk.sm <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Eukaryote", "Summer Ice", col="#4dd2d6")
rda.euk.sm[[1]] #table only
rda.euk.cr <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Eukaryote", "Cryoconite", col="#d8a4ff")
rda.euk.cr[[1]] #results
rda.euk.cr[[2]] #plot
rda.euk.cr[[3]] #pval
rda.euk.cr[[4]] #adjusted r

#micrometazoans
rda.mm.sn <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Micrometazoan", "Snow", col="#fa9f99")
rda.mm.sn[[1]] #table only
rda.mm.sp <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Micrometazoan", "Spring Ice", col="#a4c64d")
rda.mm.sp[[1]] #table only
rda.mm.sm <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Micrometazoan", "Summer Ice", col="#4dd2d6")
rda.mm.sm[[1]] #results
rda.mm.sm[[2]] #plot
rda.mm.sm[[3]] #pval
rda.mm.sm[[4]] #adjusted r
rda.mm.cr <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Micrometazoan", "Cryoconite", col="#d8a4ff")
rda.mm.cr[[1]] #results
rda.mm.cr[[2]] #plot
rda.mm.cr[[3]] #pval
rda.mm.cr[[4]] #adjusted r
```


```{r, echo=FALSE}
#table order of variables
geochem.ord = c("pH", "Conductivity", "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
                "Cl", "Br", "Fl", "Fe", "Na", "Mg", "K", "Ca", "Li", "Al", "Ti", "V", "Cr", "Mn", "Co",
                "Ni", "Cu", "Zn", "Y", "Rb", "Sr", "Mo", "Ag", "Ba", "Lu", "Pb", "Zr", "Cd", "As", "Se", "Nb", "Sn", "Re", "U", "La", "Residual")

#bind into one dataframe
names(rda.pro.sn[[1]])[2:4] = paste(names(rda.pro.sn[[1]])[2:4], "16S-Snow")
names(rda.pro.sp[[1]])[2:4] = paste(names(rda.pro.sp[[1]])[2:4], "16S-Spring")
names(rda.pro.sm[[1]])[2:4] = paste(names(rda.pro.sm[[1]])[2:4], "16S-Summer")
names(rda.pro.cr[[1]])[2:4] = paste(names(rda.pro.cr[[1]])[2:4], "16S-Cryo")

#use dplyr to join dfs
rda.pro.df = full_join(rda.pro.sn[[1]][,1:4], rda.pro.sp[[1]][,1:4], by = 'Variables')
rda.pro.df = full_join(rda.pro.df, rda.pro.sm[[1]][,1:4], by = 'Variables')
rda.pro.df = full_join(rda.pro.df, rda.pro.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(rda.pro.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, rda.pro.df$Variables)]
x = rda.pro.df[match(geo.keep, rda.pro.df$Variables),]

write.csv(x, "../results-proportional/db-rda-pro-df.csv")

#bind into one dataframe
names(rda.euk.sn[[1]])[2:4] = paste(names(rda.euk.sn[[1]])[2:4], "18S-Snow")
names(rda.euk.sp[[1]])[2:4] = paste(names(rda.euk.sp[[1]])[2:4], "18S-Spring")
names(rda.euk.sm)[2:4] = paste(names(rda.euk.sm)[2:4], "18S-Summer")
names(rda.euk.cr[[1]])[2:4] = paste(names(rda.euk.cr[[1]])[2:4], "18S-Cryo")

#use dplyr to join dfs
rda.euk.df = full_join(rda.euk.sn[[1]][,1:4], rda.euk.sp[[1]][,1:4], by = 'Variables')
rda.euk.df = full_join(rda.euk.df, rda.euk.sm[,1:4], by = 'Variables')
rda.euk.df = full_join(rda.euk.df, rda.euk.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(rda.euk.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, rda.euk.df$Variables)]
y = rda.euk.df[match(geo.keep, rda.euk.df$Variables),]

write.csv(y, "../results-proportional/db-rda-euk-df.csv")

#bind into one dataframe - micrometazoans
names(rda.mm.sn)[2:4] = paste(names(rda.mm.sn)[2:4], "MM-Snow")
names(rda.mm.sp)[2:4] = paste(names(rda.mm.sp)[2:4], "MM-Spring")
names(rda.mm.sm[[1]])[2:4] = paste(names(rda.mm.sm[[1]])[2:4], "MM-Summer")
names(rda.mm.cr[[1]])[2:4] = paste(names(rda.mm.cr[[1]])[2:4], "MM-Cryo")

#use dplyr to join dfs
rda.mm.df = full_join(rda.mm.sn[,1:4], rda.mm.sp[,1:4], by = 'Variables')
rda.mm.df = full_join(rda.mm.df, rda.mm.sm[[1]][,1:4], by = 'Variables')
rda.mm.df = full_join(rda.mm.df, rda.mm.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(rda.mm.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, rda.mm.df$Variables)]
z = rda.mm.df[match(geo.keep, rda.mm.df$Variables),]

write.csv(z, "../results-proportional/db-rda-mm-df.csv")
```


```{r, echo=FALSE}
#Which plots do I want to present?
#CCA
#prokaryotes: spring ice and cryoconite
#eukaryotes: cryoconite
#micrometazoa: cryoconite

#RDA
#prokaryotes: spring ice and cryoconite
#eukaryotes: spring ice and cryoconite
#micrometazoa: cryoconite

#get legend
#creat NMDS plot so we can extract the apporpriate legend
ps = subset_samples(ps.pro, Habitat %in% c("Spring Ice", "Cryoconite"))
#set habitat location variable
sample_data(ps)$Habitat_location <- paste0(sample_data(ps)$Habitat, " (", sample_data(ps)$Location, " Foxfonna)")
sample_data(ps)$Location <- paste0(sample_data(ps)$Location, " Foxfonna")
#perform NMDS
pro.nmds.ord <- phyloseq::ordinate(ps, "NMDS", distance="bray")
#get stress
pro.nmds.stress = round(pro.nmds.ord$stress, 3)
#extract data for plotting
positions <- pro.nmds.ord$points
#get sample data for plotting
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps))$Habitat_location, data.frame(sample_data(ps))$Habitat, data.frame(sample_data(ps))$Location, "Prokaryote")
names(data2plot) = c("Sample", "NMDS1", "NMDS2", "Habitat Location", "Habitat", "Location", "Group")
#default legend order
data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Spring Ice", "Cryoconite"))
#plot
pro.nmds.p = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Habitat, shape=Location),
             alpha=.5, size=4)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("#a4c64d", "#d8a4ff"))


#get legend
legend <- get_legend(
  pro.nmds.p  +
    #guides(color = guide_legend(nrow = 1, override.aes = list(size = 10))) +
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    #guides(color=guide_legend(nrow=2,byrow=TRUE))+
    theme(legend.box="vertical",
          axis.text=element_text(size=15),
          axis.title=element_text(size=15,face="bold"),
          legend.text = element_text(size=15),
          legend.title = element_blank(), legend.margin=margin(),
          legend.position=c(.5,.6))
  
)



#prep graphs for plotting
p1 = cca.pro.sp[[2]] + theme(legend.position="none") + ggtitle("Prokaryote")
p2 = cca.pro.cr[[2]] + theme(legend.position="none") + ggtitle("Prokaryote")
p3 = cca.euk.cr[[2]] + theme(legend.position="none") + ggtitle("Microbial eukaryote")
p4 = cca.mm.cr[[2]] + theme(legend.position="none") + ggtitle("Microfauna")
p5 = rda.pro.sp[[2]] + theme(legend.position="none") + ggtitle("Prokaryote")
p6 = rda.pro.cr[[2]] + theme(legend.position="none") + ggtitle("Prokaryote")
p7 = rda.euk.sp[[2]] + theme(legend.position="none") + ggtitle("Microbial eukaryote")
p8 = rda.euk.cr[[2]] + theme(legend.position="none") + ggtitle("Microbial eukaryote")
p9 = rda.mm.cr[[2]] + theme(legend.position="none") + ggtitle("Microfauna")


############plot grid#############
#title1 <- ggdraw() + draw_label("Prokaryote", fontface='bold')
cca.plot.1 = plot_grid(p1, p2, ncol=2)
#cca.plot.1 = plot_grid(title1, cca.plot.1, rel_heights = c(0.1, 1), ncol=1)
cca.plot.1

#title2 <- ggdraw() + draw_label("Microbial eukaryote", fontface='bold', hjust = 1.65)
# cca.plot.2 = plot_grid(p3, p4 rel_heights=c(0.1, 1), ncol=1)
# cca.plot.2
cca.plot.2 = plot_grid(p3, p4, ncol=2)
#cca.plot.2 = plot_grid(title2, cca.plot.2, rel_heights = c(0.1, 1), ncol=1)
cca.plot.2

cca.plot.3 = plot_grid(cca.plot.1, cca.plot.2, ncol=1)
cca.plot.3

#save plot
pdf("../results-proportional/cca-plot-1.pdf")
print(cca.plot.3)
dev.off()


rda.plot.1 = plot_grid(p5, p6, ncol=2)
rda.plot.1
rda.plot.2 = plot_grid(p7, p8, ncol=2)
rda.plot.2
rda.plot.3 = plot_grid(p9, NULL, ncol=2)
rda.plot.3

rda.plot.4 = plot_grid(rda.plot.1, rda.plot.2, ncol=1)
rda.plot.4

#save plot
pdf("../results-proportional/rda-plot-1.pdf")
print(rda.plot.4)
dev.off()

mm.plot = plot_grid(p9, NULL, legend, NULL, ncol=2)
mm.plot

final.plot = plot_grid(cca.plot.3, rda.plot.4, mm.plot, ncol=3,labels=c("A", "B", ""))

#final.plot = plot_grid(cca.plot.3, rda.plot.4, ncol=2,labels=c("A", "B", ""))

#save plot
pdf("../results-proportional/cca-rda-plot-all.pdf", width=15, height=6)
print(final.plot)
dev.off()
```

