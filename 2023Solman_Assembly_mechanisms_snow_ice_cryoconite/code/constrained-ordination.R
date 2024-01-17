# Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
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
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without Microfaunas
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#Microfaunas only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 


#CCA
#Methods for dealing with high VIF values
#https://quantpalaeo.wordpress.com/2014/04/14/variance-inflation-factors-and-ordination-model-selection/
  
  #Tutorial for CCA
#   https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html
# https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html

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
  
  #reduce the number of variables for Microfaunas because we have less data points
  if(group == "Microfauna"){
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
  
  #if we have variables included in our model add them to the dataframe
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
  
  #if any of those terms are significant...
  if (any(cca.final.df$p < 0.05, na.rm=TRUE)){
    
    #get the model data
    fmod <- fortify(ccamodel2) 
    
    #get site data
    data = data.frame(sample_data(s))
    data$Location = paste(data$Location, "Foxfonna")
    
    #subset to include only those samples kept in the model
    data.sub = data[data$SampleID %in% subset(fmod, score == "sites")$label,]
    
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
        geom_point(data = cbind(subset(fmod, score == "sites"), Location = data.sub$Location),
                   aes(fill = col, shape=Location, alpha=0.7), size = 4)+
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
        geom_point(data = cbind(subset(fmod, score == "sites"), Location = data.sub$Location),
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
    
    return(list(cca.final.df, NULL, pval, global_r2))
    
  }
  
  
  
}

#RDA

#Tutorial for RDA
# https://rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html
# https://r.qcbs.ca/workshop10/book-en/partial-redundancy-analysis.html

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
  
  #reduce the number of variables for Microfaunas because we have less data points
  if(group == "Microfauna"){
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
    data.sub = data[data$SampleID %in% subset(fmod, score == "sites")$label,]
    
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
        geom_point(data = cbind(subset(fmod, score == "sites"), Location = data.sub$Location),
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
        geom_point(data = cbind(subset(fmod, score == "sites"), Location = data.sub$Location),
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
    
    return(list(rda.final.df, NULL, pval, global_r2))
    
  }
  
  
  
}

######################################RUN THE CCA MODEL######################################################################

#prokaryotes
cca.pro.sn <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow", col="#fa9f99")
cca.pro.sp <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice", col="#a4c64d")
cca.pro.sm <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice", col="#4dd2d6")
cca.pro.cr <- cca_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite", col="#d8a4ff")

#eukaryotes
cca.euk.sn <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Microbial eukaryote", "Snow", col="#fa9f99")
cca.euk.sp <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Microbial eukaryote", "Spring Ice", col="#a4c64d")
cca.euk.sm <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Microbial eukaryote", "Summer Ice", col="#4dd2d6")
cca.euk.cr <- cca_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Microbial eukaryote", "Cryoconite", col="#d8a4ff")

#microfauna
cca.mm.sn <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Microfauna", "Snow", col="#fa9f99")
cca.mm.sp <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Microfauna", "Spring Ice", col="#a4c64d")
cca.mm.sm <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Microfauna", "Summer Ice", col="#4dd2d6")
cca.mm.cr <- cca_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Microfauna", "Cryoconite", col="#d8a4ff")


######################################RUN THE DB-RDA MODEL#############################

#prokaryotes
rda.pro.sn <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow", col="#fa9f99")
rda.pro.sp <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice", col="#a4c64d")
rda.pro.sm <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice", col="#4dd2d6")
rda.pro.cr <- dbRDA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite", col="#d8a4ff")

#eukaryotes
rda.euk.sn <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Eukaryote", "Snow", col="#fa9f99")
rda.euk.sp <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Eukaryote", "Spring Ice", col="#a4c64d")
rda.euk.sm <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Eukaryote", "Summer Ice", col="#4dd2d6")
rda.euk.cr <- dbRDA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Eukaryote", "Cryoconite", col="#d8a4ff")

#Microfauna
rda.mm.sn <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Microfauna", "Snow", col="#fa9f99")
rda.mm.sp <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Microfauna", "Spring Ice", col="#a4c64d")
rda.mm.sm <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Microfauna", "Summer Ice", col="#4dd2d6")
rda.mm.cr <- dbRDA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Microfauna", "Cryoconite", col="#d8a4ff")


######################################PUT ALL THE CCA RESULTS TOGETHER INTO A DATAFRAME#######################################

#table order of variables
geochem.ord = c("pH", "Conductivity", "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
                "Cl", "Br", "Fl", "Fe", "Na", "Mg", "K", "Ca", "Li", "Al", "Ti", "V", "Cr", "Mn", "Co",
                "Ni", "Cu", "Zn", "Y", "Rb", "Sr", "Mo", "Ag", "Ba", "Lu", "Pb", "Zr", "Cd", "As", "Se", "Nb", "Sn", "Re", "U", "La", "Residual")

#bind into one dataframe
names(cca.pro.sn[[1]])[2:4] = paste(names(cca.pro.sn[[1]])[2:4], "16S-Snow")
names(cca.pro.sp[[1]])[2:4] = paste(names(cca.pro.sp[[1]])[2:4], "16S-Spring")
names(cca.pro.sm[[1]])[2:4] = paste(names(cca.pro.sm[[1]])[2:4], "16S-Summer")
names(cca.pro.cr[[1]])[2:4] = paste(names(cca.pro.cr[[1]])[2:4], "16S-Cryo")

#use dplyr to join dfs
cca.pro.df = full_join(cca.pro.sn[[1]][,1:4], cca.pro.sp[[1]][,1:4], by = 'Variables')
cca.pro.df = full_join(cca.pro.df, cca.pro.sm[[1]][,1:4], by = 'Variables')
cca.pro.df = full_join(cca.pro.df, cca.pro.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.pro.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.pro.df$Variables)]
x = cca.pro.df[match(geo.keep, cca.pro.df$Variables),]

write.csv(x, "../results/cca-pro-df.csv")

#bind into one dataframe
names(cca.euk.sn[[1]])[2:4] = paste(names(cca.euk.sn[[1]])[2:4], "18S-Snow")
names(cca.euk.sp[[1]])[2:4] = paste(names(cca.euk.sp[[1]])[2:4], "18S-Spring")
names(cca.euk.sm[[1]])[2:4] = paste(names(cca.euk.sm[[1]])[2:4], "18S-Summer")
names(cca.euk.cr[[1]])[2:4] = paste(names(cca.euk.cr[[1]])[2:4], "18S-Cryo")

#use dplyr to join dfs
cca.euk.df = full_join(cca.euk.sn[[1]][,1:4], cca.euk.sp[[1]][,1:4], by = 'Variables')
cca.euk.df = full_join(cca.euk.df, cca.euk.sm[[1]][,1:4], by = 'Variables')
cca.euk.df = full_join(cca.euk.df, cca.euk.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.euk.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.euk.df$Variables)]
y = cca.euk.df[match(geo.keep, cca.euk.df$Variables),]

write.csv(y, "../results/cca-euk-df.csv")

#bind into one dataframe - Microfaunas
names(cca.mm.sn[[1]])[2:4] = paste(names(cca.mm.sn[[1]])[2:4], "MM-Snow")
names(cca.mm.sp[[1]])[2:4] = paste(names(cca.mm.sp[[1]])[2:4], "MM-Spring")
names(cca.mm.sm[[1]])[2:4] = paste(names(cca.mm.sm[[1]])[2:4], "MM-Summer")
names(cca.mm.cr[[1]])[2:4] = paste(names(cca.mm.cr[[1]])[2:4], "MM-Cryo")

#use dplyr to join dfs
cca.mm.df = full_join(cca.mm.sn[[1]][,1:4], cca.mm.sp[[1]][,1:4], by = 'Variables')
cca.mm.df = full_join(cca.mm.df, cca.mm.sm[[1]][,1:4], by = 'Variables')
cca.mm.df = full_join(cca.mm.df, cca.mm.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(cca.mm.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, cca.mm.df$Variables)]
z = cca.mm.df[match(geo.keep, cca.mm.df$Variables),]

write.csv(z, "../results/cca-mm-df.csv")



######################################PUT ALL THE DB-RDA RESULTS TOGETHER INTO A DATAFRAME#############################

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

write.csv(x, "../results/db-rda-pro-df.csv")

#bind into one dataframe
names(rda.euk.sn[[1]])[2:4] = paste(names(rda.euk.sn[[1]])[2:4], "18S-Snow")
names(rda.euk.sp[[1]])[2:4] = paste(names(rda.euk.sp[[1]])[2:4], "18S-Spring")
names(rda.euk.sm[[1]])[2:4] = paste(names(rda.euk.sm[[1]])[2:4], "18S-Summer")
names(rda.euk.cr[[1]])[2:4] = paste(names(rda.euk.cr[[1]])[2:4], "18S-Cryo")

#use dplyr to join dfs
rda.euk.df = full_join(rda.euk.sn[[1]][,1:4], rda.euk.sp[[1]][,1:4], by = 'Variables')
rda.euk.df = full_join(rda.euk.df, rda.euk.sm[[1]][,1:4], by = 'Variables')
rda.euk.df = full_join(rda.euk.df, rda.euk.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(rda.euk.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, rda.euk.df$Variables)]
y = rda.euk.df[match(geo.keep, rda.euk.df$Variables),]

write.csv(y, "../results/db-rda-euk-df.csv")

#bind into one dataframe - Microfaunas
names(rda.mm.sn[[1]])[2:4] = paste(names(rda.mm.sn[[1]])[2:4], "MM-Snow")
names(rda.mm.sp[[1]])[2:4] = paste(names(rda.mm.sp[[1]])[2:4], "MM-Spring")
names(rda.mm.sm[[1]])[2:4] = paste(names(rda.mm.sm[[1]])[2:4], "MM-Summer")
names(rda.mm.cr[[1]])[2:4] = paste(names(rda.mm.cr[[1]])[2:4], "MM-Cryo")

#use dplyr to join dfs
rda.mm.df = full_join(rda.mm.sn[[1]][,1:4], rda.mm.sp[[1]][,1:4], by = 'Variables')
rda.mm.df = full_join(rda.mm.df, rda.mm.sm[[1]][,1:4], by = 'Variables')
rda.mm.df = full_join(rda.mm.df, rda.mm.cr[[1]][,1:4], by = 'Variables')

#order the variables
setdiff(rda.mm.df$Variables, geochem.ord) #should return nothing
geo.keep = geochem.ord[!geochem.ord %in% setdiff(geochem.ord, rda.mm.df$Variables)]
z = rda.mm.df[match(geo.keep, rda.mm.df$Variables),]

write.csv(z, "../results/db-rda-mm-df.csv")

######################################REPORT RESULTS#########################################################


sink("../results/cca-db-rda.txt", type="output")
writeLines("===============================================================
CONSTRAINED ORDINATION RESULTS
===============================================================")
writeLines("Overall performance of CCA model for prokaryotes in snow")
writeLines("Adj R2")
cca.pro.sn[[4]]
writeLines("P val")
cca.pro.sn[[3]]
writeLines("Significant terms of the model")
cca.pro.sn[[1]][cca.pro.sn[[1]]$`p 16S-Snow` < 0.05 ,]

writeLines("Overall performance of CCA model for prokaryotes in spring ice")
writeLines("Adj R2")
cca.pro.sp[[4]]
writeLines("P val")
cca.pro.sp[[3]]
writeLines("Significant terms of the model")
cca.pro.sp[[1]][cca.pro.sp[[1]]$`p 16S-Spring` < 0.05 ,]

writeLines("Overall performance of CCA model for prokaryotes in summer ice")
writeLines("Adj R2")
cca.pro.sm[[4]]
writeLines("P val")
cca.pro.sm[[3]]
writeLines("Significant terms of the model")
cca.pro.sm[[1]][cca.pro.sm[[1]]$`p 16S-Summer` < 0.05 ,]

writeLines("Overall performance of CCA model for prokaryotes in cryoconite")
writeLines("Adj R2")
cca.pro.cr[[4]]
writeLines("P val")
cca.pro.cr[[3]]
writeLines("Significant terms of the model")
cca.pro.cr[[1]][cca.pro.cr[[1]]$`p 16S-Cryo` < 0.05 ,]

writeLines("Overall performance of CCA model for eukaryotes in snow")
writeLines("Adj R2")
cca.euk.sn[[4]]
writeLines("P val")
cca.euk.sn[[3]]
writeLines("Significant terms of the model")
cca.euk.sn[[1]][cca.euk.sn[[1]]$`p 18S-Snow` < 0.05 ,]

writeLines("Overall performance of CCA model for eukaryotes in spring ice")
writeLines("Adj R2")
cca.euk.sp[[4]]
writeLines("P val")
cca.euk.sp[[3]]
writeLines("Significant terms of the model")
cca.euk.sp[[1]][cca.euk.sp[[1]]$`p 18S-Spring` < 0.05 ,]

writeLines("Overall performance of CCA model for eukaryotes in summer ice")
writeLines("Adj R2")
cca.euk.sm[[4]]
writeLines("P val")
cca.euk.sm[[3]]
writeLines("Significant terms of the model")
cca.euk.sm[[1]][cca.euk.sm[[1]]$`p 18S-Summer` < 0.05 ,]

writeLines("Overall performance of CCA model for eukaryotes in cryoconite")
writeLines("Adj R2")
cca.euk.cr[[4]]
writeLines("P val")
cca.euk.cr[[3]]
writeLines("Significant terms of the model")
cca.euk.cr[[1]][cca.euk.cr[[1]]$`p 18S-Cryo` < 0.05 ,]

writeLines("Overall performance of CCA model for microfauna in snow")
writeLines("Adj R2")
cca.mm.sn[[4]]
writeLines("P val")
cca.mm.sn[[3]]
writeLines("Significant terms of the model")
cca.mm.sn[[1]][cca.mm.sn[[1]]$`p MM-Snow` < 0.05 ,]

writeLines("Overall performance of CCA model for microfauna in spring ice")
writeLines("Adj R2")
cca.mm.sp[[4]]
writeLines("P val")
cca.mm.sp[[3]]
writeLines("Significant terms of the model")
cca.mm.sp[[1]][cca.mm.sp[[1]]$`p MM-Spring` < 0.05 ,]

writeLines("Overall performance of CCA model for microfauna in summer ice")
writeLines("Adj R2")
cca.mm.sm[[4]]
writeLines("P val")
cca.mm.sm[[3]]
writeLines("Significant terms of the model")
cca.mm.sm[[1]][cca.mm.sm[[1]]$`p MM-Summer` < 0.05 ,]

writeLines("Overall performance of CCA model for microfauna in cryoconite")
writeLines("Adj R2")
cca.mm.cr[[4]]
writeLines("P val")
cca.mm.cr[[3]]
writeLines("Significant terms of the model")
cca.mm.cr[[1]][cca.mm.cr[[1]]$`p MM-Cryo` < 0.05 ,]

writeLines("Overall performance of RDA model for prokaryotes in snow")
writeLines("Adj R2")
rda.pro.sn[[4]]
writeLines("P val")
rda.pro.sn[[3]]
writeLines("Significant terms of the model")
rda.pro.sn[[1]][rda.pro.sn[[1]]$`p 16S-Snow` < 0.05 ,]

writeLines("Overall performance of RDA model for prokaryotes in spring ice")
writeLines("Adj R2")
rda.pro.sp[[4]]
writeLines("P val")
rda.pro.sp[[3]]
writeLines("Significant terms of the model")
rda.pro.sp[[1]][rda.pro.sp[[1]]$`p 16S-Spring` < 0.05 ,]

writeLines("Overall performance of RDA model for prokaryotes in summer ice")
writeLines("Adj R2")
rda.pro.sm[[4]]
writeLines("P val")
rda.pro.sm[[3]]
writeLines("Significant terms of the model")
rda.pro.sm[[1]][rda.pro.sm[[1]]$`p 16S-Summer` < 0.05 ,]

writeLines("Overall performance of RDA model for prokaryotes in cryoconite")
writeLines("Adj R2")
rda.pro.cr[[4]]
writeLines("P val")
rda.pro.cr[[3]]
writeLines("Significant terms of the model")
rda.pro.cr[[1]][rda.pro.cr[[1]]$`p 16S-Cryo` < 0.05 ,]

writeLines("Overall performance of RDA model for eukaryotes in snow")
writeLines("Adj R2")
rda.euk.sn[[4]]
writeLines("P val")
rda.euk.sn[[3]]
writeLines("Significant terms of the model")
rda.euk.sn[[1]][rda.euk.sn[[1]]$`p 18S-Snow` < 0.05 ,]

writeLines("Overall performance of rda model for eukaryotes in spring ice")
writeLines("Adj R2")
rda.euk.sp[[4]]
writeLines("P val")
rda.euk.sp[[3]]
writeLines("Significant terms of the model")
rda.euk.sp[[1]][rda.euk.sp[[1]]$`p 18S-Spring` < 0.05 ,]

writeLines("Overall performance of rda model for eukaryotes in summer ice")
writeLines("Adj R2")
rda.euk.sm[[4]]
writeLines("P val")
rda.euk.sm[[3]]
writeLines("Significant terms of the model")
rda.euk.sm[[1]][rda.euk.sm[[1]]$`p 18S-Summer` < 0.05 ,]

writeLines("Overall performance of rda model for eukaryotes in cryoconite")
writeLines("Adj R2")
rda.euk.cr[[4]]
writeLines("P val")
rda.euk.cr[[3]]
writeLines("Significant terms of the model")
rda.euk.cr[[1]][rda.euk.cr[[1]]$`p 18S-Cryo` < 0.05 ,]

writeLines("Overall performance of rda model for microfauna in snow")
writeLines("Adj R2")
rda.mm.sn[[4]]
writeLines("P val")
rda.mm.sn[[3]]
writeLines("Significant terms of the model")
rda.mm.sn[[1]][rda.mm.sn[[1]]$`p MM-Snow` < 0.05 ,]

writeLines("Overall performance of rda model for microfauna in spring ice")
writeLines("Adj R2")
rda.mm.sp[[4]]
writeLines("P val")
rda.mm.sp[[3]]
writeLines("Significant terms of the model")
rda.mm.sp[[1]][rda.mm.sp[[1]]$`p MM-Spring` < 0.05 ,]

writeLines("Overall performance of rda model for microfauna in summer ice")
writeLines("Adj R2")
rda.mm.sm[[4]]
writeLines("P val")
rda.mm.sm[[3]]
writeLines("Significant terms of the model")
rda.mm.sm[[1]][rda.mm.sm[[1]]$`p MM-Summer` < 0.05 ,]

writeLines("Overall performance of rda model for microfauna in cryoconite")
writeLines("Adj R2")
rda.mm.cr[[4]]
writeLines("P val")
rda.mm.cr[[3]]
writeLines("Significant terms of the model")
rda.mm.cr[[1]][rda.mm.cr[[1]]$`p MM-Cryo` < 0.05 ,]

sink()




######################################PLOT#########################################################


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

pro.nmds.p


#get legend
leg <- get_legend(
  pro.nmds.p  +
    #guides(color = guide_legend(nrow = 1, override.aes = list(size = 10))) +
    guides(fill=guide_legend(override.aes=list(shape=22)))+
    #guides(color=guide_legend(nrow=2,byrow=TRUE))+
    theme_bw()+
    theme(legend.box="horizontal",
          legend.text = element_text(size=12),
          legend.title = element_blank(), legend.margin=margin(),
          legend.position=c(.5,.6))
  
)



#prep graphs for plotting
p1 = rda.pro.sp[[2]] + theme(legend.position="none", plot.title = element_text(face="bold", size=10)) + xlim(-1.5,2.5) + ylim(-2,2) + ggtitle("Prokaryote")
p2 = rda.pro.cr[[2]] + theme(legend.position="none") + xlim(-1.5,2.5) + ylim(-2,2)
p3 = rda.euk.sp[[2]] + theme(legend.position="none", plot.title = element_text(face="bold", size=10)) + xlim(-1.5,2.5) + ylim(-2,2) + ggtitle("Microbial Eukaryote")
p4 = rda.euk.cr[[2]] + theme(legend.position="none") + xlim(-1.5,2.5) + ylim(-2,2)


final.p = (p1 + p2) / (p3 + p4) / plot_spacer() / leg + plot_layout(heights = c(1, 1, 0.05, 0.1))
final.p

#save plot
pdf("../results/rda-plot.pdf", width=6, height=7)
print(final.p)
dev.off()



