# 1. Clear workspace and load packages
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
library(ggvenn)
library(ggplot2)
library(ggforce)
library(patchwork)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without Microfaunas
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#Microfaunas only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#VPA + Mantel tests
#Excellent tutorial: https://www.davidzeleny.net/anadat-r/doku.php/en:varpart_examples

# ps = ps.euk
# hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID
# group = "Eukaryote"
# habitat = "Cryoconite"

VPA_func <- function(ps, hab, group, habitat){
  
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
  
  #remove columns with less than 50% values != 0 or NA
  c <- vars
  c[c > 0 | c < 0] <- 1 #give values greater than 0 or less than 0 the value of 1
  c[is.na(c)] <- 0 #make NA 0
  var.trim <- vars[,colSums(c) >= round(nrow(vars)*0.50)] #only keep those with at least 50% values != 0 or NA
  #only keep complete cases
  var.comp <- var.trim[complete.cases(var.trim),]
  
  #reduce the number of variables for Microfaunas because we have less data points
  if(group == "Microfauna"){
    var.comp = var.comp[,names(var.comp) %in% c("pH", "Conductivity", "TC", "Br", "Fe", "Na", "Mg", "K", "Ca")]
  }
  
  #z-score transform the environmental data
  var.z = data.frame(scale(var.comp))
  
  #make sure our samples are the same for comm and env
  comm = comm[rownames(comm) %in% rownames(var.z),]
  
  #############CARRY OUT DISTANCE BASED RDA TO CHECK FOR COLLINEAR ENVIRONMENTAL VARIABLES AND REMOVE THEM####################
  
  #perform db-RDA
  rdamodel1 <- capscale(comm ~., var.z, dist="bray", add=TRUE) #model with all explanatory variables
  
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
  
  #WHICH ENVIRONMENTAL VARIABLES ARE SIGNIFICANT?
  rdamodel2 <- capscale(comm ~., var.z.trim, dist="bray", add=TRUE)
  sig = anova.cca(rdamodel2, by="terms")
  #grab significant terms from model anova results
  env.keep = rownames(sig[sig$`Pr(>F)` < 0.05,])
  env.keep = env.keep[env.keep != "NA"]
  env.sel = data.frame(var.z.trim[,names(var.z.trim) %in% env.keep]) #subset our data to the significant terms
  #names(env.sel) = env.keep
  
  #GET OUR GEOGRAPHIC FACTORS
  #get distance data
  lat.long = data.frame(sample_data(s))[,c(7:8)]
  #make sure our samples are the same as comm samples
  geo = lat.long[rownames(lat.long) %in% rownames(comm),]
  
  #Calculate PCNMs from a Euclidean distance matrix of sample coordinates and extract scores associated with these new variables.
  geo.pcnm <- as.data.frame(scores(pcnm(dist(geo))))
  
  # Use ordistep to select significant PCNM axes from db-RDA
  fwd.sel.geo <- ordistep(capscale(comm ~ 1, geo.pcnm, dist="bray", add=TRUE),
                          scope = formula(capscale(comm ~ ., geo.pcnm, dist="bray", add=TRUE)),
                          direction="forward")
  
  #Look at the new model with forward selected variables
  fwd.sel.geo$call
  
  #grab significant terms from model anova results
  axis.keep = rownames(fwd.sel.geo$anova)
  axis.keep2 = sub('..', '', axis.keep) #remove bits we don't need
  geo.pcnm.sel = data.frame(geo.pcnm[,names(geo.pcnm) %in% axis.keep2]) #subset our pcnm axis to the significant terms
  names(geo.pcnm.sel) = axis.keep2
  
  #######################PERFORM VPA#######################
  
  if (length(axis.keep2) > 0 & ncol(env.sel) > 0){ #if we have significant spatial terms and env terms
    
    vpa <- varpart(vegdist(comm, method="bray"), env.sel, geo.pcnm.sel) #carry out VPA
    
    #total dataframe for env and geo varibles
    tot.df = cbind(geo.pcnm.sel, env.sel)
    
    ###################CONDITIONAL EFFECT OF ENVIRONMENTAL VARIABLES###########################
    #how much variance does environmental factors explain, controlling for geography?
    env.only.var.r2 = round(vpa$part$indfract$Adj.R.squared[[1]], 2)
    
    #is this variance significant?
    #do the same analysis but using capscale to test for significance
    #define partial ordination models
    db.rda.env.cont.geo = capscale(as.formula(sprintf("%s ~ %s + Condition(%s)", "comm", paste(env.sel, collapse = " + "), paste(geo.pcnm.sel, collapse = " + "))), data=tot.df, dist="bray", add=TRUE)
    round(RsquareAdj(db.rda.env.cont.geo)$adj.r.squared,2) == env.only.var.r2 #should equal TRUE because it's the same R2 result
    env.only.var.sig = anova(db.rda.env.cont.geo)$`Pr(>F)`[1]
    
    ###################SIMPLE EFFECT OF ENVIRONMENTAL VARIABLES###########################
    #how much variance do environmental factors explain without controlling for geography?
    env.var.r2 = round(vpa$part$fract$Adj.R.squared[[1]],2)
    
    #is this variance significant?
    #do the same analysis but using capscale to test for significance
    db.rda.env = capscale(comm ~ ., env.sel, dist="bray", add=TRUE)
    round(RsquareAdj(db.rda.env)$adj.r.squared,2)  == env.var.r2 #should be TRUE
    env.var.sig = anova(db.rda.env)$`Pr(>F)`[1]
    
    ###################CONDITIONAL EFFECT OF GEOGRAPHIC VARIABLES###########################
    #how much variance does geographic factors explain, controlling for environment?
    geo.only.var.r2 = round(vpa$part$indfract$Adj.R.squared[[2]], 2)
    
    #is this variance significant?
    #do the same analysis but using capscale to test for significance
    #define partial ordination models
    db.rda.geo.cont.env = capscale(as.formula(sprintf("%s ~ %s + Condition(%s)", "comm", paste(geo.pcnm.sel, collapse = " + "), paste(env.sel, collapse = " + "))), data=tot.df, dist="bray", add=TRUE)
    round(RsquareAdj(db.rda.geo.cont.env)$adj.r.squared,2)  == geo.only.var.r2 #should equal TRUE because it's the same R2 result
    geo.only.var.sig = anova(db.rda.geo.cont.env)$`Pr(>F)`[1]
    
    ###################SIMPLE EFFECT OF GEOGRAPHIC VARIABLES###########################
    #how much variance do geographic factors explain without controlling for environmental factors?
    geo.var.r2 = round(vpa$part$fract$Adj.R.squared[[2]],2)
    
    #is this variance significant?
    #do the same analysis but using capscale to test for significance
    db.rda.geo = capscale(comm ~ ., geo.pcnm.sel, dist="bray", add=TRUE)
    round(RsquareAdj(db.rda.geo)$adj.r.squared,2)  == geo.var.r2 #should be TRUE
    geo.var.sig = anova(db.rda.geo)$`Pr(>F)`[1]
    
    ###################SHARED EFFECT OF ENVIRONMENTAL AND GEOGRAPHCI VARIABLES###########################
    #how much overlapping variance do environmental and geographic factors explain? This cannot be tested for significance
    env.geo.overlap.var.r2 = round(vpa$part$indfract$Adj.R.squared[[3]], 2)
    
    ###################TOTAL EFFECT OF ENVIRONMENTAL AND GEOGRAPHCI VARIABLES###########################
    #how much variance do environmental factors and geographic factors explain together?
    total.var.r2 = env.only.var.r2 + geo.only.var.r2 + env.geo.overlap.var.r2
    
    #is this variance significant?
    #do the same analysis but using capscale to test for significance
    db.rda.env.geo = capscale(comm ~ ., tot.df, dist="bray", add=TRUE)
    round(RsquareAdj(db.rda.env.geo)$adj.r.squared,2) == round(total.var.r2, 2)
    tot.var.sig = anova(db.rda.env.geo)$`Pr(>F)`[1]
    
    ##########################UNEXPLAINED VARIANCE############################
    #how much variance is unexplained?
    unex.var = 1-total.var.r2
    
    
    ##########################RESULTS TABLE############################
    
    #Results into a table
    res = data.frame(Term = c("Global", "Environmental (Marginal)", "Geographic (Marginal)", "Environmental (Partial)", "Geographic (Partial)"), 
                     AdjR2 = c(total.var.r2, env.var.r2, geo.var.r2, env.only.var.r2, geo.only.var.r2), 
                     PVal = c(tot.var.sig, env.var.sig, geo.var.sig, env.only.var.sig, geo.only.var.sig))
    res$Group = group
    res$Habitat = habitat
    
    #################ADDITIONAL MANTEL TESTS###############################
    
    rownames(comm) == rownames(env.sel) #should all be TRUE
    rownames(comm) == rownames(geo) #should all be TRUE
    
    #Get Bray-Curtis dissimilarity matrix for community data
    commdist<-vegdist(comm, method="bray")
    
    #Get get distance matrix for significant environmental factors
    #chemdist<-vegdist(var.z[,names(var.z) %in% env.keep], method="euclidean")
    #chemdist<-vegdist(var.z, method="euclidean")
    chemdist<-vegdist(env.sel, method="euclidean")
    
    #Get distance matrix for spatial factors
    geodist<-vegdist(geo, method="euclidean")
    
    
    # vegan::mantel(commdist, chemdist, method="pearson", permutations=999)
    # vegan::mantel.partial(commdist, chemdist, geodist, method="pearson", permutations=999)
    # 
    # mantel(vegdist(comm,'bray'),
    #        vegdist(scale(var.comp[,names(var.comp) %in% env.keep]),'euclidian'), 
    #        method="pearson", permutations=999)
    # 
    # mantel.partial(vegdist(comm,'bray'),
    #                vegdist(scale(var.comp[,names(var.comp) %in% env.keep]),'euclidian'), 
    #                vegdist(geo,'euclidian'), 
    #                method="pearson", permutations=999)
    # 
    # 
    # vegan::mantel(commdist, geodist, method="pearson", permutations=999)
    # mantel(vegdist(comm,'bray'),vegdist(geo,'euclidian'), method="pearson", permutations=999)
    
    #Mantel test for correlation between community composition and environmental factors
    comm_env <- vegan::mantel(commdist, chemdist, method="pearson", permutations=999)
    comm_env
    comm_env.r2 = round(comm_env$statistic, 2)
    comm_env.p = comm_env$signif
    
    #Partial Mantel test for correlation between community composition and environmental factors, controlling for spatial factors
    comm_env_geo <- vegan::mantel.partial(commdist, chemdist, geodist, method="pearson", permutations=999)
    comm_env_geo
    comm_env_geo.r2 = round(comm_env_geo$statistic, 2)
    comm_env_geo.p = comm_env_geo$signif
    
    #Mantel test for correlation between community composition and spatial factors.
    comm_geo <- vegan::mantel(commdist, geodist, method="pearson", permutations=999)
    comm_geo.r2 = round(comm_geo$statistic, 2)
    comm_geo.p = comm_geo$signif
    
    #Partial Mantel test for correlation between community composition and spatial factors, controlling for environmental factors
    comm_geo_env <- vegan::mantel.partial(commdist, geodist, chemdist, method="pearson", permutations=999)
    comm_geo_env.r2 = round(comm_geo_env$statistic, 2)
    comm_geo_env.p = comm_geo_env$signif
    
    #results dataframe
    res2 <- data.frame(Test = c("Environmental (Marginal)", "Geographic (Marginal)", "Environmental (Partial)", "Geographic (Partial)"), 
                       MantelR = c(comm_env.r2, comm_geo.r2, comm_env_geo.r2, comm_geo_env.r2), 
                       PVal = c(comm_env.p, comm_geo.p, comm_env_geo.p, comm_geo_env.p))
    res2$Group = group
    res2$Habitat = habitat
    
    return(list(res, vpa, res2))
    
  } else if (length(axis.keep2) == 0 & ncol(env.sel) > 0){ #if we have sig env terms only
    
    #we can look at our environmental variables only
    
    #test fraction A (environmental)
    db.rda.env = capscale(comm ~ ., env.sel, dist="bray", add=TRUE)
    a = RsquareAdj(db.rda.env)$adj.r.squared
    
    #TEST IF THESE VARIANCES ARE SIGNIFICANT
    
    #Simple (marginal) effect of env
    anova(db.rda.env)
    
    #Results into a table
    res = data.frame(Term = c("Environmental (Marginal)"), AdjR2 = c(a), PVal = c(anova(db.rda.env)$`Pr(>F)`[1]))
    res$AdjR2 = round(res$AdjR2, 3)
    res$Group = group
    res$Habitat = habitat
    
    return(list(res, res))
    
  } else if (length(axis.keep2) > 0 & ncol(env.sel) < 1){ #just test the geographic component
    
    #test geographic variables only
    db.rda.geo = capscale(comm ~ ., geo.pcnm.sel, dist="bray", add=TRUE)
    b = RsquareAdj(db.rda.geo)$adj.r.squared
    
    #TEST IF THESE VARIANCES ARE SIGNIFICAN
    anova(db.rda.geo)
    
    #Results into a table
    res = data.frame(Term = c("Geographic (Marginal)"), AdjR2 = c(b), PVal = c(anova(db.rda.geo)$`Pr(>F)`[1]))
    res$AdjR2 = round(res$AdjR2, 3)
    res$Group = group
    res$Habitat = habitat
    
    #################ADDITIONAL MANTEL TESTS###############################
    
    #Get Bray-Curtis dissimilarity matrix for community data
    commdist<-vegdist(comm, method="bray")
    
    #Get distance matrix for spatial factors
    geodist<-vegdist(geo, method="euclidean")
    
    #plot our data
    plot(geodist,commdist,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
    
    #Mantel test for correlation between community composition and spatial factors.
    comm_geo <- vegan::mantel(commdist, geodist, method="pearson", permutations=999)
    comm_geo
    
    #results dataframe
    res2 <- data.frame(Test = c("Environmental (Marginal)", "Geographic (Marginal)", "Environmental (Partial)", "Geographic (Partial)"), MantelR = c(0, round(comm_geo$statistic, 3), 0, 0), PVal = c(NA, round(comm_geo$signif, 3), NA, NA))
    res2$Group = group
    res2$Habitat = habitat
    
    return(list(res, res, res2))
    
  } else if (length(axis.keep2) < 1 & ncol(env.sel) < 1){
    
    #################ADDITIONAL MANTEL TESTS###############################
    
    #Get Bray-Curtis dissimilarity matrix for community data
    commdist<-vegdist(comm, method="bray")
    
    #Get distance matrix for spatial factors
    geodist<-vegdist(geo, method="euclidean")
    
    #Mantel test for correlation between community composition and spatial factors.
    comm_geo <- vegan::mantel(commdist, geodist, method="pearson", permutations=999)
    comm_geo.r2 = round(comm_geo$statistic, 2)
    comm_geo.p = comm_geo$signif
    
    #results dataframe
    res2 <- data.frame(Test = c("Environmental (Marginal)", "Geographic (Marginal)", "Environmental (Partial)", "Geographic (Partial)"), 
                       MantelR = c(NA, comm_geo.r2, NA, NA), 
                       PVal = c(NA, comm_geo.p, NA, NA))
    res2$Group = group
    res2$Habitat = habitat
    
    return(list(res2))
    
  }
  
}


#Prokaryotes

#Prokaryotes
pro.sn <- VPA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
pro.sn[[1]] #VPA
plot(pro.sn[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
pro.sn[[3]] #Mantel

pro.sp <- VPA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice")
pro.sp[[1]]
plot(pro.sp[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))

pro.sm <- VPA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice")
pro.sm[[1]]
plot(pro.sm[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))

pro.cr <- VPA_func(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite")
pro.cr[[1]]
plot(pro.cr[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))


#Snow

#plot
#for spliting axis labels over two lines
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

#define dataframe for venn diagram
df.venn <- data.frame(x = c(3, 1),y = c(1, 1),labels = c('A', 'B'))

##################Snow VENN DIAGRAM#########################
#pro snow
pro.sn[[1]]
pro.sn[[2]]
env = paste0("<0", "%") #less than zero so changing to zero
geo = paste0(round(pro.sn[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%")
combi = paste0(round(pro.sn[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(pro.sn[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

pro.sn.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE, aes(fill = labels)) +
  coord_fixed()+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  scale_color_manual(values=c("#88d2d6", "#f8c4be"))+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Snow")
pro.sn.plot

############SNOW BARPLOT####################
#snow mantel test barchart
df = pro.sn[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
pro.sn.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"), c(1, 0.7, 1, 0.7)))+
  ggtitle("Snow")
print(pro.sn.barplot)

#put the plots together
# pro.snow.plot.final = pro.sn.plot + pro.sn.barplot 
# # pro.snow.plot.final + ggtitle("Prokaryotes: Snow")
# pro.snow.plot.final = pro.snow.plot.final + plot_annotation(title = "Prokaryotes: Snow") & 
#   theme(plot.title = element_text(hjust = 0, size=10), plot.tag = element_text(size = 15, face = "bold"))
# print(pro.snow.plot.final)

#Spring Ice

#plot
##################SPRING ICE VENN DIAGRAM#########################
#pro spring
env = paste0(round(pro.sp[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%") #less than zero so changing to zero
geo = paste0(round(pro.sp[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%*")
combi = paste0(round(pro.sp[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(pro.sp[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

pro.sp.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Spring Ice")
pro.sp.plot

############SPRING ICE BARPLOT####################
#spring ice mantel test barchart
df = pro.sp[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")


#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
pro.sp.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Spring Ice")
print(pro.sp.barplot)

#put the plots together
# pro.sp.plot.final = pro.sp.plot + pro.sp.barplot 
# pro.sp.plot.final = pro.sp.plot.final + plot_annotation(title = "Prokaryotes: Spring Ice") & 
#   theme(plot.title = element_text(hjust = 0, size=10), plot.tag = element_text(size = 15, face = "bold"))
# print(pro.sp.plot.final)

# pro.snow.plot.final2 = pro.snow.plot.final + theme(axis.text.x = element_blank())
# pro.snow.plot.final2 / pro.sp.plot.final 
# 
# title1 <- grid::textGrob(label = "Prokaryotes: Snow")
# title2 <- grid::textGrob(label = "Prokaryotes: Spring Ice")
# 
# combine_1 =  (wrap_elements(panel = title1) / (pro.snow.plot.final)) + plot_layout(heights = c(1, 10))
# combine_2 =  (wrap_elements(panel = title2) / (pro.sp.plot.final)) + plot_layout(heights = c(1, 10))
# 
# (combine_1 / combine_2) +  plot_layout(heights = c(1, 10, 11))
# 
# plot_grid(pro.snow.plot.final2,pro.sp.plot.final, nrow=2, align="hv")


#Summer Ice

#plot
##################Summer ICE VENN DIAGRAM#########################
#pro Summer
pro.sm[[1]]
pro.sm[[2]]
env = paste0(round(pro.sm[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%*") #less than zero so changing to zero
geo = paste0(round(pro.sm[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%*")
combi = paste0(round(pro.sm[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(pro.sm[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

pro.sm.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Summer Ice")
pro.sm.plot

############Summer ICE BARPLOT####################
#Summer ice mantel test barchart
df = pro.sm[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
pro.sm.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Summer Ice")
print(pro.sm.barplot)

#put the plots together
# pro.sm.plot.final = pro.sm.plot + pro.sm.barplot 
# pro.sm.plot.final = pro.sm.plot.final + plot_annotation(title = "Prokaryotes: Summer Ice") & 
#   theme(plot.title = element_text(hjust = 0, size=10), plot.tag = element_text(size = 15, face = "bold"))
# print(pro.sm.plot.final)


#Cryoconite

#plot
##################cryoconite VENN DIAGRAM#########################
#pro Summer
pro.cr[[1]]
env = paste0(round(pro.cr[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%*")
geo = paste0(round(pro.cr[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%*")
combi = paste0(round(pro.cr[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(pro.cr[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

pro.cr.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Cryoconite")
pro.cr.plot

############cryoconite BARPLOT####################
#cryoconite mantel test barchart
df = pro.cr[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
pro.cr.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Cryoconite")
print(pro.cr.barplot)

#put the plots together
# pro.cr.plot.final = pro.cr.plot + pro.cr.barplot 
# pro.cr.plot.final = pro.cr.plot.final + plot_annotation(title = "Prokaryotes: Cryoconite") & 
#   theme(plot.title = element_text(hjust = 0, size=10), plot.tag = element_text(size = 10, face = "bold"),
#         plot.margin = unit(c(-7, 1, -7, 1), "cm"))
# print(pro.cr.plot.final)

# pro.cr.plot.final = pro.cr.plot + pro.cr.barplot 
# pro.cr.plot.final = pro.cr.plot.final + plot_annotation(title = "Prokaryotes: Cryoconite") & 
#   theme(plot.title = element_text(hjust = 0, size=10), plot.tag = element_text(size = 10, face = "bold"))
# print(pro.cr.plot.final)

#Combine

#VPA 
pro.p1 = pro.sn.plot + pro.sp.plot + pro.sm.plot + pro.cr.plot + pro.sn.barplot + pro.sp.barplot + pro.sm.barplot + pro.cr.barplot + plot_layout(ncol=2) + plot_annotation(tag_levels=list(c("A", "", "", "", "B", "", "", "")), 
                                                                                                                                                                           title = 'Prokaryotes')
pro.p1

# pro.p2 = pro.sn.plot + pro.sp.plot + pro.sm.plot + pro.cr.plot + pro.sn.barplot + pro.sp.barplot + pro.sm.barplot + pro.cr.barplot + plot_layout(ncol=2) + plot_annotation(title="Prokaryotes")
# pro.p2

pdf("../results-proportional/vpa-mantel-plots-prokaryote.pdf", width=5, height=8)
print(pro.p1)
dev.off()

# pdf("../results-proportional/vpa-mantel-plots-prokaryote2.pdf", width=5, height=8)
# print(pro.p2)
# dev.off()

#Microbial eukaryotes

#Eukaryotes
euk.sn <- VPA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Eukaryote", "Snow")
euk.sn[[1]]
plot(euk.sn[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
euk.sn[[3]]

euk.sp <- VPA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Eukaryote", "Spring Ice")
euk.sp[[1]]
euk.sp[[2]]
plot(euk.sp[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
euk.sp[[3]]

euk.sm <- VPA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Eukaryote", "Summer Ice")
euk.sm[[1]]
#plot(euk.sm[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
euk.sm[[3]]

euk.cr <- VPA_func(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Eukaryote", "Cryoconite")
euk.cr[[1]]
plot(euk.cr[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
euk.cr[[3]]


#Snow

#plot
##################Snow VENN DIAGRAM#########################
#euk snow
euk.sn[[1]]
env = paste0("<0", "%") #less than zero so changing to zero
geo = paste0("<0", "%") #less than zero so changing to zero
combi = paste0(round(euk.sn[[2]]$part$indfract$Adj.R.squared[[3]], 2)*100, "%")
resid = paste0(round(euk.sn[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

euk.sn.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Snow")
euk.sn.plot

############SNOW BARPLOT####################
#snow mantel test barchart
df = euk.sn[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
euk.sn.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Snow")
print(euk.sn.barplot)

#Spring Ice
#plot
##################SPRING ICE VENN DIAGRAM#########################
#euk spring
euk.sp[[1]]
env = paste0(round(euk.sp[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%")
geo = paste0(round(euk.sp[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%*")
combi = paste0(round(euk.sp[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(euk.sp[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

euk.sp.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Spring Ice")
euk.sp.plot

############SPRING ICE BARPLOT####################
#spring ice mantel test barchart
df = euk.sp[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")


#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
euk.sp.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Spring Ice")
print(euk.sp.barplot)


#Summer Ice
#plot
##################Summer ICE VENN DIAGRAM#########################
euk.sm[[1]]
env = paste0("<0", "%") #less than zero so changing to zero
geo = paste0(round(euk.sm[[1]]$AdjR2, 2)*100, "%*")
combi = paste0(round(0, 2)*100, "%")
resid = paste0(round(100 - (euk.sm[[1]]$AdjR2*100), 0), "%")

euk.sm.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Summer Ice")
euk.sm.plot

############Summer ICE BARPLOT####################
#Summer ice mantel test barchart
df = euk.sm[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
euk.sm.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Summer Ice")
print(euk.sm.barplot)

#Cryoconite
#plot
##################cryoconite VENN DIAGRAM#########################
euk.cr[[1]]
env = paste0("<0", "%") #less than zero so changing to zero
geo = paste0(round(euk.cr[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%*")
combi = paste0(round(euk.cr[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(euk.cr[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

euk.cr.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Cryoconite")
euk.cr.plot

############cryoconite BARPLOT####################
#cryoconite mantel test barchart
df = euk.cr[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
euk.cr.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Cryoconite")
print(euk.cr.barplot)

#Combine

#VPA 
euk.p1 = euk.sn.plot + euk.sp.plot + euk.sm.plot + euk.cr.plot + euk.sn.barplot + euk.sp.barplot + euk.sm.barplot + euk.cr.barplot + plot_layout(ncol=2) + plot_annotation(tag_levels=list(c("A", "", "", "", "B", "", "", "")), 
                                                                                                                                                                           title = 'Eukaryotes')
euk.p1

pdf("../results-proportional/vpa-mantel-plots-eukaryote.pdf", width=5, height=8)
print(euk.p1)
dev.off()

#Microfauna
#mmkaryotes
mm.sn <- VPA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Microfauna", "Snow")
mm.sn[[1]]

mm.sp <- VPA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Microfauna", "Spring Ice")
mm.sp[[1]]
plot(mm.sp[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
mm.sp[[3]]

mm.sm <- VPA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Microfauna", "Summer Ice")
mm.sm[[1]]

# plot(mm.sm[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
mm.cr <- VPA_func(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Microfauna", "Cryoconite")
mm.cr[[1]]
plot(mm.cr[[2]], digits=2, Xnames=c('Environmental', 'Spatial'), id.size=0.75, bg = c('navy', 'tomato'))
mm.cr[[3]]

#Snow
##################Snow VENN DIAGRAM#########################
#mm snow
mm.sn[[1]]
env = paste0(0, "%") #less than zero so changing to zero
geo = paste0(0, "%")
combi = paste0(0, "%")
resid = paste0(100, "%")

mm.sn.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Snow")
mm.sn.plot

############SNOW BARPLOT####################
#snow mantel test barchart
df = mm.sn[[1]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
mm.sn.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Snow")
print(mm.sn.barplot)

#Spring Ice
#plot
##################SPRING ICE VENN DIAGRAM#########################
#mm spring
mm.sp[[1]]
env = paste0(round(mm.sp[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%*")
geo = paste0(round(mm.sp[[2]]$part$indfract$Adj.R.squared[2], 2)*100, "%")
combi = paste0(round(mm.sp[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(mm.sp[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

mm.sp.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Spring Ice")
mm.sp.plot

############SPRING ICE BARPLOT####################
#spring ice mantel test barchart
df = mm.sp[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")


#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
mm.sp.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Spring Ice")
print(mm.sp.barplot)

#Summer Ice
#plot
##################Summer ICE VENN DIAGRAM#########################
#mm Summer
mm.sm[[1]]
env = paste0(0, "%")
geo = paste0(0, "%")
combi = paste0(0, "%")
resid = paste0(100, "%")

mm.sm.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Summer Ice")
mm.sm.plot

############Summer ICE BARPLOT####################
#Summer ice mantel test barchart
df = mm.sm[[1]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")
df$PVal = NA

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
mm.sm.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.title.y = element_text(size=7), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Summer Ice")
print(mm.sm.barplot)

#Cryoconite
#plot
##################cryoconite VENN DIAGRAM#########################
#mm Summer
mm.cr[[1]]
env = paste0(round(mm.cr[[2]]$part$indfract$Adj.R.squared[1], 2)*100, "%") 
geo = paste0("<0", "%")
combi = paste0(round(mm.cr[[2]]$part$indfract$Adj.R.squared[3], 2)*100, "%")
resid = paste0(round(mm.cr[[2]]$part$indfract$Adj.R.squared[[4]], 2)*100, "%")

mm.cr.plot <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 2, colour = 'grey', show.legend = FALSE ) +
  coord_fixed()+
  scale_fill_manual(values=c("#88d2d6", "#f8c4be"))+
  annotate("text", x = c(3, 1, 2, 3, 0, 4), y = c(1, 1, 1, -0.5, 2.5, 2.5), 
           label = c(geo, env, combi, paste('Residual :', resid), "Env", "Geo"), 
           size = 2)+
  theme_bw()+
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.title = element_text(hjust = 0.5, size=8), plot.tag = element_text(size = 8))+
  ggtitle("Cryoconite")
mm.cr.plot

############cryoconite BARPLOT####################
#cryoconite mantel test barchart
df = mm.cr[[3]]
df$Test = c("Env (Marginal)", "Geo (Marginal)", "Env (Partial)", "Geo (Partial)")

#add significance indicator 
df$significant <- dplyr::case_when(
  df$PVal < 0.05 & df$MantelR > 0 ~ TRUE,
  TRUE ~ FALSE
)

#df$MantelR[df$MantelR < 0] = 0
mm.cr.barplot = ggplot(df, aes(x = Test, y = MantelR, fill=Test))+
  geom_bar(stat="identity")+
  ylab("Correlation Coefficient")+
  geom_text(aes(label = ifelse(significant, "*", ""), group = Test), 
            position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt)+
  ylim(-0.2, 0.85)+
  scale_x_discrete(breaks=unique(df$Test), labels=addline_format(c(df$Test)))+
  theme_bw()+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=8), axis.text.x = element_text(size=5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.tag = element_text(size = 8))+
  scale_fill_manual(values=ggplot2::alpha(c("#f8c4be", "#f8c4be", "#88d2d6", "#88d2d6"),c(1, 0.7, 1, 0.7)))+
  ggtitle("Cryoconite")
print(mm.cr.barplot)

#Combine
#VPA 
mm.p1 = mm.sn.plot + mm.sp.plot + mm.sm.plot + mm.cr.plot + mm.sn.barplot + mm.sp.barplot + mm.sm.barplot + mm.cr.barplot + plot_layout(ncol=2) + plot_annotation(tag_levels=list(c("A", "", "", "", "B", "", "", "")), 
                                                                                                                                                                  title = 'Microfauna')
mm.p1

pdf("../results-proportional/vpa-mantel-plots-microfauna.pdf", width=5, height=8)
print(mm.p1)
dev.off()


#all three plots together
all.p = pro.p1 | euk.p1  | mm.p1 

# Set theme for annotations
thm <- theme(plot.title = element_text(face = 2, size = 10))
first_plot      <- wrap_elements(pro.p1 + plot_annotation(title = "Prokaryote", theme = thm))
second_plot   <- wrap_elements(euk.p1 + plot_annotation(title = "Microbial Eukaryote", theme = thm))
third_plot   <- wrap_elements(mm.p1 + plot_annotation(title = "Microfauna", theme = thm))

all.p = first_plot | second_plot | third_plot


pdf("../results-proportional/vpa-mantel-plots-all.pdf", width=13, height=8)
print(all.p)
dev.off()

#Table of results
#for mantel test we report 
man.res = rbind(pro.sn[[3]], pro.sp[[3]], pro.sm[[3]], pro.cr[[3]],
                euk.sn[[3]], euk.sp[[3]], euk.sm[[3]], euk.cr[[3]])
names(man.res) = c("Test", "Mantel_AdjR2", "Mantel_PVal", "Group", "Habitat")

write.csv(man.res, "../results-proportional/mantel-partial-results.csv")

vpa.res = rbind(pro.sn[[1]], pro.sp[[1]], pro.sm[[1]], pro.cr[[1]],
                euk.sn[[1]], euk.sp[[1]], euk.sm[[1]], euk.cr[[1]])
names(vpa.res) = c("Test", "VPA_AdjR2", "VPA_PVal", "Group", "Habitat")

write.csv(vpa.res, "../results-proportional/variance-partitioning-results.csv")

#micrometazoa results
man.res.mm = rbind(mm.sn[[1]], mm.sp[[3]], mm.sm[[1]], mm.cr[[3]])

names(man.res.mm) = c("Test", "Mantel_AdjR2", "Mantel_PVal", "Group", "Habitat")

write.csv(man.res.mm, "../results-proportional/mantel-partial-results-Microfauna.csv")

vpa.res.mm = rbind(mm.sp[[1]], mm.cr[[1]])

write.csv(vpa.res.mm, "../results-proportional/variance-partitioning-results-microfauna.csv")

#