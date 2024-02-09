# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

# library(ggpubr)
#install.packages("ggstatsplot")
# library(ggstatsplot)
source("00-solman-functions.R")
# library(ggcorrplot)
#library(fdrtool)
# library(psych)
# library(cowplot)
# #devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)
library(stringr) #for subsetting strings in lat long function
library(fossil) #for earth.dist function
library(vegan) #for community dissimilarity matrix
library(tibble) #for rownames to column
library(dplyr) #coalesce function + rows update function
library(reshape2) #for melt function
library(funrar)
library(phyloseq)
library(ggplot2)
#install.packages("ggpmisc")
library(ggpmisc)
library(picante) #for faith's pd
library(cowplot)

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

#Arctic prokaryotes
pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds") 
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds") 

#Antarctic prokaryotes
pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun.rds") 
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int.rds") 
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare.rds") 

#Arctic eukaryotes
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds") 
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds") 

#Antarctic eukaryotes
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun.rds") 
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int.rds") 
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare.rds") 

# Cryoconite hole area from phyloseq object

# hole_func <- function(phylo){
#   
#   #put our data into a data frame
# #Get cryoconite areas
# samp_data <- data.frame(sample_data(phylo))
#   
# #get cryoconite hole areas
# area1 <- pi*samp_data$Radius^2
# area2 <- pi*(samp_data$NS/2)*(samp_data$EW/2)
# Areas <- as.data.frame(coalesce(area1,area2))
# row.names(Areas) = row.names(samp_data)
# names(Areas) = "Area"
# 
# return(Areas)
# 
# }


#Get ASV richness and Faith's PD from a phyloseq object
# phylo = pro
# plot_richness(pro, measures = c("Observed"))
how.rich <- function(phylo){
  
  comm = data.frame(t(otu_table(phylo)), check.names=FALSE)
  tree = phy_tree(phylo)
  alpha.df = pd(comm, tree, include.root = TRUE)
  
  return(alpha.df)
}

#Bind data for analysis
# phylo = pro.ant.abun
# subcom = "Abundant"
# group = "Prokaryote"

area_and_rich <- function(phylo, subcom, group){
  
  # area = hole_func(phylo)
  area = data.frame(sample_data(phylo))$Area
  rich = how.rich(phylo)
  
  #bind into data frame
  data = data.frame(cbind(area, rich$PD, rich$SR))
  data$Subcommunity = subcom
  data$Group = group
  
  #remove rows with na
  data.comp = data.frame(data[complete.cases(data),])
  
  names(data.comp) = c("area", "PD", "SR", "Subcommunity", "Group")
  
  return(data.comp)
}

#get data for analysis
pro_abun_data = area_and_rich(pro.ant.abun, "Abundant", "Prokaryote")
pro_int_data = area_and_rich(pro.ant.int, "Intermediate", "Prokaryote")
pro_rare_data = area_and_rich(pro.ant.rare, "Rare", "Prokaryote")
euk_abun_data = area_and_rich(euk.ant.abun, "Abundant", "Eukaryote")
euk_int_data = area_and_rich(euk.ant.int, "Intermediate", "Eukaryote")
euk_rare_data = area_and_rich(euk.ant.rare, "Rare", "Eukaryote")

total.data = rbind(pro_abun_data, pro_int_data, pro_rare_data, euk_abun_data, euk_int_data, euk_rare_data)

# df = total.data
# subcom = "Abundant"
# group = "Prokaryote"
# alpha = "SR"

lm_mod_func <- function(df, subcom, group, alpha){
  
  test.data = df[which(df$Subcommunity == subcom & df$Group == group),]
  test.data = test.data[test.data$SR > 0,]
  
  if (alpha == "PD"){
    
    #remove outliers
    # Q <- quantile(test.data$PD, probs=c(.25, .75), na.rm = FALSE)
    # iqr <- IQR(test.data$PD)
    # #find the cutoff ranges beyond which data points are outliers
    # up <-  Q[2]+1.5*iqr # Upper Range
    # low<- Q[1]-1.5*iqr # Lower Range
    # eliminated<- subset(test.data, test.data$PD > (Q[1] - 1.5*iqr) & test.data$PD < (Q[2]+1.5*iqr))
    
    #log transform our data before fitting the linear model
    test.data$area <- log10(test.data$area)
    test.data$PD <- log10(test.data$PD)
    
    mod = lm(PD ~ area, data = test.data)
    mod
    x = summary(mod)
    
    
  } else if (alpha == "SR"){
    
    #remove outliers
    # Q <- quantile(test.data$SR, probs=c(.25, .75), na.rm = FALSE)
    # iqr <- IQR(test.data$SR)
    # #find the cutoff ranges beyond which data points are outliers
    # up <-  Q[2]+1.5*iqr # Upper Range
    # low<- Q[1]-1.5*iqr # Lower Range
    # eliminated<- subset(test.data, test.data$SR > (Q[1] - 1.5*iqr) & test.data$SR < (Q[2]+1.5*iqr))
    
    #log transform our data before fitting the linear model
    test.data$area <- log10(test.data$area)
    test.data$SR <- log10(test.data$SR)
    
    mod = lm(SR ~ area, data = test.data)
    mod
    x = summary(mod)
  }
  
  
  #Linear regression results into a data frame
  res.df = data.frame(Group = group, Subcommunity = subcom, Metric = alpha, Slope = round(mod$coefficients[[2]], 3), 
                      R2 = round(x$r.squared, 3), 
                      PVal = round(lmp(mod), 4))
  
  return(list(res.df, test.data))
}

pro_abun_lm_mod_PD = lm_mod_func(total.data, "Abundant", "Prokaryote", "PD")
pro_abun_lm_mod_SR = lm_mod_func(total.data, "Abundant", "Prokaryote", "SR")
pro_int_lm_mod_PD = lm_mod_func(total.data, "Intermediate", "Prokaryote", "PD")
pro_int_lm_mod_SR = lm_mod_func(total.data, "Intermediate", "Prokaryote", "SR")
pro_rare_lm_mod_PD = lm_mod_func(total.data, "Rare", "Prokaryote", "PD")
pro_rare_lm_mod_SR = lm_mod_func(total.data, "Rare", "Prokaryote", "SR")
euk_abun_lm_mod_PD = lm_mod_func(total.data, "Abundant", "Eukaryote", "PD")
euk_abun_lm_mod_SR = lm_mod_func(total.data, "Abundant", "Eukaryote", "SR")
euk_int_lm_mod_PD = lm_mod_func(total.data, "Intermediate", "Eukaryote", "PD")
euk_int_lm_mod_SR = lm_mod_func(total.data, "Intermediate", "Eukaryote", "SR")
euk_rare_lm_mod_PD = lm_mod_func(total.data, "Rare", "Eukaryote", "PD")
euk_rare_lm_mod_SR = lm_mod_func(total.data, "Rare", "Eukaryote", "SR")

#bind results
lm_mod.res = rbind(pro_abun_lm_mod_PD[[1]], pro_abun_lm_mod_SR[[1]], 
                   pro_int_lm_mod_PD[[1]], pro_int_lm_mod_SR[[1]], 
                   pro_rare_lm_mod_PD[[1]], pro_rare_lm_mod_SR[[1]], 
                   euk_abun_lm_mod_PD[[1]], euk_abun_lm_mod_SR[[1]], 
                   euk_int_lm_mod_PD[[1]], euk_int_lm_mod_SR[[1]], 
                   euk_rare_lm_mod_PD[[1]], euk_rare_lm_mod_SR[[1]])

#save results
write.csv(lm_mod.res, "../results/TAR-results.csv")

total.data.final = rbind(pro_abun_lm_mod_SR[[2]], 
                         pro_int_lm_mod_SR[[2]], 
                         pro_rare_lm_mod_SR[[2]], 
                         euk_abun_lm_mod_SR[[2]], 
                         euk_int_lm_mod_SR[[2]], 
                         euk_rare_lm_mod_SR[[2]])

total.data.final = total.data.final %>%
  mutate(across(Group, factor, levels=c("Prokaryote","Eukaryote")))%>%
  mutate(across(Subcommunity, factor, levels=c("Rare","Intermediate", "Abundant")))

p1 = ggplot(total.data.final, aes(x=area, y=SR, fill=Subcommunity, shape=Subcommunity)) +
  geom_point(data=total.data.final, aes(x=area, y=SR, fill=Subcommunity),
             alpha=.5, size=3) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, show.legend = FALSE, aes(color=Subcommunity))+
  facet_wrap(~Group, ncol=1, scale="free")+
  theme_bw()+
  # xlab("log(area cm2)")+
  # ylab("log(ASV richness)")+
  # labs(x = bquote('log(area cm',^2,')', y = "y axis")
  labs(x = expression(paste("log(area cm"^2,")")), y = "log(ASV richness)")+
  scale_shape_manual(values = c(21, 21, 21))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        legend.text = element_text(size=15), 
        legend.title = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 7, alpha=.8) ) )

p1

p2 = ggplot(total.data.final, aes(x=area, y=PD, fill=Subcommunity, shape=Subcommunity)) +
  geom_point(data=total.data.final, aes(x=area, y=PD, fill=Subcommunity),
             alpha=.5, size=3) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, show.legend = FALSE, aes(color=Subcommunity))+
  facet_wrap(~Group, ncol=1, scale="free")+
  theme_bw()+
  # xlab("log(area cm2)")+
  # ylab("log(ASV richness)")+
  # labs(x = bquote('log(area cm',^2,')', y = "y axis")
  labs(x = expression(paste("log(area cm"^2,")")), y = "log(Faith's PD)")+
  scale_shape_manual(values = c(21, 21, 21))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  #scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        legend.text = element_text(size=15), 
        legend.title = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 7, alpha=.8) ) )

p2

pdf("../results/TAR-ASV-richness.pdf")
print(p1)
dev.off()

pdf("../results/TAR-Faiths-PD.pdf")
print(p2)
dev.off()
