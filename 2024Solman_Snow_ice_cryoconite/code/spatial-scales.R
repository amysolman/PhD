# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

source("00-solman-functions.R")
library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(cowplot) #plot_grid
library(reshape2) #for function melt
library(ggpmisc) #for adding lm stats to the plot
library(broom) #for adding stats to the plot
library(fossil) #for earth.dist function
library(patchwork)
library(ggpubr)
library(dplyr)
library(plyr)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#Looking at different spatial scales

spatial_scales_data <- function(ps, sam, hab){
  
  #Separate habitat data
  h <- prune_samples(sam, ps)
  
  #remove empty samples or ASVs
  h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
  h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)
  
  #Calculate Bray-Curtis dissimilarities between sites
  bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")
  
  #Calculate geographic distance between sites
  ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
  geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers
  
  #Calculate environmental dissimilarity between sites
  chem <- data.frame(sample_data(h.f))
  chem <- chem[,c(9,10,13:ncol(chem))]
  #z-score transform the environmental data
  chem.z = data.frame(scale(chem))
  #Get get distance matrix for environmental factors
  chem.dist<-vegdist(chem.z, method="euclidean", na.rm=TRUE)
  
  #transform these data into columns
  comm.dist = melt(as.matrix(bc))
  names(comm.dist) = c("Sample1", "Sample2", "Beta")
  comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]
  
  geo.dist = melt(as.matrix(geo.dist))
  names(geo.dist) = c("Sample1", "Sample2", "km")
  geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]
  
  chem.dist = melt(as.matrix(chem.dist))
  names(chem.dist) = c("Sample1", "Sample2", "ChemDist")
  chem.dist.trim = chem.dist[chem.dist['Sample1'] != chem.dist['Sample2'],]
  
  #make into dataframe
  data = comm.dist.trim 
  data$km = geo.dist.trim$km
  data$m = data$km*1000
  data$ChemDist = chem.dist.trim$ChemDist
  
  #spearmans correlations tests
  cor(data$m, data$Beta, method="spearman") #spearman's correlation
  cor(data$m, data$ChemDist, method="spearman") #spearman's correlation
  cor(data$ChemDist, data$Beta, method="spearman") #spearman's correlation
  
  df_55 = data[data$m <= 55,]
  df_55$Group = "<55m"
  df_200 = data[data$m > 55 & data$m <= 200,]
  df_200$Group = "55-200m"
  df_300 = data[data$m > 200 & data$m <= 300,]
  df_300$Group = "200-300m"
  df_400 = data[data$m > 300 & data$m <= 400,]
  df_400$Group = "300-400m"
  df_500 = data[data$m > 400 & data$m <= 500,]
  df_500$Group = "400-500m"
  df_750 = data[data$m > 500 & data$m <= 750,]
  df_750$Group = "500-750m"
  df_1000 = data[data$m > 750 & data$m <= 1000,]
  df_1000$Group = "750-1000m"
  df_1250 = data[data$m > 1000 & data$m <= 1250,]
  df_1250$Group = "1000-1250m"
  df_1500 = data[data$m > 1250 & data$m <= 1500,]
  df_1500$Group = "1250-1500m"
  df_2200 = data[data$m > 1500 & data$m <= 2200,]
  df_2200$Group = "1500-2200m"
  
  #DATA FOR TESTING AND PLOTTING
  test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
  test.data$Habitat = hab
  
  return(test.data)
  
}


####FUNCTION FOR PLOTTING THE DATA

#df = pro.sum.stats

scales_plot_function <- function(df){
  
  
  ####A DIFFERENT PLOT
  snow.df = df[df$Habitat == "Snow",]
  
  scale=20
  p1 = ggplot(snow.df, aes(x=Group, y=Beta_mean)) +
    geom_line(aes(group = 1), color="#fa9f99") +
    geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
    scale_y_continuous(limits=c(0,1)) +
    expand_limits(y = 0)+
    theme_bw()+
    geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#fa9f99") +
    geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
    geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#fa9f99")+
    geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x=element_blank())+
    ylab("Community Bray-Curtis Dissimilarity")+
    ggtitle("Snow")
  p1
  
  
  sp.df = df[df$Habitat == "Spring Ice",]
  
  p2 = ggplot(sp.df, aes(x=Group, y=Beta_mean)) +
    geom_line(aes(group = 1), color="#a4c64d") +
    geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
    scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
    expand_limits(y = 0)+
    theme_bw()+
    geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#a4c64d") +
    geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
    geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#a4c64d")+
    geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank())+
    ylab("")+
    ggtitle("Spring Ice")
  p2
  
  
  sum.df = df[df$Habitat == "Summer Ice",]
  p3 = ggplot(sum.df, aes(x=Group, y=Beta_mean)) +
    geom_line(aes(group = 1), color="#4dd2d6") +
    geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
    scale_y_continuous(limits=c(0,1)) +
    expand_limits(y = 0)+
    theme_bw()+
    geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#4dd2d6") +
    geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
    geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#4dd2d6")+
    geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1, size=10))+
    ylab("Community Bray-Curtis Dissimilarity")+
    ggtitle("Summer Ice")
  p3
  
  cryo.df = df[df$Habitat == "Cryoconite",]
  p4 = ggplot(cryo.df, aes(x=Group, y=Beta_mean)) +
    geom_line(aes(group = 1), color="#d8a4ff") +
    geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
    scale_y_continuous(limits=c(0,1), 
                       sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
    expand_limits(y = 0)+
    theme_bw()+
    geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean), width = 0.2, color="#d8a4ff") +
    geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
    geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#d8a4ff")+
    geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
    theme(legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust=1, size=10),
          axis.text.y.left = element_blank(),
          axis.ticks.y.left = element_blank())+
    ylab("")+
    ggtitle("Cryoconite")
  p4
  
  line.p = p1 + p2 + p3 + p4
  line.p
  
  
  return(line.p)
}

###PROKARYOTE

#get test data
pro.snow.df = spatial_scales_data(ps = ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, hab = "Snow")
pro.sp.df = spatial_scales_data(ps = ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, hab = "Spring Ice")
pro.sm.df = spatial_scales_data(ps = ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, hab = "Summer Ice")
pro.cr.df = spatial_scales_data(ps = ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, hab = "Cryoconite")

#put all our dataframes together
pro.all.habs = rbind(pro.snow.df, pro.sp.df, pro.sm.df, pro.cr.df)

#STEP ONE ORDER THE LEVELS TO BE TESTED
pro.all.habs$Group <- ordered(pro.all.habs$Group,
                              levels = c("<55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                         "1000-1250m", "1250-1500m", "1500-2200m"))

pro.all.habs$Habitat <- ordered(pro.all.habs$Habitat,
                                levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))


# STEP TWO COMPUTE SUMMARY STATISTICS
pro.sum.stats = group_by(pro.all.habs, Group, Habitat) %>%
  dplyr::summarise(
    count = n(),
    Beta_mean = mean(Beta, na.rm = TRUE),
    Beta_sd = sd(Beta, na.rm = TRUE),
    Chem_mean = mean(ChemDist, na.rm = TRUE),
    Chem_sd = sd(ChemDist, na.rm = TRUE))

pro.line.p = scales_plot_function(pro.sum.stats)
pro.line.p

pdf("../results/prokaryote-line-distance.pdf", width=20, height=10)
print(pro.line.p)
dev.off()

#####TEST FOR DIFFERENCES
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?

#snow
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.snow.df[pro.snow.df$m < 55,]$ChemDist
y = pro.snow.df[pro.snow.df$m > 1500,]$ChemDist

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#conduct t.test (which will be a welch test which assumes normality but not equal variances)
snow.euc = t.test(x,y,alternative = "two.sided", var.equal=FALSE) #two-sided = greater or less, var.equal = are the variances equal? If FALSE then a Welch test is used.
#answer = no!

#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.snow.df[pro.snow.df$m < 55,]$Beta
y = pro.snow.df[pro.snow.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p > 0.05 so variances are not significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
pro.snow.com = wilcox.test(x, y, alternative = "two.sided") # p > 0.05 so not significantly different


#spring ice
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.sp.df[pro.sp.df$m < 55,]$ChemDist
y = pro.sp.df[pro.sp.df$m > 1500,]$ChemDist

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
sp.euc = wilcox.test(x, y, alternative = "two.sided") # p < 0.05  significantly different

#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.sp.df[pro.sp.df$m < 55,]$Beta
y = pro.sp.df[pro.sp.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
pro.sp.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different



#summer ice
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.sm.df[pro.sm.df$m < 55,]$ChemDist
y = pro.sm.df[pro.sm.df$m > 1500,]$ChemDist

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
sm.euc = wilcox.test(x, y, alternative = "two.sided") # p < 0.05  significantly different

#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.sm.df[pro.sm.df$m < 55,]$Beta
y = pro.sm.df[pro.sm.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
pro.sm.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different


#cryoconite
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.cr.df[pro.cr.df$m < 55,]$ChemDist
y = pro.cr.df[pro.cr.df$m > 1500,]$ChemDist

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
cr.euc = wilcox.test(x, y, alternative = "two.sided") # p < 0.05  significantly different

#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = pro.cr.df[pro.cr.df$m < 55,]$Beta
y = pro.cr.df[pro.cr.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
pro.cr.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different

#EUKARYOTE

#get test data
euk.snow.df = spatial_scales_data(ps = ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, hab = "Snow")
euk.sp.df = spatial_scales_data(ps = ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, hab = "Spring Ice")
euk.sm.df = spatial_scales_data(ps = ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, hab = "Summer Ice")
euk.cr.df = spatial_scales_data(ps = ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, hab = "Cryoconite")

#put all our dataframes together
euk.all.habs = rbind(euk.snow.df, euk.sp.df, euk.sm.df, euk.cr.df)

#STEP ONE ORDER THE LEVELS TO BE TESTED
euk.all.habs$Group <- ordered(euk.all.habs$Group,
                              levels = c("<55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                         "1000-1250m", "1250-1500m", "1500-2200m"))

euk.all.habs$Habitat <- ordered(euk.all.habs$Habitat,
                                levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))


# STEP TWO COMPUTE SUMMARY STATISTICS
euk.sum.stats = dplyr::group_by(euk.all.habs, Group, Habitat) %>%
  dplyr::summarise(
    count = n(),
    Beta_mean = mean(Beta, na.rm = TRUE),
    Beta_sd = sd(Beta, na.rm = TRUE),
    Chem_mean = mean(ChemDist, na.rm = TRUE),
    Chem_sd = sd(ChemDist, na.rm = TRUE))

euk.line.p = scales_plot_function(euk.sum.stats)
euk.line.p

pdf("../results/eukaryote-line-distance.pdf", width=20, height=10)
print(euk.line.p)
dev.off()

#####TEST FOR DIFFERENCES
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?

#snow
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = euk.snow.df[euk.snow.df$m < 55,]$Beta
y = euk.snow.df[euk.snow.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p > 0.05 so variances are not significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
euk.snow.com = wilcox.test(x, y, alternative = "two.sided") # p > 0.05 so not significantly different


#spring ice
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = euk.sp.df[euk.sp.df$m < 55,]$Beta
y = euk.sp.df[euk.sp.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#perform Welch's test for normally distributed data with unequal variances
euk.sp.com = t.test(x,y,alternative = "two.sided", var.equal=FALSE)


#summer ice
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = euk.sm.df[euk.sm.df$m < 55,]$Beta
y = euk.sm.df[euk.sm.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
euk.sm.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different


#cryoconite
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = euk.cr.df[euk.cr.df$m < 55,]$Beta
y = euk.cr.df[euk.cr.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p > 0.05 so variances are not significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
euk.cr.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different

#MICROFAUNA

#get test data
mm.snow.df = spatial_scales_data(ps = ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, hab = "Snow")
mm.sp.df = spatial_scales_data(ps = ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, hab = "Spring Ice")
mm.sm.df = spatial_scales_data(ps = ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, hab = "Summer Ice")
mm.cr.df = spatial_scales_data(ps = ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, hab = "Cryoconite")

#put all our dataframes together
mm.all.habs = rbind(mm.snow.df, mm.sp.df, mm.sm.df, mm.cr.df)

#STEP ONE ORDER THE LEVELS TO BE TESTED
mm.all.habs$Group <- ordered(mm.all.habs$Group,
                              levels = c("<55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                         "1000-1250m", "1250-1500m", "1500-2200m"))

mm.all.habs$Habitat <- ordered(mm.all.habs$Habitat,
                                levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))


# STEP TWO COMPUTE SUMMARY STATISTICS
mm.sum.stats = dplyr::group_by(mm.all.habs, Group, Habitat) %>%
  dplyr::summarise(
    count = n(),
    Beta_mean = mean(Beta, na.rm = TRUE),
    Beta_sd = sd(Beta, na.rm = TRUE),
    Chem_mean = mean(ChemDist, na.rm = TRUE),
    Chem_sd = sd(ChemDist, na.rm = TRUE))

mm.line.p = scales_plot_function(mm.sum.stats)
mm.line.p

pdf("../results/microfauna-line-distance.pdf", width=20, height=10)
print(mm.line.p)
dev.off()

#####TEST FOR DIFFERENCES
#Is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?

#snow
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = mm.snow.df[mm.snow.df$m < 55,]$Beta
y = mm.snow.df[mm.snow.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p > 0.05 so variances are not significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
mm.snow.com = wilcox.test(x, y, alternative = "two.sided") # p > 0.05 so not significantly different


#spring ice
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = mm.sp.df[mm.sp.df$m < 55,]$Beta
y = mm.sp.df[mm.sp.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. not significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. not significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
mm.sp.com = wilcox.test(x, y, alternative = "two.sided")


#summer ice
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = mm.sm.df[mm.sm.df$m < 55,]$Beta
y = mm.sm.df[mm.sm.df$m > 1300,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p > 0.05 so the data is normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
mm.sm.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different


#cryoconite
#Is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?
x = mm.cr.df[mm.cr.df$m < 55,]$Beta
y = mm.cr.df[mm.cr.df$m > 1500,]$Beta

#check t.test assumptions
# Shapiro-Wilk normality test for <55m distance
shapiro.test(x) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
shapiro.test(y) #p < 0.05 so the data is not normally distributed i.e. significant different from normal distribution
#Do they have the same variances
var.test(x, y) #p < 0.05 so variances are significantly different
#we will use a Wilcoxon test instead (non-parametric alternative to t.test)
mm.cr.com = wilcox.test(x, y, alternative = "two.sided") # p < 0.05 so significantly different

###################################################################################################
###################################################################################################

#Report RESULTS
sink("../results/distance-decay-relationship-results.txt", type="output")
writeLines("===============================================================
DISTANCE DECAY RESULTS
===============================================================")
writeLines("For snow, is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?")
snow.euc
writeLines("For prokayrotes in snow, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
pro.snow.com
writeLines("For eukayrotes in snow, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
euk.snow.com
writeLines("For microfauna in snow, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
mm.snow.com

writeLines("For spring ice, is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?")
sp.euc
writeLines("For prokayrotes in spring ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
pro.sp.com
writeLines("For eukayrotes in spring ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
euk.sp.com
writeLines("For microfauna in spring ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
mm.sp.com

writeLines("For summer ice, is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?")
sp.euc
writeLines("For prokayrotes in summer ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
pro.sm.com
writeLines("For eukayrotes in summer ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
euk.sm.com
writeLines("For microfauna in summer ice, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
mm.sm.com

writeLines("For cryoconite, is the mean chemical dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?")
cr.euc
writeLines("For prokayrotes in cryoconite, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
pro.cr.com
writeLines("For eukayrotes in cryoconite, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
euk.cr.com
writeLines("For microfauna in cryoconite, is the mean community dissimilarity significantly higher between samples 1500-2200m apart compared to <55m apart?") 
mm.cr.com

sink()

###################################################################################################
###################################################################################################

#Combine Plots + save results tables
write.csv(pro.sum.stats, "../results/prokaryote-distance-stats.csv")
write.csv(euk.sum.stats, "../results/eukaryote-distance-stats.csv")
write.csv(mm.sum.stats, "../results/microfauna-distance-stats.csv")

#plot
thm <- theme(plot.title = element_text(face = "bold", size = 30))
first_plot      <- wrap_elements(pro.line.p + plot_annotation(title = "Prokaryote", theme = thm))
second_plot   <- wrap_elements(euk.line.p + plot_annotation(title = "Microbial Eukaryote", theme = thm))
third_plot   <- wrap_elements(mm.line.p + plot_annotation(title = "Microfauna", theme = thm))

all.p = first_plot / second_plot / third_plot
all.p

pdf("../results/environmental-community-dissimilarities.pdf", width=12, height=25)
print(all.p)
dev.off()



#Combined plot

#put all data together
pro.sum.stats$Taxa = "Prokaryote"
euk.sum.stats$Taxa = "Microbial Eukaryote"
mm.sum.stats$Taxa = "Microfauna"
data2plot = rbind(pro.sum.stats, euk.sum.stats, mm.sum.stats)

#plot

#Snow
snow.df = data2plot[data2plot$Habitat == "Snow",]
chem_data = snow.df[snow.df$Taxa == "Prokaryote",]
#scale=12

#scale the data
#Taken from https://stackoverflow.com/questions/58774705/y-limits-for-ggplot-with-sec-axis
ylim.prim <- c(0, 1)
ylim.sec <- c(6, 12)
EUCL <- snow.df$Chem_mean #needed for coherent normalisation
# This is quite hacky, but it works if you want to set a boundary for the secondary y-axis
fit = lm(b ~ . + 0, 
         tibble::tribble(
           ~a, ~s,  ~b,
           1,  (ylim.sec[1] - mean(EUCL))/sd(EUCL),  ylim.prim[1],
           1,  (ylim.sec[2] - mean(EUCL))/sd(EUCL), ylim.prim[2]))
a <- fit$coefficients['a']
s <- fit$coefficients['s']

snow.df$Taxa <- factor(snow.df$Taxa, levels=c("Prokaryote", "Microbial Eukaryote", "Microfauna"))
snow.p <-ggplot(snow.df, aes(x=Group, y=Beta_mean, group=Taxa)) +
  geom_line(aes(color="#fa9f99"))+
  geom_point(aes(color="#fa9f99", shape=Taxa, size=4))+
  geom_point(data = chem_data, aes(x = Group, y =(a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s), size=4), color="#000000")+
  geom_line(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s)), color="#000000")+
  scale_shape_manual(values=c(10,9,8))+
  scale_y_continuous(limits=c(0,1)) +
  expand_limits(y = 0)+
  theme_bw()+
  # geom_errorbar(aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Taxa), width = 0.2) +
  # geom_errorbar(data = chem_data, aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000")
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        plot.title = element_text(face="bold"))+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Snow")
  
snow.p
  

#Spring Ice
sp.df = data2plot[data2plot$Habitat == "Spring Ice",]
chem_data = sp.df[sp.df$Taxa == "Prokaryote",]

#scale the data
#Taken from https://stackoverflow.com/questions/58774705/y-limits-for-ggplot-with-sec-axis
ylim.prim <- c(0, 1)
ylim.sec <- c(6, 12)
EUCL <- sp.df$Chem_mean #needed for coherent normalisation
# This is quite hacky, but it works if you want to set a boundary for the secondary y-axis
fit = lm(b ~ . + 0, 
         tibble::tribble(
           ~a, ~s,  ~b,
           1,  (ylim.sec[1] - mean(EUCL))/sd(EUCL),  ylim.prim[1],
           1,  (ylim.sec[2] - mean(EUCL))/sd(EUCL), ylim.prim[2]))
a <- fit$coefficients['a']
s <- fit$coefficients['s']

sp.p <-ggplot(sp.df, aes(x=Group, y=Beta_mean, group=Taxa)) +
  geom_line(aes(color="#a4c64d"), color="#a4c64d")+
  geom_point(aes(color="#a4c64d", shape=Taxa, size=4), color="#a4c64d")+
  geom_point(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s), size=4), color="#000000")+
  geom_line(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s)), color="#000000")+
  scale_shape_manual(values=c(10,9,8))+
  scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~ (. - a) / s * sd(EUCL) + mean(EUCL), name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  # geom_errorbar(aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Taxa), width = 0.2) +
  # geom_errorbar(data = chem_data, aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000")
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(face="bold"))+
  ylab(" ")+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Spring Ice")

sp.p


#Summer Ice
sm.df = data2plot[data2plot$Habitat == "Summer Ice",]
chem_data = sm.df[sm.df$Taxa == "Prokaryote",]

#scale the data
#Taken from https://stackoverflow.com/questions/58774705/y-limits-for-ggplot-with-sec-axis
ylim.prim <- c(0, 1)
ylim.sec <- c(6, 12)
EUCL <- sm.df$Chem_mean #needed for coherent normalisation
# This is quite hacky, but it works if you want to set a boundary for the secondary y-axis
fit = lm(b ~ . + 0, 
         tibble::tribble(
           ~a, ~s,  ~b,
           1,  (ylim.sec[1] - mean(EUCL))/sd(EUCL),  ylim.prim[1],
           1,  (ylim.sec[2] - mean(EUCL))/sd(EUCL), ylim.prim[2]))
a <- fit$coefficients['a']
s <- fit$coefficients['s']

sm.p <-ggplot(sm.df, aes(x=Group, y=Beta_mean, group=Taxa)) +
  geom_line(aes(color="#4dd2d6"), color="#4dd2d6")+
  geom_point(aes(color="#4dd2d6", shape=Taxa, size=4), color="#4dd2d6")+
  geom_point(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s), size=4), color="#000000")+
  geom_line(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s)), color="#000000")+
  scale_shape_manual(values=c(10,9,8))+
  scale_y_continuous(limits=c(0,1)) +
  expand_limits(y = 0)+
  theme_bw()+
  # geom_errorbar(aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Taxa), width = 0.2) +
  # geom_errorbar(data = chem_data, aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000")
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, size=10),
        plot.title = element_text(face="bold"))+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Summer Ice")

sm.p

#Cryoconite
cr.df = data2plot[data2plot$Habitat == "Cryoconite",]
chem_data = cr.df[cr.df$Taxa == "Prokaryote",]

#scale the data
#Taken from https://stackoverflow.com/questions/58774705/y-limits-for-ggplot-with-sec-axis
ylim.prim <- c(0, 1)
ylim.sec <- c(6, 12)
EUCL <- sm.df$Chem_mean #needed for coherent normalisation
# This is quite hacky, but it works if you want to set a boundary for the secondary y-axis
fit = lm(b ~ . + 0, 
         tibble::tribble(
           ~a, ~s,  ~b,
           1,  (ylim.sec[1] - mean(EUCL))/sd(EUCL),  ylim.prim[1],
           1,  (ylim.sec[2] - mean(EUCL))/sd(EUCL), ylim.prim[2]))
a <- fit$coefficients['a']
s <- fit$coefficients['s']

cr.p <-ggplot(cr.df, aes(x=Group, y=Beta_mean, group=Taxa)) +
  geom_line(aes(color="#d8a4ff"), color="#d8a4ff")+
  geom_point(aes(color="#d8a4ff", shape=Taxa, size=4), color="#d8a4ff")+
  geom_point(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s), size=4), color="#000000")+
  geom_line(data = chem_data, aes(x = Group, y = (a + ((Chem_mean - mean(EUCL))/sd(EUCL)) * s)), color="#000000")+
  scale_shape_manual(values=c(10,9,8))+
  scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~ (. - a) / s * sd(EUCL) + mean(EUCL), name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  # geom_errorbar(aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Taxa), width = 0.2) +
  # geom_errorbar(data = chem_data, aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000")
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1, size=10),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.title.y.left = element_blank(),
        plot.title = element_text(face="bold"))+
  ggtitle("Cryoconite")

cr.p


#put the plots together

#get legend

leg1 = get_legend(snow.p + theme(legend.position = "bottom", 
                                legend.title=element_blank(),
                                legend.text = element_text(size=15)) + 
                   guides(color = "none", size="none", 
                          shape=guide_legend(override.aes = list(size=6))))

#get a second legend by making a plot that will plot the euclidean distances on the primary axis
plot.df = cr.df 
plot.df$ChemShape = "Environmental Euclidean Distance"
plot4leg <-ggplot(plot.df, aes(x=Group, y=Chem_mean)) +
  geom_line(aes(color="#000000"), color="#000000")+
  geom_point(aes(color="#000000", shape = ChemShape), color="#000000", size=4)+
  theme_bw()+
  theme(legend.title = element_blank(), legend.text = element_text(size=15))+
  guides(shape=guide_legend(override.aes = list(size=6)))
plot4leg

leg2 = get_legend(plot4leg)

full.p = (snow.p + sp.p + sm.p + cr.p) / leg1 / leg2 + plot_layout(heights = c(1, 0.1, 0.1))
full.p

pdf("../results/coms-and-euclidean.pdf", width=10, height=8)
print(full.p)
dev.off()
