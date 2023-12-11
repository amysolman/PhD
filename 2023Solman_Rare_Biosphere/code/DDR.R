#Methods

#The rate of distance-decay was calculated as the slope of the least squares linear regression on log-transformed geographic distance + 0.001 and log-transformed Bray-Curtis dissimilarity. Geographic distance in kilometers was calculated from latitude and longitude coordinates using the earth.dist function of the fossil package. 


# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Calculate geographic distance between sites
# 4. Calculate Bray-Curtis dissimilarities between sites
# 5. Carry out least-squares linear regression on log-transformed geographic distance and log-transformed community dissimilarity
# 6. Carry out a metric permutation test to test statistical significance of distance-decay slope.

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


# phylo1 = pro.ant
# phylo2 = pro.arc
# phylo = pro
# subcom = "Total"
# group = "Prokaryote"

# phylo = euk
# phylo1 = euk.ant.int
# phylo2 = euk.arc.int
# subcom = "Intermediate"
# group = "Eukaryote"

lm_func_DDR <- function(phylo, phylo1, phylo2, subcom, group){
  
  # 3. Calculate Bray-Curtis dissimilarities between sites
  
  #get out count tables
  x1 = data.frame(otu_table(phylo1), check.names=FALSE)
  x2 = data.frame(otu_table(phylo2), check.names=FALSE)
  
  #create empty dataframe with ASVs as rows and sample names as columns
  nams = c(colnames(x1), colnames(x2)) #sample names
  row.nams = taxa_names(phylo) #ASV IDs
  final.df <- data.frame(matrix(ncol = length(nams), nrow = length(row.nams)))
  colnames(final.df) = nams
  row.names(final.df) = row.nams
  final.df[is.na(final.df)] <- 0
  final.df <- tibble::rownames_to_column(final.df, "ASV") #make rownames a column
  x1 <- tibble::rownames_to_column(x1, "ASV")
  x2 <- tibble::rownames_to_column(x2, "ASV")
  
  #add in the count data
  final.df.plus = final.df %>% 
    rows_update(x1, by = "ASV")
  final.df.plus = final.df.plus %>% 
    rows_update(x2, by = "ASV")
  
  final.df.plus2 <- final.df.plus[,-1]  #remove ASV ID column
  rownames(final.df.plus2) <- final.df.plus[,1] #make IDs rownames
  
  sum(colSums(final.df.plus2)) == sum(colSums(x1[,-1])) + sum(colSums(x2[,-1])) #this should be true
  
  #remove rows of columns with zero counts
  comm.out = final.df.plus2[rowSums(final.df.plus2[])>0,]
  comm.out = comm.out[,colSums(comm.out)>0]
  min(rowSums(comm.out)) #should be at least 1
  min(colSums(comm.out)) #should be at least 1
  #comm.out = t(comm.out) #make samples rows and ASVs as columns
  
  #make into relative abundance before calculating bray-curtis as this assumes the same sample size
  comm.rel <- make_relative(t(comm.out))
  rowSums(comm.rel) #should all equal 1
  
  comm.dist<-vegdist(comm.rel, method="bray") #compute dissimilarities
  
  # 4. Calculate geographic distance between sites
  #data = data.frame(sample_data(phylo))
  geo1 = get_my_long_lat(phylo1)
  geo2 = get_my_long_lat(phylo2)
  geo.comb = rbind(geo1, geo2)
  
  #make sure we've got the same samples in our geo data and com dissimilarity data
  geo.comb = geo.comb[rownames(geo.comb) %in% names(comm.out),]
  
  nrow(geo.comb) == nrow(t(comm.out)) #should be true
  rownames(geo.comb) == rownames(t(comm.out)) #should all be true
  
  geo.dist <- earth.dist(geo.comb, dist=TRUE) #distance in kilometers
  
  # check.geo = matrix(geo.dist)
  
  # 5. Carry out least-squares linear regression on log-transformed geographic distance and log-transformed community dissimilarity
  
  # comm.dist = comm.dist
  # loggeo.dist = geo.dist
  
  #transform these data into columns
  comm.dist.df = melt(as.matrix(comm.dist))
  names(comm.dist.df) = c("Sample1", "Sample2", "Beta")
  comm.dist.df.trim = comm.dist.df[comm.dist.df['Sample1'] != comm.dist.df['Sample2'],]
  
  geo.dist.df = melt(as.matrix(geo.dist))
  names(geo.dist.df) = c("Sample1", "Sample2", "km")
  geo.dist.df.trim = geo.dist.df[geo.dist.df['Sample1'] != geo.dist.df['Sample2'],]
  
  #make into dataframe
  data = geo.dist.df.trim
  data$Beta = comm.dist.df.trim$Beta
  data$Subcommunity = subcom
  data$Group = group
  
  plot(data$km, data$Beta) #we can see a somehwat linear relationship
  cor(data$km, data$Beta, method="spearman") #spearman's correlation
  
  #check for outliers and remove
  #boxplot(data$Beta) #visualise with boxplot
  
  #statistical methods
  #using the quantile() function to find the 25th and the 75th percentile of the dataset
  #use the IQR() function to get the difference of the 75th and 25th percentiles. 
  # Q <- quantile(data$Beta, probs=c(.25, .75), na.rm = FALSE)
  # iqr <- IQR(data$Beta)
  # #find the cutoff ranges beyond which data points are outliers
  # up <-  Q[2]+1.5*iqr # Upper Range  
  # low<- Q[1]-1.5*iqr # Lower Range
  # #remove the outliers
  # eliminated<- subset(data, data$Beta > (Q[1] - 1.5*iqr) & data$Beta < (Q[2]+1.5*iqr))
  # 
  # #visualise the data without outliers
  # #ggbetweenstats(eliminated, Beta, km, outlier.tagging = TRUE) 
  # 
  # #alternatively we can identify the outliers as such...
  # outliers <- boxplot(data$Beta, plot=FALSE)$out
  # #and then remove them
  # x<-data
  # x<- x[-which(data$Beta %in% outliers),]
  # 
  # #check plot and cor again
  # plot(x$km, x$Beta) #we can see a somewhat linear relationship
  # cor(x$km, x$Beta, method="spearman") #spearman's correlation
  
  # plot(eliminated$km, eliminated$Beta) #we can see a somewhat linear relationship
  # cor(eliminated$km, eliminated$Beta, method="spearman") #spearman's correlation
  
  #log transform our data before fitting the linear model
  data$km <- log10(data$km+0.001)
  data$Beta <- log10(data$Beta)
  
  
  #fit the model  
  fit <- lm(Beta ~ km, data = data)
  # #model results
  #summary(fit)
  
  res.list = list(data, fit)
  
  return(res.list)
  
}


pro.abun.res = lm_func_DDR(pro, pro.ant.abun, pro.arc.abun, "Abundant", "Prokaryote")
pro.int.res = lm_func_DDR(pro, pro.ant.int, pro.arc.int, "Intermediate", "Prokaryote")
pro.rare.res = lm_func_DDR(pro, pro.ant.rare, pro.arc.rare, "Rare", "Prokaryote")
euk.abun.res = lm_func_DDR(euk, euk.ant.abun, euk.arc.abun, "Abundant", "Eukaryote")
euk.int.res = lm_func_DDR(euk, euk.ant.int, euk.arc.int, "Intermediate", "Eukaryote")
euk.rare.res = lm_func_DDR(euk, euk.ant.rare, euk.arc.rare, "Rare", "Eukaryote")

plot.df = rbind(pro.abun.res[[1]], pro.int.res[[1]], pro.rare.res[[1]], 
                euk.abun.res[[1]], euk.int.res[[1]], euk.rare.res[[1]])

summary(pro.abun.res[[2]])
summary(pro.int.res[[2]])
summary(pro.rare.res[[2]])
summary(euk.abun.res[[2]])
summary(euk.int.res[[2]])
summary(euk.rare.res[[2]])

abun.pro.tab <- data.frame(Subcommunity = "Abundant", Slope = round(pro.abun.res[[2]]$coefficients[[2]], 3), R2 = round(summary(pro.abun.res[[2]])[[8]], 3), P=anova(pro.abun.res[[2]])$'Pr(>F)'[1])
int.pro.tab <- data.frame(Subcommunity = "Intermediate", Slope = round(pro.int.res[[2]]$coefficients[[2]], 3), R2 = round(summary(pro.int.res[[2]])[[8]], 3), P=anova(pro.int.res[[2]])$'Pr(>F)'[1])
rare.pro.tab <- data.frame(Subcommunity = "Rare", Slope = round(pro.rare.res[[2]]$coefficients[[2]], 3), R2 = round(summary(pro.rare.res[[2]])[[8]], 3), P=anova(pro.rare.res[[2]])$'Pr(>F)'[1])
abun.euk.tab <- data.frame(Subcommunity = "Abundant", Slope = round(euk.abun.res[[2]]$coefficients[[2]], 3), R2 = round(summary(euk.abun.res[[2]])[[8]], 3), P=anova(euk.abun.res[[2]])$'Pr(>F)'[1])
int.euk.tab <- data.frame(Subcommunity = "Intermediate", Slope = round(euk.int.res[[2]]$coefficients[[2]], 3), R2 = round(summary(euk.int.res[[2]])[[8]], 3), P=anova(euk.int.res[[2]])$'Pr(>F)'[1])
rare.euk.tab <- data.frame(Subcommunity = "Rare", Slope = round(euk.rare.res[[2]]$coefficients[[2]], 3), R2 = round(summary(euk.rare.res[[2]])[[8]], 3), P=anova(euk.rare.res[[2]])$'Pr(>F)'[1])

pro.res.tab = rbind(abun.pro.tab, int.pro.tab, rare.pro.tab)
euk.res.tab = rbind(abun.euk.tab, int.euk.tab, rare.euk.tab)

all_res = rbind(pro.res.tab, euk.res.tab)
all_res = cbind(Group=c(rep("Prokaryote", 3), rep("Eukaryote", 3)), all_res)
write.csv(all_res, "../results/DDR-results.csv")

plot.df = plot.df %>%
  mutate(across(Group, factor, levels=c("Prokaryote","Eukaryote"))) %>%
  mutate(across(Subcommunity, factor, levels=c("Rare","Intermediate", "Abundant")))

p = ggplot(plot.df, aes(x=km, y=Beta, fill=Subcommunity, shape=Subcommunity)) +
  geom_point(data=plot.df, aes(x=km, y=Beta, fill=Subcommunity), alpha=.5, size=3)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, show.legend = FALSE, aes(color=Subcommunity))+
  facet_wrap(~Group)+
  theme_bw()+
  xlab("log10(km+0.001)")+
  ylab("log10(community dissimilarity)")+
  scale_shape_manual(values = c(21, 21, 21))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.text = element_text(size=15), 
        legend.title = element_blank())+
  guides(fill = guide_legend(override.aes = list(size = 7, alpha=.8) ) )

print(p)

pdf("../results/DDR.pdf")
print(p)
dev.off()
