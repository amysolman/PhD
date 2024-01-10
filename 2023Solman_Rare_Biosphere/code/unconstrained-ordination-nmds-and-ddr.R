
rm(list=ls())
source("00-solman-functions.R")
library(ggplot2)
library(phyloseq)
library(patchwork)
library(ggpubr)
library(glue)
library(tibble) # for add_column
library(reshape2)
library(cowplot)
library(betapart)
library(dplyr)
library(stringr) #for subsetting strings in lat long function
library(fossil) #for earth.dist function
library(vegan) #for community dissimilarity matrix
library(funrar)
library(ggpmisc)
library(picante) #for faith's pd

# 3. Import data and separate into Arctic and Antarctic samples
# 2. Import data
pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds") 

###################################################################################################
###################################################################################################

#Ordination by group (Prokaryotes/Eukaryotes) with NMDS and bray curtis dissimilarities

set.seed(666)

#PROKARYOTES

ps = pro

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on bray curtis distances
ordBR <- ordinate(ps.RA, method="NMDS", distance="bray")
pro_stress = round(ordBR$stress,3) #stress 0.09

#extract data for plotting
positions <- ordBR$points[,1:2]
colnames(positions) <- c("NMDS1", "NMDS2")

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "NMDS1", "NMDS2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))

data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p.nmds.1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill= "none",
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20))+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-bray-nmds-both-poles.pdf", width=8, height=6)
print(p.nmds.1)
dev.off()

#save our data for putting plots together
nmds1.data2plot = data2plot

set.seed(666)
#EUKARYOTES

ps = euk

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on bray curtis distances
ordBR <- ordinate(ps.RA, method="NMDS", distance="bray")
euk_stress = round(ordBR$stress,3) #stress 0.09

#extract data for plotting
positions <- ordBR$points[,1:2]
colnames(positions) <- c("NMDS1", "NMDS2")

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "NMDS1", "NMDS2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
#plot the data
p.nmds.2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill= "none",
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20))+
  ggtitle("Eukaryote")
#geom_text(x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))

pdf("../results/eukaryote-bray-nmds-both-poles.pdf", width=8, height=6)
print(p.nmds.2)
dev.off()

#save our data for putting plots together
nmds2.data2plot = data2plot

#put the plots together
nmds1.data2plot$Group = "Prokaryote"
nmds2.data2plot$Group = "Eukaryote"
nmds3.data2plot = as.data.frame(rbind(nmds1.data2plot, nmds2.data2plot))
nmds3.data2plot$Group = factor(nmds3.data2plot$Group, levels=c("Prokaryote", "Eukaryote"))

p.nmds.3 = ggplot(data=nmds3.data2plot) +
  geom_point(data=nmds3.data2plot, aes(x=NMDS1, y=NMDS2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  facet_wrap(~Group, scales="free")+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill= "none",
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=15))
#geom_text(x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))

pdf("../results/pro-euk-bray-nmds-both-poles.pdf", width=8, height=6)
print(p.nmds.3)
dev.off()

###################################################################################################
###################################################################################################

#Report NMDS Stress Results
sink("../results/nmds-stress.txt", type="output")
writeLines("===============================================================
NMDS STRESS RESULTS
===============================================================")
writeLines("NMDS stress values below 0.05 indicate a good fit, 0.05-0.1 are a fair fit, 0.1-0.2 is a poor fit, 0.2-0.3 is a failed fit.") 
writeLines("NMDS stress for prokaryote dataset with Bray-Curtis Dissimilarities:")
pro_stress
writeLines("NMDS stress for eukaryote dataset with Bray-Curtis Dissimilarities:")
euk_stress
sink()

###################################################################################################
###################################################################################################

#Combine with boxplot

#BETA DIVERSITY BOXPLOTS
#for saving our dataframes
keep.df = data.frame()


#####PROKARYOTES ARCTIC##########

#get bray-curtis distance matrix and melt
sub = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get bray dissimilarities matrix and melt
sub = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get bray dissimilarities matrix and melt
sub = subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

#####PROKARYOTES ANTARCTIC##########

#get bray dissimilarities matrix and melt
sub = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get bray dissimilarities matrix and melt
sub = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get bray dissimilarities matrix and melt
sub = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ARCTIC##########

#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ANTARCTIC##########

#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get bray dissimilarities matrix and melt
sub = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
rel = transform_sample_counts(sub, function(x) x/sum(x))
dis <- phyloseq::distance(sub, method = "bray") 
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

keep.df$group = factor(keep.df$group, levels = c("Prokaryote", "Eukaryote"))
keep.df$pole = factor(keep.df$pole, levels = c("Antarctic", "Arctic"))


#Combine plots

plot.a = ggarrange(p.nmds.1 + annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5, size=100,label=paste0("Stress:", round(ordBR$stress, 2))),
                   p.nmds.2 + annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5, size=100,label=paste0("Stress:", round(ordBR$stress, 2))),
                   ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                   font.label=list(size=20))

print(plot.a)


keep.df$Subcommunity = factor(keep.df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

#beta diversity boxplot
plot.b = ggplot(keep.df, aes(x=Subcommunity, y=Beta, fill = Subcommunity), alpha=0.8) +
  geom_boxplot()+
  #facet_grid(pole ~ group, margins=TRUE)
  facet_wrap(~ group + pole, nrow = 1, scales = "free_x")+
  ylab("Bray-Curtis Dissimilarities")+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 18))
print(plot.b)

# plot.c = plot.a / plot.b + plot_annotation(tag_levels = "A") & 
# theme(plot.tag = element_text(size = 25))
# print(plot.c)

plot.c = plot_grid(plot.a, plot.b, nrow=2, labels="AUTO", label_size = 20)
print(plot.c)

pdf("../results/bray-NMDS-and-beta-diversity-boxplot.pdf", height=8, width=10)
print(plot.c)
dev.off()

###################################################################################################
###################################################################################################

#DISTANCE DECAY RELATIONSHIP

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

plot.d = ggplot(plot.df, aes(x=km, y=Beta, fill=Subcommunity, shape=Subcommunity)) +
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

print(plot.d)

pdf("../results/DDR.pdf")
print(plot.d)
dev.off()


#put all the plots together
pa = p.nmds.3 + theme(strip.text = element_text(size=10))
pb = plot.b + theme(strip.text.x = element_text(size=10))
pc = plot.d + theme(strip.text = element_text(size=10)) + guides(fill = guide_legend(override.aes = list(size = 7, alpha=.8, shape=22) ) )

#make titles
title1 <- ggdraw() + draw_label("NMDS Ordinations", size=14, hjust=2.69, fontface = "bold")
title2 <- ggdraw() + draw_label("Bray-Curtis Dissimilarities", size=14, hjust=1.85,fontface = "bold")
title3 <- ggdraw() + draw_label("Distance-Decay Relationships", size=14, hjust=1.62,fontface = "bold")

final.plot = title1 / pa / plot_spacer() /title2 / pb / plot_spacer() /title3 / pc + plot_layout(heights = c(0.1,1,0.02,0.1,1.2,0.03,0.1,1.2))

pdf("../results/beta-diversity-plots.pdf", height=14, width=10)
final.plot
dev.off()