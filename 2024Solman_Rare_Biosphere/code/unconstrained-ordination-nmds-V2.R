
rm(list=ls())
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

# 3. Import data and separate into Arctic and Antarctic samples
# 2. Import data
# pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
# euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds") 

pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds") 

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
  add_column(data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole, data.frame(sample_data(ps.RA))$Glacier)
names(data2plot) = c("Sample", "NMDS1", "NMDS2",  "Region", "Pole", "Glacier")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))

#data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
# p.nmds.1 = ggplot(data=data2plot) +
#   geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Subcommunity, shape=Pole),
#              alpha=.5, size=4)+
#   annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
#   theme_bw()+
#   scale_shape_manual(values = c(22,24))+
#   scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
#   guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
#          shape = guide_legend(override.aes=list(size=5), order = 2), 
#          color = "none")+
#   stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Subcommunity),type = "norm")+
#   scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
#   theme(legend.position="bottom", legend.title = element_blank(),
#         legend.text = element_text(size=20))+
#   ggtitle("Prokaryote")
# 
# pdf("../results/prokaryote-bray-nmds-both-poles.pdf", width=8, height=6)
# print(p.nmds.1)
# dev.off()

#plot the data
p.nmds.1a = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Glacier, shape=Pole),
             alpha=.5, size=4)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  #scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  #stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Glacier),type = "norm")+
  #scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20))+
  ggtitle("Prokaryote")

p.nmds.1a

pdf("../results/prokaryote-bray-nmds-regions.pdf", width=8, height=6)
print(p.nmds.1a)
dev.off()

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
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
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
        strip.text.x = element_text(size = 20))
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


