
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
pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds") 


###################################################################################################
###################################################################################################

#PCOA USING UNWEIGHTED UNIFRAC
#PROKARYOTES - ARCTIC 

ps = subset_samples(pro, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23, 24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-arctic-unweighted-unifrac-pcoa.pdf", width=8, height=6)
print(p1)
dev.off()

#PROKARYOTES - ANTARCTIC 

ps = subset_samples(pro, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")


pdf("../results/prokaryote-antarctic-unweighted-unifrac-pcoa.pdf", width=8, height=6)
print(p2)
dev.off()


#EUKARYOTES - ARCTIC 

ps = subset_samples(euk, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-arctic-unweighted-unifrac-pcoa.pdf", width=8, height=6)
print(p3)
dev.off()


#EuKARYOTES - ANTARCTIC 

ps = subset_samples(euk, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-antarctic-unweighted-unifrac-pcoa.pdf", width=8, height=6)
print(p4)
dev.off()

#put all the plots together

p1a = p1 + theme(legend.position = "none", legend.title = element_blank())
p2a = p2 +theme(legend.position = "none", legend.title = element_blank())
p3a = p3 + theme(legend.position = "bottom")
p4a = p4 + theme(legend.position = "bottom")

###

p13 = ggarrange(p1a, p3a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
p24 = ggarrange(p2a, p4a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

final.p = p13 / p24 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
print(final.p)

pdf("../results/separate-poles-uunifrac-pcoa.pdf", height=8, width=8)
print(final.p)
dev.off()

###################################################################################################
###################################################################################################

#PCOA USING WEIGHTED UNIFRAC

#PROKARYOTES - ARCTIC 

ps = subset_samples(pro, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23, 24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-arctic-weighted-unifrac-pcoa.pdf", width=8, height=6)
print(p1)
dev.off()

#PROKARYOTES - ANTARCTIC 

ps = subset_samples(pro, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-antarctic-weighted-unifrac-pcoa.pdf", width=8, height=6)
print(p2)
dev.off()


#EUKARYOTES - ARCTIC 

ps = subset_samples(euk, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-arctic-weighted-unifrac-pcoa.pdf", width=8, height=6)
print(p3)
dev.off()


#EuKARYOTES - ANTARCTIC 

ps = subset_samples(euk, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-antarctic-weighted-unifrac-pcoa.pdf", width=8, height=6)
print(p4)
dev.off()

#put all the plots together

p1a = p1 + theme(legend.position = "none", legend.title = element_blank())
p2a = p2 +theme(legend.position = "none", legend.title = element_blank())
p3a = p3 + theme(legend.position = "bottom")
p4a = p4 + theme(legend.position = "bottom")

###

p13 = ggarrange(p1a, p3a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
p24 = ggarrange(p2a, p4a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

final.p = p13 / p24 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
print(final.p)

pdf("../results/separate-poles-wunifrac-pcoa.pdf", height=8, width=8)
print(final.p)
dev.off()

###################################################################################################
###################################################################################################

#PCOA USING BRAY CURTIS

#PROKARYOTES - ARCTIC 

ps = subset_samples(pro, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23, 24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-arctic-bray-pcoa.pdf", width=8, height=6)
print(p1)
dev.off()

#PROKARYOTES - ANTARCTIC 

ps = subset_samples(pro, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-antarctic-bray-pcoa.pdf", width=8, height=6)
print(p2)
dev.off()


#EUKARYOTES - ARCTIC 

ps = subset_samples(euk, Pole == "Arctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,23,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-arctic-bray-pcoa.pdf", width=8, height=6)
print(p3)
dev.off()


#EuKARYOTES - ANTARCTIC 

ps = subset_samples(euk, Pole == "Antarctic")

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Southern Victoria Land", "Queen Maud Land"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Region),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-antarctic-bray-pcoa.pdf", width=8, height=6)
print(p4)
dev.off()

#put all the plots together

p1a = p1 + theme(legend.position = "none", legend.title = element_blank())
p2a = p2 +theme(legend.position = "none", legend.title = element_blank())
p3a = p3 + theme(legend.position = "bottom")
p4a = p4 + theme(legend.position = "bottom")

###

p13 = ggarrange(p1a, p3a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
p24 = ggarrange(p2a, p4a, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

final.p = p13 / p24 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
print(final.p)

pdf("../results/separate-poles-bray-pcoa.pdf", height=8, width=8)
print(final.p)
dev.off()

###################################################################################################
###################################################################################################

#Ordination by group (Prokaryotes/Eukaryotes) with PCoA and unweighted UniFrac distances

#PROKARYOTES 

ps = pro

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p5 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-unweighted-unifrac-pcoa-both-poles.pdf", width=8, height=6)
print(p5)
dev.off()


#EUKARYOTES 

ps = euk

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p6 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-unweighted-unifrac-pcoa-both-poles.pdf", width=8, height=6)
print(p6)
dev.off()


#put the plots together

final.p2 = ggarrange(p5, p6, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = "AUTO",
                     font.label=list(size=20))
print(final.p2)

pdf("../results/polar-uunifrac-pcoa.pdf", height=5, width=15)
print(final.p2)
dev.off()

###################################################################################################
###################################################################################################

#Ordination by group (Prokaryotes/Eukaryotes) with PCoA and weighted UniFrac distances

#PROKARYOTES 

ps = pro

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p7 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-weighted-unifrac-pcoa-both-poles.pdf", width=8, height=6)
print(p7)
dev.off()


#EUKARYOTES 

ps = euk

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordUU <- ordinate(ps.RA, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ordUU$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordUU$values$Relative_eig[1]), 100*sum(ordUU$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p8 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-weighted-unifrac-pcoa-both-poles.pdf", width=8, height=6)
print(p8)
dev.off()

#put the plots together

final.p3 = ggarrange(p7, p8, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = "AUTO", font.label=list(size=20))
print(final.p3)

pdf("../results/polar-wunifrac-pcoa.pdf", height=5, width=15)
print(final.p3)
dev.off()

###################################################################################################
###################################################################################################

#Ordination by group (Prokaryotes/Eukaryotes) with PCoA and Bray-Curtis dissimilarities

#PROKARYOTES 

ps = pro

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p7 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Prokaryote")

pdf("../results/prokaryote-bray-pcoa-both-poles.pdf", width=8, height=6)
print(p7)
dev.off()


#EUKARYOTES 

ps = euk

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )


#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="PCoA", distance="bray")

#extract data for plotting
positions <- ordBR$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ordBR$values$Relative_eig[1]), 100*sum(ordBR$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.RA))$Subcommunity, data.frame(sample_data(ps.RA))$Region, data.frame(sample_data(ps.RA))$Pole)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Subcommunity", "Region", "Pole")

#default legend order
data2plot$Region <- factor(data2plot$Region, levels=c("Greenland", "Svalbard", "Sweden", "Southern Victoria Land", "Queen Maud Land"))
data2plot$Pole <- factor(data2plot$Pole, levels=c("Arctic", "Antarctic"))
data2plot$Subcommunity <- factor(data2plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#plot the data
p8 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=4)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(22,24))+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  stat_ellipse(aes(x=pcoa1, y=pcoa2, color=Subcommunity),type = "norm")+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20),
        plot.title = element_text(size=20))+
  ggtitle("Eukaryote")

pdf("../results/eukaryote-bray-pcoa-both-poles.pdf", width=8, height=6)
print(p8)
dev.off()

#put the plots together

final.p3 = ggarrange(p7, p8, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = "AUTO", font.label=list(size=20))
print(final.p3)

pdf("../results/polar-bray-pcoa.pdf", height=5, width=15)
print(final.p3)
dev.off()

###################################################################################################
###################################################################################################