
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

#Ordination by group (Prokaryotes/Eukaryotes) with NMDS and bray curtis dissimilarities

set.seed(666)
#PROKARYOTES

ps = pro

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="NMDS", distance="bray")
ordBR #stress 0.09

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
             alpha=.5, size=3)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  scale_shape_manual(values = c(23,24))+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  guides(fill=guide_legend(override.aes=list(shape=21, size=5, alpha=0.8), order = 1),
         shape = guide_legend(override.aes=list(size=5), order = 2), 
         color = "none")+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=Subcommunity),type = "norm")+
  scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text = element_text(size=20))+
  ggtitle("Prokaryote")

print(p.nmds.1)


#EUKARYOTES

ps = euk

#remove taxa or samples with zero counts
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)

#normalise the data by transforming into relative abundance because our subcommunity sizes will be different
ps.RA = transform_sample_counts(ps, function(x) x / sum(x) )

#perform ordination on unweighted unifrac distances
ordBR <- ordinate(ps.RA, method="NMDS", distance="bray")

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
             alpha=.5, size=3)+
  annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5,label=paste0("Stress:", round(ordBR$stress, 2)))+
  theme_bw()+
  scale_shape_manual(values = c(23,24))+
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

print(p.nmds.2)

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
keep.df$pole = factor(keep.df$pole, levels = c("Arctic", "Antarctic"))


#Combine plots

plot.a = ggarrange(p.nmds.1 + annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5, size=100,label=paste0("Stress:", round(ordBR$stress, 2))),
                   p.nmds.2 + annotate("text",x=-Inf,y=-Inf,hjust=-7.5,vjust=-0.5, size=100,label=paste0("Stress:", round(ordBR$stress, 2))),
                   ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                   font.label=list(size=20))

print(plot.a)

p.nmds.1
p.nmds.2
keep.df$Subcommunity = factor(keep.df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))
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

#plot the data
p7 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=3)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(23,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Prokaryote")

print(p7)


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

#plot the data
p8 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=Subcommunity, shape=Pole),
             alpha=.5, size=3)+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(23,24))+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())+
  ggtitle("Eukaryote")

print(p8)

#BETA DIVERSITY BOXPLOTS
#for saving our dataframes
keep.df = data.frame()


#####PROKARYOTES ARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

#####PROKARYOTES ANTARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ANTARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Abundant" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Rare" & Pole == "Antarctic"))
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

keep.df$group = factor(keep.df$group, levels = c("Prokaryote", "Eukaryote"))
keep.df$pole = factor(keep.df$pole, levels = c("Arctic", "Antarctic"))

p<-ggplot(keep.df, aes(x=Subcommunity, y=Beta, fill = Subcommunity)) +
  geom_boxplot()+
  #facet_grid(pole ~ group, margins=TRUE)
  facet_wrap(~ group + pole, scales = "free_x")+
  ylab("Unweighted Unifrac Distances")+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 20))
p


pdf("../results/uweighted-beta-diversity-boxplot.pdf", height=8, width=10)
print(p)
dev.off()


#Combine plots

plot.a = ggarrange(p5, p6, ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                   font.label=list(size=20))
print(plot.a)

keep.df$Subcommunity = factor(keep.df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))
plot.b = ggplot(keep.df, aes(x=Subcommunity, y=Beta, fill = Subcommunity), alpha=0.8) +
  geom_boxplot()+
  #facet_grid(pole ~ group, margins=TRUE)
  facet_wrap(~ group + pole, nrow = 1, scales = "free_x")+
  ylab("Unweighted Unifrac Distances")+
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

pdf("../results/unweighted-PCoA-and-beta-diversity-boxplot.pdf", height=8, width=10)
print(plot.c)
dev.off()

#for saving our dataframes
keep.df = data.frame()


#####PROKARYOTES ARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

#####PROKARYOTES ANTARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Prokaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Arctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)


#####EUKARYOTES ANTARCTIC##########

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Abundant" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Abundant"

keep.df = rbind(keep.df, df)


#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Intermediate"

keep.df = rbind(keep.df, df)

#get unweighted unifrac distance matrix and melt
dis = UniFrac(subset_samples(euk, Subcommunity == "Rare" & Pole == "Antarctic"), weighted=TRUE)
dis.melt = melt(as.matrix(dis))
names(dis.melt) = c("Sample1", "Sample2", "Beta")
#remove samples being compared to themseves
df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
df$group = "Eukaryote"
df$pole = "Antarctic"
df$Subcommunity = "Rare"

keep.df = rbind(keep.df, df)

keep.df$group = factor(keep.df$group, levels = c("Prokaryote", "Eukaryote"))
keep.df$pole = factor(keep.df$pole, levels = c("Arctic", "Antarctic"))

p<-ggplot(keep.df, aes(x=Subcommunity, y=Beta, fill = Subcommunity)) +
  geom_boxplot()+
  #facet_grid(pole ~ group, margins=TRUE)
  facet_wrap(~ group + pole, scales = "free_x")+
  ylab("Weighted Unifrac Distances")+
  scale_fill_manual(values=c("#ef7a76","#74a9d8", "#82ca81"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 20))
p


pdf("../results/weighted-beta-diversity-boxplot.pdf", height=8, width=10)
print(p)
dev.off()

#Combine plots

plot.a = ggarrange(p7, p8, ncol=2, nrow=1, common.legend = TRUE, legend="bottom",
                   font.label=list(size=20))
print(plot.a)

keep.df$Subcommunity = factor(keep.df$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))
plot.b = ggplot(keep.df, aes(x=Subcommunity, y=Beta, fill = Subcommunity), alpha=0.8) +
  geom_boxplot()+
  #facet_grid(pole ~ group, margins=TRUE)
  facet_wrap(~ group + pole, nrow = 1, scales = "free_x")+
  ylab("Weighted Unifrac Distances")+
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

pdf("../results/weighted-PCoA-and-beta-diversity-boxplot.pdf", height=8, width=10)
print(plot.c)
dev.off()

#Beta diversity partitioning analysis 
#calculate total beta diversity and beta diversity due to richness and turnover and plot

# phylo = pro_res[[7]] #arctic abundant community
# subcom = "Rare"
# data = "Prokaryote"
# pole = "Arctic"

beta_part <- function(phylo, subcom, data, pole){
  
  set.seed(666)
  
  #normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)
  
  #for each subcommunity get the count table with rows as samples and species as columns
  counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)
  
  #max(rowSums(counts)) - min(rowSums(counts))
  
  #beta partition analysis
  out = bray.part(counts)
  
  #convert the matrices into a dataframe to plot
  out.turn = melt(as.matrix(out$bray.bal))
  names(out.turn) = c("Sample1", "Sample2", "Beta")
  out.turn$Part = "Turnover"
  out.rich = melt(as.matrix(out$bray.gra))
  names(out.rich) = c("Sample1", "Sample2", "Beta")
  out.rich$Part = "Richness"
  out.whole = melt(as.matrix(out$bray))
  names(out.whole) = c("Sample1", "Sample2", "Beta")
  out.whole$Part = "Whole"
  
  #bind data frame
  three.df = rbind(out.whole, out.turn, out.rich)
  
  #remove samples being compared to themseves
  three.df.trim = three.df[three.df['Sample1'] != three.df['Sample2'],]
  
  #add abundance class
  three.df.trim$Subcommunity = subcom
  
  #and dataset
  three.df.trim$Data = data
  three.df.trim$Pole = pole
  
  #get proportions
  two.df = rbind(out.turn, out.rich)
  #remove samples being compared to themseves
  two.df.trim = two.df[two.df['Sample1'] != two.df['Sample2'],]
  #get proportions
  two.df.trim$Proportion = two.df.trim$Beta/sum(two.df.trim$Beta)
  sum(two.df.trim$Proportion)
  two.df.trim$Subcommunity = subcom
  two.df.trim$Data = data
  two.df.trim$Pole = pole
  
  res.list = list(three.df.trim, two.df.trim)
  
  return(res.list)
  
}

ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.abun.beta.part.arc = beta_part(ps, "Abundant", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.int.beta.part.arc = beta_part(ps, "Intermediate", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.rare.beta.part.arc = beta_part(ps, "Rare", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.abun.beta.part.ant = beta_part(ps, "Abundant", "Prokaryote", "Antarctic")

ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.int.beta.part.ant = beta_part(ps, "Intermediate", "Prokaryote", "Antarctic")

ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.rare.beta.part.ant = beta_part(ps, "Rare", "Prokaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.abun.beta.part.arc = beta_part(ps, "Abundant", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.int.beta.part.arc = beta_part(ps, "Intermediate", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.rare.beta.part.arc = beta_part(ps, "Rare", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.abun.beta.part.ant = beta_part(ps, "Abundant", "Eukaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.int.beta.part.ant = beta_part(ps, "Intermediate", "Eukaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.rare.beta.part.ant = beta_part(ps, "Rare", "Eukaryote", "Antarctic")

#make into one dataframe
final.beta = rbind(pro.abun.beta.part.arc[[1]], pro.int.beta.part.arc[[1]], pro.rare.beta.part.arc[[1]],
                   pro.abun.beta.part.ant[[1]], pro.int.beta.part.ant[[1]], pro.rare.beta.part.ant[[1]],
                   euk.abun.beta.part.arc[[1]], euk.int.beta.part.arc[[1]], euk.rare.beta.part.arc[[1]],
                   euk.abun.beta.part.ant[[1]], euk.int.beta.part.ant[[1]], euk.rare.beta.part.ant[[1]])

final.beta.barplot = rbind(pro.abun.beta.part.arc[[2]], pro.int.beta.part.arc[[2]], pro.rare.beta.part.arc[[2]],
                           pro.abun.beta.part.ant[[2]], pro.int.beta.part.ant[[2]], pro.rare.beta.part.ant[[2]],
                           euk.abun.beta.part.arc[[2]], euk.int.beta.part.arc[[2]], euk.rare.beta.part.arc[[2]],
                           euk.abun.beta.part.ant[[2]], euk.int.beta.part.ant[[2]], euk.rare.beta.part.ant[[2]])


# ggplot() + 
#    geom_histogram(aes(Petal.Width))+ 
#    facet_grid(Species~.)

#using KRuskal-Wallist Test are there significant differences in the proportion of dissimilarity account for Turnover and Richness between the subcommunities?
#Test Turnover
turn.df.pro.arc = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(turn.df.pro.arc$Beta, turn.df.pro.arc$Subcommunity,
                     p.adjust.method = "BH")
turn.df.pro.ant = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(turn.df.pro.ant$Beta, turn.df.pro.ant$Subcommunity,
                     p.adjust.method = "BH")
turn.df.euk.arc = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(turn.df.euk.arc$Beta, turn.df.euk.arc$Subcommunity,
                     p.adjust.method = "BH")
turn.df.euk.ant = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(turn.df.euk.ant$Beta, turn.df.euk.ant$Subcommunity,
                     p.adjust.method = "BH")

rich.df.pro.arc = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(rich.df.pro.arc$Beta, rich.df.pro.arc$Subcommunity,
                     p.adjust.method = "BH")
rich.df.pro.ant = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(rich.df.pro.ant$Beta, rich.df.pro.ant$Subcommunity,
                     p.adjust.method = "BH")
rich.df.euk.arc = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(rich.df.euk.arc$Beta, rich.df.euk.arc$Subcommunity,
                     p.adjust.method = "BH")
rich.df.euk.ant = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(rich.df.euk.ant$Beta, rich.df.euk.ant$Subcommunity,
                     p.adjust.method = "BH")

final.beta %>%
  filter(Part == "Whole") %>%
  mutate(across(Data, factor, levels=c("Prokaryote","Eukaryote"))) %>%
  mutate(across(Pole, factor, levels=c("Arctic","Antarctic"))) %>%
  ggplot(aes(x=Subcommunity, y=Beta, fill=Part)) +
  geom_boxplot(position=position_dodge(1))+
  facet_wrap(~ Data + Pole, scales = "free_x")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Bray-Curtis Dissimilarity")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        legend.text = element_text(size=10), 
        legend.title = element_blank())

p = final.beta.barplot %>%
  mutate(across(Data, factor, levels=c("Prokaryote","Eukaryote"))) %>%
  mutate(across(Pole, factor, levels=c("Arctic","Antarctic"))) %>%
  mutate(across(Subcommunity, factor, levels=c("Rare", "Intermediate", "Abundant"))) %>%
  group_by(Part, Subcommunity, Pole, Data) %>%
  dplyr::summarise(across(c(Proportion), sum)) %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.y=element_text(size=10,face="bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())

pdf("../results/beta-diversity-partitioning.pdf")
print(p)
dev.off()


#get percentages from final.beta.barplot dataframe

#merge proportion by part, subcommunity, data and pole

df = final.beta.barplot %>% group_by(Part, Subcommunity, Data, Pole) %>% summarise_each(funs(sum))

#round(df[df$Part == "Richness" & df$Pole == "Arctic" & df$Subcommunity == "Abundant" & df$Data == "Prokaryote",]$Proportion*100, 2)

x = data.frame(cbind(df$Part, df$Subcommunity, df$Data, df$Pole, round(df$Proportion*100, 2)))
names(x) = c("Partition", "Subcommunity", "Dataset", "Pole", "Proportion")

rownames(x) <- NULL
knitr::kable(x, caption = 'Table 7. Beta-diversity partitioning analysis results.')

###Beta diversity partitioning analysis

df.keep = data.frame()

#Abundant Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Arctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())






#Abundant Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Antarctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.x=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())






#Abundant Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Eukaryote"
pole = "Arctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Eukaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Eukaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)

x = df.keep[df.keep$Data == "Eukaryote" & df.keep$Part == "Richness",]

#barchart

p = df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())
pdf("../results/test-barplot.pdf", height = 30)
print(p)
dev.off()


#Abundant Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Antarctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.x=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())