
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

# 2. Import data
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 


#PROKARYOTE PCOA USING UNWEIGHTED UNIFRAC

#SNOW

#subset samples
sub = subset_samples(pro, Habitat == "Snow")
#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())
  #+ggtitle("Snow")

print(p1)


#SPRING ICE

#subset samples
sub = subset_samples(pro, Habitat == "Spring Ice")

#x = data.frame(sample_data(sub))

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())
  #+ggtitle("Spring Ice")

print(p2)



#SUMMER ICE

#subset samples
sub = subset_samples(pro, Habitat == "Summer Ice")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())
  #+ggtitle("Summer Ice")

print(p3)



#CRYOCONITE

#subset samples
sub = subset_samples(pro, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())
  #+ggtitle("Cryoconite")

print(p4)

#put all the plots together
pro.final.p = (p1 + p2) / (p3 + p4)

pdf("../results/pma-prokaryote-uunifrac-pcoa.pdf", height=8, width=10)
print(pro.final.p)
dev.off()


#EUKARYOTE PCOA USING UNWEIGHTED UNIFRAC

#EUKARYOTES

#SNOW

#subset samples
sub = subset_samples(euk, Habitat == "Snow")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p1)


#SPRING ICE

#subset samples
sub = subset_samples(euk, Habitat == "Spring Ice")

#x = data.frame(sample_data(sub))

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p2)



#SUMMER ICE

#subset samples
sub = subset_samples(euk, Habitat == "Summer Ice")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p3)



#CRYOCONITE

#subset samples
sub = subset_samples(euk, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="unifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p4)

#put all the plots together
euk.final.p = (p1 + p2) / (p3 + p4)

pdf("../results/pma-eukaryote-uunifrac-pcoa.pdf", height=8, width=10)
print(euk.final.p)
dev.off()

#PROKARYOTE PCOA USING WEIGHTED UNIFRAC

#SNOW

#subset samples
sub = subset_samples(pro, Habitat == "Snow")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p1)


#SPRING ICE

#subset samples
sub = subset_samples(pro, Habitat == "Spring Ice")

#x = data.frame(sample_data(sub))

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p2)

#SUMMER ICE

#subset samples
sub = subset_samples(pro, Habitat == "Summer Ice")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p3)

#CRYOCONITE

#subset samples
sub = subset_samples(pro, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p4)

#put all the plots together
pro.final.p = (p1 + p2) / (p3 + p4)

pdf("../results/pma-prokaryote-wunifrac-pcoa.pdf", height=8, width=10)
print(pro.final.p)
dev.off()

#EUKARYOTE PCOA USING WEIGHTED UNIFRAC

#SNOW

#subset samples
sub = subset_samples(euk, Habitat == "Snow")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p1)


#SPRING ICE

#subset samples
sub = subset_samples(euk, Habitat == "Spring Ice")

#x = data.frame(sample_data(sub))

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p2)



#SUMMER ICE

#subset samples
sub = subset_samples(euk, Habitat == "Summer Ice")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p3)



#CRYOCONITE

#subset samples
sub = subset_samples(euk, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="wunifrac")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p4)

#put all the plots together
euk.final.p = (p1 + p2) / (p3 + p4)

pdf("../results/pma-eukaryote-wunifrac-pcoa.pdf", height=8, width=10)
print(euk.final.p)
dev.off()

#PROKARYOTE PCOA USING BRAY CURTIS

#SNOW

#subset samples
sub = subset_samples(pro, Habitat == "Snow")
#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p1)


#SPRING ICE

#subset samples
sub = subset_samples(pro, Habitat == "Spring Ice")

#x = data.frame(sample_data(sub))

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p2)



#SUMMER ICE

#subset samples
sub = subset_samples(pro, Habitat == "Summer Ice")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p3)



#CRYOCONITE

#subset samples
sub = subset_samples(pro, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p4)

#look at our plots
p1
p2
p3
p4

#put all the plots together
pro.final.p = (p1 + ggtitle("Prokaryote") + p2 ) / (p3 + p4)
pro.final.p

pdf("../results/pma-prokaryote-bray-pcoa.pdf", height=8, width=10)
print(pro.final.p)
dev.off()


#EUKARYOTE PCOA USING BRAY CURTIS

#SNOW

#subset samples
sub = subset_samples(euk, Habitat == "Snow")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#default legend order
#data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#plot the data
p1 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p1)


#SPRING ICE
#subset samples
sub = subset_samples(euk, Habitat == "Spring Ice")
#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p2 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p2)



#SUMMER ICE

#subset samples
sub = subset_samples(euk, Habitat == "Summer Ice")
#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p3 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p3)



#CRYOCONITE

#subset samples
sub = subset_samples(euk, Habitat == "Cryoconite")

#remove taxa with zero counts
sub = filter_taxa(sub, function(x) sum(x) > 0, TRUE)

#remove samples with zero counts
sub = prune_samples(sample_sums(sub)>=1, sub)

#perform ordination
ord <- ordinate(sub, method="PCoA", distance="bray")

#extract data for plotting
positions <- ord$vectors[,1:2]
colnames(positions) <- c("pcoa1", "pcoa2")

#get percentages explained by the first 2 axis
percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

#create axis labels
labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
            glue("PCoA Axis 2 ({pretty_pe[2]}%)"))

#give sample names as row names and add regions and poles
x = data.frame(sample_data(sub))
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(sub))$Habitat, data.frame(sample_data(sub))$Treatment, data.frame(sample_data(sub))$SampleID)
names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "Habitat", "Treatment", "SampleID")

#plot the data
p4 = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleID, shape=Treatment),
             alpha=.5, size=3)+
  facet_grid(cols = vars(Habitat))+
  labs(x=labels[1], y=labels[2])+
  theme_bw()+
  scale_shape_manual(values = c(21,23))+
  #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
         shape = guide_legend(order = 2))+
  theme(legend.position="bottom", legend.title = element_blank())

print(p4)

p1
p2
p3
p4

#put all the plots together
euk.final.p = (p1 + ggtitle("Eukaryote") + p2 ) / (p3 + p4)
euk.final.p

pdf("../results/pma-eukaryote-bray-pcoa.pdf", height=8, width=10)
print(euk.final.p)
dev.off()

#Combine Prokaryote and eukaryote bray-curtis PCoA plots
# pro.final.p
# euk.final.p
# final.p = pro.final.p / euk.final.p

final.p = plot_grid(pro.final.p,
                    euk.final.p,
                    nrow=2)
final.p

pdf("../results/pma-pro-and-euk-bray-pcoa.pdf", height=15, width=10)
print(final.p)
dev.off()