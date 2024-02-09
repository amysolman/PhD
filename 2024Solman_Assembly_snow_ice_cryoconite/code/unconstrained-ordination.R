# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(RColorBrewer)
# library(tidyr) #trans data between wide and long format
library(ggplot2)
library(vegan) #anosim function
# source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
# library(dplyr)
library(cowplot)
# library(reshape) 
# library(funrar) #for make relative
# library(stringr) #mutate function, to split columns
# library(gridExtra) #for exporting as pdf
# library(scales)
# library(microbiome) #for summarize_phyloseq
library(tibble) # for add_column
library(patchwork)
#library(glue)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#set seed for nmds
set.seed(666)

#NMDS Prokaryotes
#set habitat location variable
sample_data(ps.pro)$Habitat_location <- paste0(sample_data(ps.pro)$Habitat, " (", sample_data(ps.pro)$Location, " Foxfonna)")
sample_data(ps.pro)$Location <- paste0(sample_data(ps.pro)$Location, " Foxfonna")
#perform NMDS
pro.nmds.ord <- phyloseq::ordinate(ps.pro, "NMDS", distance="bray")
#get stress
pro.nmds.stress = round(pro.nmds.ord$stress, 3)
#extract data for plotting
positions <- pro.nmds.ord$points
#get sample data for plotting
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.pro))$Habitat_location, data.frame(sample_data(ps.pro))$Habitat, data.frame(sample_data(ps.pro))$Location, "Prokaryote")
names(data2plot) = c("Sample", "NMDS1", "NMDS2", "Habitat Location", "Habitat", "Location", "Group")
#get colours for plotting
brewer.pal(n = 4, name = "Set1")
#default legend order
data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
#plot
pro.nmds.p = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, fill=Habitat, shape=Location),
             alpha=.5, size=4)+
  scale_shape_manual(values = c(21,24))+
  xlim(-1.5, 4.5)+
  theme_bw()+
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle(expression(bold("Prokaryote")))

print(pro.nmds.p)


#NMDS Eukaryotes
#set habitat location variable
sample_data(ps.euk)$Habitat_location <- paste0(sample_data(ps.euk)$Habitat, " (", sample_data(ps.euk)$Location, " Foxfonna)")
sample_data(ps.euk)$Location <- paste0(sample_data(ps.euk)$Location, " Foxfonna")
#perform NMDS
euk.nmds.ord <- phyloseq::ordinate(ps.euk, "NMDS", distance="bray")
#get stress
euk.nmds.stress = round(euk.nmds.ord$stress, 3)
#extract data for plotting
positions <- euk.nmds.ord$points
#get sample data for plotting
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.euk))$Habitat_location, data.frame(sample_data(ps.euk))$Habitat, data.frame(sample_data(ps.euk))$Location, "Eukaryote")
names(data2plot) = c("Sample", "NMDS1", "NMDS2", "Habitat Location", "Habitat", "Location", "Group")
#get colours for plotting
brewer.pal(n = 4, name = "Set1")
#default legend order
data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
#plot
euk.nmds.p = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, shape=Location, fill= Habitat),
             alpha=.5, size=4)+
  scale_shape_manual(values = c(21,24))+
  xlim(-1.5, 4.5)+
  guides(fill=guide_legend(override.aes=list(shape=22)))+
  theme_bw()+
  theme(legend.position=c(.8,.55),
        legend.box="vertical", legend.margin=margin(),
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
ggtitle(expression(bold("Microbial Eukaryote")))

print(euk.nmds.p)


#NMDS Microfauna
#set habitat location variable
sample_data(ps.mm)$Habitat_location <- paste0(sample_data(ps.mm)$Habitat, " (", sample_data(ps.mm)$Location, " Foxfonna)")
sample_data(ps.mm)$Location <- paste0(sample_data(ps.mm)$Location, " Foxfonna")
#perform NMDS
mm.nmds.ord <- phyloseq::ordinate(ps.mm, "NMDS", distance="bray")
#get stress
mm.nmds.stress = round(mm.nmds.ord$stress, 3)
#extract data for plotting
positions <- mm.nmds.ord$points
#get sample data for plotting
data2plot = positions %>%
  as_tibble(rownames = "samples") %>%
  add_column(data.frame(sample_data(ps.mm))$Habitat_location, data.frame(sample_data(ps.mm))$Habitat, data.frame(sample_data(ps.mm))$Location, "Microfauna")
names(data2plot) = c("Sample", "NMDS1", "NMDS2", "Habitat Location", "Habitat", "Location", "Group")
#get colours for plotting
brewer.pal(n = 4, name = "Set1")
#default legend order
data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
#plot
mm.nmds.p = ggplot(data=data2plot) +
  geom_point(data=data2plot, aes(x=NMDS1, y=NMDS2, shape=Location, fill= Habitat),
             alpha=.5, size=4)+
  scale_shape_manual(values = c(21,24))+
  xlim(-1.5, 4.5)+
  theme_bw()+
  theme(legend.position="none")+
ggtitle(expression(bold("Microfauna")))


print(mm.nmds.p)

#plot with patchwork
nmds.all.p1 = pro.nmds.p / euk.nmds.p / mm.nmds.p
nmds.all.p1

pdf("../results/nmds-plot-all-groups.pdf", height=8, width=8)
print(nmds.all.p1)
dev.off()

# #PCoA Prokaryotes
# #set habitat location variable
# sample_data(ps.pro)$Habitat_location <- paste0(sample_data(ps.pro)$Habitat, " (", sample_data(ps.pro)$Location, " Foxfonna)")
# #perform PCoA
# pro.PCoA.ord <- phyloseq::ordinate(ps.pro, "PCoA", distance="bray")
# #get percentages explained by the axis
# percent_explained <- 100 * pro.PCoA.ord$values$Eigenvalues / sum(pro.PCoA.ord$values$Eigenvalues)
# pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
# #create axis labels
# labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
#             glue("PCoA Axis 2 ({pretty_pe[2]}%)"))
# #extract data for plotting
# positions <- pro.PCoA.ord$vectors[,1:2]
# #get sample data for plotting
# data2plot = positions %>%
#   as_tibble(rownames = "samples") %>%
#   add_column(data.frame(sample_data(ps.pro))$Habitat_location, 
#              data.frame(sample_data(ps.pro))$Habitat, 
#              data.frame(sample_data(ps.pro))$Location, "Prokaryote")
# names(data2plot) = c("Sample", "PCoA1", "PCoA2", "Habitat Location", "Habitat", "Location", "Group")
# #get colours for plotting
# brewer.pal(n = 4, name = "Set1")
# #default legend order
# data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
# #plot
# pro.PCoA.p = ggplot(data=data2plot) +
#   geom_point(data=data2plot, aes(x=PCoA1, y=PCoA2, color=Location, shape=Habitat, fill= Location),
#              alpha=.5, size=4)+
#   scale_colour_manual(name = "",
#                       labels = c("Upper", "Lower"),
#                       values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   scale_fill_manual(name = "",
#                     labels = c("Upper", "Lower"),
#                     values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   
#   scale_shape_manual(name = "",
#                      labels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"),
#                      values = c(22, 24, 25, 23))+
#   labs(x=labels[1], y=labels[2])+
#   theme_bw()+
#   theme(legend.position="none")+
#   ggtitle("Prokaryotes")
# 
# print(pro.PCoA.p)
# 
# 
# #PCoA Eukaryotes
# #set habitat location variable
# sample_data(ps.euk)$Habitat_location <- paste0(sample_data(ps.euk)$Habitat, " (", sample_data(ps.euk)$Location, " Foxfonna)")
# #perform PCoA
# euk.PCoA.ord <- phyloseq::ordinate(ps.euk, "PCoA", distance="bray")
# #get percentages explained by the axis
# percent_explained <- 100 * euk.PCoA.ord$values$Eigenvalues / sum(euk.PCoA.ord$values$Eigenvalues)
# pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
# #create axis labels
# labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
#             glue("PCoA Axis 2 ({pretty_pe[2]}%)"))
# #extract data for plotting
# positions <- euk.PCoA.ord$vectors[,1:2]
# #get sample data for plotting
# data2plot = positions %>%
#   as_tibble(rownames = "samples") %>%
#   add_column(data.frame(sample_data(ps.euk))$Habitat_location, data.frame(sample_data(ps.euk))$Habitat, data.frame(sample_data(ps.euk))$Location, "Eukaryote")
# names(data2plot) = c("Sample", "PCoA1", "PCoA2", "Habitat Location", "Habitat", "Location", "Group")
# #get colours for plotting
# brewer.pal(n = 4, name = "Set1")
# #default legend order
# data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
# #plot
# euk.PCoA.p = ggplot(data=data2plot) +
#   geom_point(data=data2plot, aes(x=PCoA1, y=PCoA2, color=Location, shape=Habitat, fill= Location),
#              alpha=.5, size=4)+
#   scale_colour_manual(name = "",
#                       labels = c("Upper", "Lower"),
#                       values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   scale_fill_manual(name = "",
#                     labels = c("Upper", "Lower"),
#                     values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   
#   scale_shape_manual(name = "",
#                      labels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"),
#                      values = c(22, 24, 25, 23))+
#   labs(x=labels[1], y=labels[2])+
#   theme_bw()+
#   theme(legend.position="none")+
#   ggtitle("Eukaryotes")
# 
# print(euk.PCoA.p)
# 
# 
# #PCoA Microfauna
# #set habitat location variable
# sample_data(ps.mm)$Habitat_location <- paste0(sample_data(ps.mm)$Habitat, " (", sample_data(ps.mm)$Location, " Foxfonna)")
# #perform PCoA
# mm.PCoA.ord <- phyloseq::ordinate(ps.mm, "PCoA", distance="bray")
# #get percentages explained by the axis
# percent_explained <- 100 * mm.PCoA.ord$values$Eigenvalues / sum(mm.PCoA.ord$values$Eigenvalues)
# pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
# #create axis labels
# labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
#             glue("PCoA Axis 2 ({pretty_pe[2]}%)"))
# #extract data for plotting
# positions <- mm.PCoA.ord$vectors[,1:2]
# #get sample data for plotting
# data2plot = positions %>%
#   as_tibble(rownames = "samples") %>%
#   add_column(data.frame(sample_data(ps.mm))$Habitat_location, data.frame(sample_data(ps.mm))$Habitat, data.frame(sample_data(ps.mm))$Location, "Microfauna")
# names(data2plot) = c("Sample", "PCoA1", "PCoA2", "Habitat Location", "Habitat", "Location", "Group")
# #get colours for plotting
# brewer.pal(n = 4, name = "Set1")
# #default legend order
# data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
# #plot
# mm.PCoA.p = ggplot(data=data2plot) +
#   geom_point(data=data2plot, aes(x=PCoA1, y=PCoA2, color=Location, shape=Habitat, fill= Location),
#              alpha=.5, size=4)+
#   scale_colour_manual(name = "",
#                       labels = c("Upper", "Lower"),
#                       values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   scale_fill_manual(name = "",
#                     labels = c("Upper", "Lower"),
#                     values = alpha(c(brewer.pal(n = 2, name = "Set1")), .6)) +
#   
#   scale_shape_manual(name = "",
#                      labels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"),
#                      values = c(22, 24, 25, 23))+
#   labs(x=labels[1], y=labels[2])+
#   theme_bw()+
#   theme(legend.position="bottom")+
#   ggtitle("Microfauna")
# 
# print(mm.PCoA.p)
# 
# #plot with patchwork
# PCoA.all.p1 = pro.PCoA.p / euk.PCoA.p / mm.PCoA.p
# PCoA.all.p1
# 
# pdf("../results/PCoA-plot-all-groups.pdf")
# print(PCoA.all.p1)
# dev.off()
# 
