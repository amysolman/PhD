
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
library(glue)
library(dendextend)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#Code developed from: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#Prokaryotes
#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps.pro))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#set habitat location variable
sample_data(ps.pro)$Habitat_location <- paste0(sample_data(ps.pro)$Habitat, " (", sample_data(ps.pro)$Location, " Foxfonna)")
#get colours for plotting
brewer.pal(n = 8, name = "Paired")
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps.pro))
#colorCode <- c(`Upper` = "red", `Lower` = "blue")
# colorCode <- c(`Snow (Lower Foxfonna)` = "#A6CEE3", 
#                `Snow (Upper Foxfonna)`= "#1F78B4", 
#                `Spring Ice (Lower Foxfonna)` =  "#B2DF8A",
#                `Spring Ice (Upper Foxfonna)` =  "#33A02C",
#                `Summer Ice (Lower Foxfonna)` = "#FB9A99",
#                `Summer Ice (Upper Foxfonna)` = "#E31A1C",
#                `Cryoconite (Lower Foxfonna)` = "#FDBF6F",
#                `Cryoconite (Upper Foxfonna)` = "#FF7F00")
colorCode <- c(`Snow`= "#1F78B4", 
               `Spring Ice` =  "#33A02C",
               `Summer Ice` = "#E31A1C",
               `Cryoconite` = "#FF7F00")
# labels_colors(ward) <- colorCode[meta$Location][order.dendrogram(ward)]
# labels_colors(ward) <- colorCode[meta$Habitat_location][order.dendrogram(ward)]
labels_colors(ward) <- colorCode[meta$Habitat][order.dendrogram(ward)]
#Plot
plot(ward)

#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps.pro))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
#Save as dendrogram
ward1 <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#set habitat location variable
sample_data(ps.pro)$Habitat_location <- paste0(sample_data(ps.pro)$Habitat, " (", sample_data(ps.pro)$Location, " Foxfonna)")
#get colours for plotting
#brewer.pal(n = 8, name = "Paired")
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps.pro))
colorCode <- c(`Upper` = "red", `Lower` = "blue")
labels_colors(ward1) <- colorCode[meta$Location][order.dendrogram(ward1)]
#Plot
pdf("../results/hierachical-clustering-prokaryotes.pdf", width=14, height=6)
plot(ward1)
dev.off()

plot(ward)


####EUKARYOTES

#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps.euk))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
#Save as dendrogram
ward2 <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#set habitat location variable
sample_data(ps.euk)$Habitat_location <- paste0(sample_data(ps.euk)$Habitat, " (", sample_data(ps.euk)$Location, " Foxfonna)")
#get colours for plotting
#brewer.pal(n = 8, name = "Paired")
#provide color codes
meta <- data.frame(phyloseq::sample_data(ps.euk))
colorCode <- c(`Upper` = "red", `Lower` = "blue")
labels_colors(ward2) <- colorCode[meta$Location][order.dendrogram(ward2)]
#Plot
pdf("../results/hierachical-clustering-eukaryotes.pdf", width=14, height=6)
plot(ward2)
dev.off()


####MICROFAUNA

#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps.mm))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
#Save as dendrogram
ward3 <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#set habitat location variable
sample_data(ps.mm)$Habitat_location <- paste0(sample_data(ps.mm)$Habitat, " (", sample_data(ps.mm)$Location, " Foxfonna)")
#get colours for plotting
#brewer.pal(n = 8, name = "Paired")
#provide color codes
meta <- data.frame(phyloseq::sample_data(ps.mm))
colorCode <- c(`Upper` = "red", `Lower` = "blue")
labels_colors(ward3) <- colorCode[meta$Location][order.dendrogram(ward3)]
#Plot
pdf("../results/hierachical-clustering-microfauna.pdf", width=14, height=6)
plot(ward3)
dev.off()

final.p = plot(ward1) / plot(ward2)
final.p
