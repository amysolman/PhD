
# Taken from paper: https://link.springer.com/article/10.1007/s00248-020-01616-4
# Permutational multivariate analysis of variance (PERMANOVA) using the adonis function and non-metric multidimensional scaling (NMDS) using the metaMDS function were applied from the vegan package and calculated with Bray-Curtis distance. In particular for PERMANOVA, the pairwise multiple comparison (post hoc) was further carried out with Bonferroni method in the pairwisea.adonis function. 
# 
# Code developed from: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#   

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
# Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#Prokaryote
#set habitat location variable
sample_data(ps.pro)$Habitat_location <- paste0(sample_data(ps.pro)$Habitat, " (", sample_data(ps.pro)$Location, " Foxfonna)")
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.pro, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.pro)$Habitat_location) #significant differences between groups

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Habitat_location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
pro.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Habitat_location)
pro.pair.out$Group = "Prokaryote"
write.csv(pro.pair.out, "../results/PERMANOVA-prokaryotes.csv")

#Eukaryote
#set habitat location variable
sample_data(ps.euk)$Habitat_location <- paste0(sample_data(ps.euk)$Habitat, " (", sample_data(ps.euk)$Location, " Foxfonna)")
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.euk, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.euk)$Habitat_location) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Habitat_location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
euk.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Habitat_location)
euk.pair.out$Group = "Eukaryote"
write.csv(euk.pair.out, "../results/PERMANOVA-eukaryotes.csv")

#Microfauna
#set habitat location variable
sample_data(ps.mm)$Habitat_location <- paste0(sample_data(ps.mm)$Habitat, " (", sample_data(ps.mm)$Location, " Foxfonna)")
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.mm, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.mm)$Habitat_location) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Habitat_location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
mm.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Habitat_location)
mm.pair.out$Group = "Microfauna"
write.csv(mm.pair.out, "../results/PERMANOVA-microfauna.csv")

#merge into one dataframe and round results

#combine location PERMANOVA
full.loc_hab.perma = rbind(pro.pair.out, euk.pair.out, mm.pair.out)

#round
full.loc_hab.perma$SumsOfSqs = round(full.loc_hab.perma$SumsOfSqs, 2)
full.loc_hab.perma$F.Model = round(full.loc_hab.perma$F.Model, 2)
full.loc_hab.perma$R2 = round(full.loc_hab.perma$R2, 2)

#export the results
write.csv(full.loc_hab.perma, "../results/PERMANOVA-final-loc-and-hab.csv")

################DIFFERENCES BETWEEN HABITATS##########################

#Prokaryote
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.pro, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.pro)$Habitat) #significant differences between groups

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Habitat)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
pro.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Habitat)
pro.pair.out$Group = "Prokaryote"

#Eukaryote
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.euk, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.euk)$Habitat) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Habitat)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
euk.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Habitat)
euk.pair.out$Group = "Eukaryote"

#Microfauna
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.mm, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.mm)$Habitat) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Habitat)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates no significant difference in dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
mm.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Habitat)
mm.pair.out$Group = "Microfauna"

#combine habitat PERMANOVA
full.hab.perma = rbind(pro.pair.out, euk.pair.out, mm.pair.out)


#############DIFFERENCES BETWEEN UPPER AND LOWER FOXFONNA###############

#Prokaryote
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.pro, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.pro)$Location) #significant differences between groups

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
pro.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.pro)$Location)
pro.pair.out$Group = "Prokaryote"

#Eukaryote
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.euk, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.euk)$Location) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates significantly different dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
euk.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.euk)$Location)
euk.pair.out$Group = "Eukaryote"

#Microfauna
#Generate distance matrix
bray_dist_matrix <- phyloseq::distance(ps.mm, method = "bray") 
#ADONIS test for differences between our groups 
vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(ps.mm)$Location) #significant differences between groups
#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
#Dispersion test and plot
dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr) #indicates no significant difference in dispersion between the groups

#look for significant differences in the dissimilarities of our habitats/locations
mm.pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(ps.mm)$Location)
mm.pair.out$Group = "Microfauna"

#combine location PERMANOVA
full.loc.perma = rbind(pro.pair.out, euk.pair.out, mm.pair.out)

#combine all results and round
final.res = rbind(full.hab.perma, full.loc.perma)

#round
final.res$SumsOfSqs = round(final.res$SumsOfSqs, 2)
final.res$F.Model = round(final.res$F.Model, 2)
final.res$R2 = round(final.res$R2, 2)

#export the results
write.csv(final.res, "../results/PERMANOVA-final.csv")
