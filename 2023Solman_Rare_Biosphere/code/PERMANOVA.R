
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

pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds")


# We can test for differences between group of communities using two different methods, ANOSIM and PERMANOVA.
# ANOSIM - tests whether distances between groups (between abundant and rare taxa) are greater than distances within groups.
# PERMANOVA - test whether distances differ between groups.
# 
# A test by Anderson and Walsh (2013) found PERMANOVA was more rebust when using ecological data but the assumptions of each method must be tested.
# 
# MY QUESTIONS:
#   IS THE PHYLOGENETIC STRUCTURE DIFFERENT BETWEEN POLES FOR EACH SUBCOMMUNITY? (E.g. are abundant groups different between the poles - and how different?)
# IS THE PHYLOGENETIC STRUCTURE DIFFERENT BETWEEN EACH SUBCOMMUNITY? (E.g. are abundant groups phylogenetically different from intermediate and rare groups - and how different?)
# 
# PERMANOVA
# 
# Taken from paper: https://link.springer.com/article/10.1007/s00248-020-01616-4
# Permutational multivariate analysis of variance (PERMANOVA) using the adonis function and non-metric multidimensional scaling (NMDS) using the metaMDS function were applied from the vegan package and calculated with UniFrac distance. In particular for PERMANOVA, the pairwise multiple comparison (post hoc) was further carried out with Bonferroni method in the pairwisea.adonis function. 

# Code developed from: https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
  
###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN ABUNDANT PROKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps1 = subset_samples(pro, Subcommunity == "Abundant")

#remove taxa or samples with zero counts
ps1 = filter_taxa(ps1, function(x) sum(x) >= 1, TRUE)
ps1 = prune_samples(sample_sums(ps1)>=1, ps1)

#transform data by relative abundance
ps1.rel = transform_sample_counts(ps1, function(x) x/sum(x))
colSums(data.frame(otu_table(ps1.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps1.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps1.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps1.rel)$Pole)
permutest(dispr) #indicates significantly different dispersion between the groups

#Should we not use this test because the dispersions are different?
# testing null hypothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#in this example Pole explained little of the difference in dispersion (0.33) compared to the residual variation (2.25) so I won't worry about differences in dispersion too much

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
pro.abun.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps1.rel)$Pole)

write.csv(pro.abun.pair.out, "../results/PERMANOVA-abundant-prokaryotes.csv")



###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN INTERMEDIATE PROKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps2 = subset_samples(pro, Subcommunity == "Intermediate")

#remove taxa or samples with zero counts
ps2 = filter_taxa(ps2, function(x) sum(x) >= 1, TRUE)
ps2 = prune_samples(sample_sums(ps2)>=1, ps2)

#transform data by relative abundance
ps2.rel = transform_sample_counts(ps2, function(x) x/sum(x))
colSums(data.frame(otu_table(ps2.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps2.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps2.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps2.rel)$Pole)
permutest(dispr) #indicates no significantly different dispersion between the groups

# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor


#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
pro.int.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps2.rel)$Pole)

write.csv(pro.int.pair.out, "../results/PERMANOVA-intermediate-prokaryotes.csv")

###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN RARE PROKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps3 = subset_samples(pro, Subcommunity == "Rare")

#remove taxa or samples with zero counts
ps3 = filter_taxa(ps3, function(x) sum(x) >= 1, TRUE)
ps3 = prune_samples(sample_sums(ps3)>=1, ps3)

#transform data by relative abundance
ps3.rel = transform_sample_counts(ps3, function(x) x/sum(x))
colSums(data.frame(otu_table(ps3.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps3.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps3.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps3.rel)$Pole)
permutest(dispr) #indicates no significantly different dispersion between the groups

# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
pro.rare.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps3.rel)$Pole)

write.csv(pro.rare.pair.out, "../results/PERMANOVA-rare-prokaryotes.csv")

###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN ABUNDANT EUKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps4 = subset_samples(euk, Subcommunity == "Abundant")

#remove taxa or samples with zero counts
ps4 = filter_taxa(ps4, function(x) sum(x) >= 1, TRUE)
ps4 = prune_samples(sample_sums(ps4)>=1, ps4)

#transform data by relative abundance
ps4.rel = transform_sample_counts(ps4, function(x) x/sum(x))
colSums(data.frame(otu_table(ps4.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps4.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps4.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps4.rel)$Pole)
permutest(dispr) #indicates no significantly different dispersion between the groups
# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
euk.abun.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps4.rel)$Pole)

write.csv(euk.abun.pair.out, "../results/PERMANOVA-abundant-eukaryotes.csv")


###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN INTERMEDIATE EUKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps5 = subset_samples(euk, Subcommunity == "Intermediate")

#remove taxa or samples with zero counts
ps5 = filter_taxa(ps5, function(x) sum(x) >= 1, TRUE)
ps5 = prune_samples(sample_sums(ps5)>=1, ps5)

#transform data by relative abundance
ps5.rel = transform_sample_counts(ps5, function(x) x/sum(x))
colSums(data.frame(otu_table(ps5.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps5.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps5.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps5.rel)$Pole)
permutest(dispr) #indicates no significantly different dispersion between the groups

# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor


#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
euk.int.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps5.rel)$Pole)

write.csv(euk.int.pair.out, "../results/PERMANOVA-intermediate-eukaryotes.csv")

###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN RARE EUKARYOTES COMMUNITIES FROM THE ARCTIC AND ANTARCTIC?

#subset to abundant coms only
ps6 = subset_samples(euk, Subcommunity == "Rare")

#remove taxa or samples with zero counts
ps6 = filter_taxa(ps6, function(x) sum(x) >= 1, TRUE)
ps6 = prune_samples(sample_sums(ps6)>=1, ps6)

#transform data by relative abundance
ps6.rel = transform_sample_counts(ps6, function(x) x/sum(x))
colSums(data.frame(otu_table(ps6.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps6.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps6.rel)$Pole) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps6.rel)$Pole)
permutest(dispr) #indicates no significantly different dispersion between the groups

# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our habitats/locations
euk.rare.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps6.rel)$Pole)

write.csv(euk.rare.pair.out, "../results/PERMANOVA-rare-eukaryotes.csv")

###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN PROKARYOTE SUBCOMMUNITIES?

#transform data by relative abundance
ps_pro.rel = transform_sample_counts(pro, function(x) x/sum(x))
colSums(data.frame(otu_table(ps_pro.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps_pro.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_pro.rel)$Subcommunity) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_pro.rel)$Subcommunity)
permutest(dispr) #indicates significantly different dispersion between the groups
# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#in this example Subcommunity explained little of the difference in dispersion (1.37) compared to the residual variation (2.88) so I won't worry about differences in dispertion too much

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our subcommunities
pro.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps_pro.rel)$Subcommunity)

write.csv(pro.pair.out, "../results/PERMANOVA-prokaryotes.csv")


###################################################################################################
###################################################################################################

#ARE THERE DIFFERENCES BETWEEN EUKARYOTE SUBCOMMUNITIES?

#remove taxa or samples with zero counts
ps_euk = filter_taxa(euk, function(x) sum(x) >= 1, TRUE)
ps_euk = prune_samples(sample_sums(ps_euk)>=1, ps_euk)

#transform data by relative abundance
ps_euk.rel = transform_sample_counts(ps_euk, function(x) x/sum(x))
colSums(data.frame(otu_table(ps_euk.rel))) == 1 #should all be TRUE

#Generate distance matrix
dist_matrix <- phyloseq::distance(ps_euk.rel, method = "bray") 

#ADONIS test for differences between our groups (Poles)
vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps_euk.rel)$Subcommunity) 
#indicates significant differences between groups P < 0.05

#ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this

#Dispersion test and plot
dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps_euk.rel)$Subcommunity)
permutest(dispr) #indicates significantly different dispersion between the groups

# testing null hyothesis that means are the same against alternative hypothesis that there is a difference.
#we know what our p value means but what about sum of squares (SS) and mean sum of squares (MSS)?
#MSS = SS / degrees of freedom (number of data points 234 minus the number of factors 2 for Arctic and Antarctic).
#Total SS = square root of the difference between the datapoint and the overall mean value (how different are our data points from the mean, how much variation is there in the data)
#Residual SS = square root of the difference between the datapoint and the mean for that Pole (how different are our data points from the Pole mean, how much variation is there in the data)
#Pole SS = total SS - residual ss = variation explained by the pole factor

#in this example Subcommunity explained little of the difference in dispersion (0.81) compared to the residual variation (2.45) so I won't worry about differences in dispertion too much

#plot 
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")

#look for significant differences in the dissimilarities of our subcommunities
euk.pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps_euk.rel)$Subcommunity)

write.csv(euk.pair.out, "../results/PERMANOVA-eukaryotes.csv")

###################################################################################################
###################################################################################################

#Combine all tables

res = rbind(pro.abun.pair.out, pro.int.pair.out, pro.rare.pair.out, 
      euk.abun.pair.out, euk.int.pair.out, euk.rare.pair.out,
      pro.pair.out, euk.pair.out)
group = c("Prokaryotes", "Prokaryotes", "Prokaryotes", 
          "Eukaryotes", "Eukaryotes", "Eukaryotes", 
          "Prokaryotes", "Prokaryotes", "Prokaryotes", 
          "Eukaryotes", "Eukaryotes", "Eukaryotes")
sub = c("Abundant", "Intermediate", "Rare",
        "Abundant", "Intermediate", "Rare",
        NA, NA, NA, NA, NA, NA)
res2 = cbind(group, sub, res)

#round our results
res2$SumsOfSqs = round(res2$SumsOfSqs, 3)
res2$F.Model = round(res2$F.Model, 3)
res2$R2 = round(res2$R2, 3)

#export table
write.csv(res2, "../results/PERMANOVA.csv")
