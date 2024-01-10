
rm(list=ls())
library(phyloseq)
library(vegan) #anosim function
library(dendextend)
library(pairwiseAdonis)

# 2. Import data
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 

#PERMANOVA function
phylo = pro
habitat = gsub(" ", "", data.frame(sample_data(subset_samples(pro, Habitat == "Snow")))$R_Names)
gene = "16S"
hab_nam = "Snow"

PERMANOVA_func <- function(phylo, habitat, gene, hab_nam){
  
  set.seed(6666)
  
  #subset to habitat
  sub <- prune_samples(habitat, phylo)
  
  #Generate distance matrix
  bray_dist_matrix <- phyloseq::distance(sub, method = "bray") 
  
  #ADONIS test (a.k.a. PERMANOVA) for differences between ALL GROUPS 
  group.out = vegan::adonis2(bray_dist_matrix ~ phyloseq::sample_data(sub)$Treatment)
  group.out
  #p > 0.05 no significant differences between groups
  #p < 0.05 significant differences between groups
  
  #ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
  #Dispersion test
  dispr <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(sub)$Treatment)
  print(permutest(dispr))
  #p > 0.05 no significant differences between groups
  #p < 0.05 significant differences between groups
  
  #dispersion plots
  plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
  boxplot(dispr, main = "", xlab = "")
  
  #look for significant differences in the dissimilarities between PAIRS
  pair.out = pairwise.adonis(bray_dist_matrix, phyloseq::sample_data(sub)$Treatment)
  
  #prep group and pair permanova results
  group.df = cbind(gene, hab_nam, group.out)
  pair.df = cbind(gene, hab_nam, pair.out)
  
  return(list(group.df, pair.df))
  
}


#Run function

#Pro snow
pro.perm.snow <- PERMANOVA_func(pro, 
                                gsub(" ", "", data.frame(sample_data(subset_samples(pro, Habitat == "Snow")))$R_Names),
                                "16S", "Snow")
#Pro spring ice
pro.perm.sp <- PERMANOVA_func(pro, 
                              gsub(" ", "", data.frame(sample_data(subset_samples(pro, Habitat == "Spring Ice")))$R_Names),
                              "16S", "Spring Ice")
#Pro summer ice
pro.perm.sum <- PERMANOVA_func(pro, 
                               gsub(" ", "", data.frame(sample_data(subset_samples(pro, Habitat == "Summer Ice")))$R_Names),
                               "16S", "Summer Ice")
#Pro cryoconite
pro.perm.cry <- PERMANOVA_func(pro, 
                               gsub(" ", "", data.frame(sample_data(subset_samples(pro, Habitat == "Cryoconite")))$R_Names),
                               "16S", "Cryoconite")
#euk snow
euk.perm.snow <- PERMANOVA_func(euk, 
                                gsub(" ", "", data.frame(sample_data(subset_samples(euk, Habitat == "Snow")))$R_Names),
                                "18S", "Snow")
#euk spring ice
euk.perm.sp <- PERMANOVA_func(euk, 
                              gsub(" ", "", data.frame(sample_data(subset_samples(euk, Habitat == "Spring Ice")))$R_Names),
                              "18S", "Spring Ice")
#euk summer ice
euk.perm.sum <- PERMANOVA_func(euk, 
                               gsub(" ", "", data.frame(sample_data(subset_samples(euk, Habitat == "Summer Ice")))$R_Names),
                               "18S", "Summer Ice")
#euk cryoconite
euk.perm.cry <- PERMANOVA_func(euk, 
                               gsub(" ", "", data.frame(sample_data(subset_samples(euk, Habitat == "Cryoconite")))$R_Names),
                               "18S", "Cryoconite")

#make group results into single dataframe
final.df = rbind(pro.perm.snow[[1]],
                 pro.perm.sp[[1]],
                 pro.perm.sum[[1]],
                 pro.perm.cry[[1]],
                 euk.perm.snow[[1]],
                 euk.perm.sp[[1]],
                 euk.perm.sum[[1]],
                 euk.perm.cry[[1]])

#round our results
final.df$SumOfSqs = round(final.df$SumOfSqs, 3)
final.df$F = round(final.df$F, 3)
final.df$R2 = round(final.df$R2, 3)
final.df = cbind(rownames(final.df), final.df)

write.csv(final.df, paste("../results/pma-PERMANOVA.csv"))
