
rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan) #anosim function
library(dendextend)
library(pairwiseAdonis)

# 2. Import data
ps.pro <- readRDS("../data/16S-phylo-object-rel.rds") 
ps.euk <- readRDS("../data/18S-phylo-object-rel.rds") 
#remove blanks
ps.euk = subset_samples(ps.euk, SampleType != "Control")
ps.euk = filter_taxa(ps.euk, function(x) sum(x) > 0, TRUE)

#remove the 2 most abundant taxa from ps.euk
#which ASVs have the most counts?
tax = data.frame(tax_table(ps.euk))
x = data.frame(sort(taxa_sums(ps.euk), decreasing = FALSE))
names(x) = "value"
y = x %>% 
  top_n(2, value)
z = tax[rownames(tax) %in% rownames(y),]
allTaxa = taxa_names(ps.euk)
allTaxa <- allTaxa[!(allTaxa %in% rownames(z))]
ps.euk.rm = prune_taxa(allTaxa, ps.euk)

#PERMANOVA function
# ps = ps.pro
# gene = "16S"

PERMANOVA_func <- function(ps, gene, m){
  
  set.seed(6666)
  
  #Generate distance matrix
  dist_matrix <- phyloseq::distance(ps, method = m) 
  
  #ADONIS test (a.k.a. PERMANOVA) for differences between ALL GROUPS 
  group.out = vegan::adonis2(dist_matrix ~ phyloseq::sample_data(ps)$SampleType)
  group.out
  #p > 0.05 no significant differences between groups
  #p < 0.05 significant differences between groups
  
  #ADONIS test can be confounded by differences in dispersion (spread of the data) so we need to test this
  #Dispersion test
  dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(ps)$SampleType)
  print(permutest(dispr))
  #p > 0.05 no significant differences between groups
  #p < 0.05 significant differences between groups
  
  #dispersion plots
  plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
  boxplot(dispr, main = "", xlab = "")
  
  #look for significant differences in the dissimilarities between PAIRS
  pair.out = pairwise.adonis(dist_matrix, phyloseq::sample_data(ps)$SampleType)
  
  #prep group and pair permanova results
  group.df = cbind(gene, group.out)
  pair.df = cbind(gene, pair.out)
  
  return(list(group.df, pair.df))
  
}


#Run function

#Pro
pro <- PERMANOVA_func(ps.pro, "16S", "wunifrac")
pro[[1]]
pro[[2]]

euk <- PERMANOVA_func(ps.euk, "18S", "wunifrac")
euk[[1]]
euk[[2]]

# #Pro
# pro <- PERMANOVA_func(ps.pro, "16S", "bray")
# pro[[1]]
# pro[[2]]
# 
# euk <- PERMANOVA_func(ps.euk, "18S", "bray")
# euk[[1]]
# euk[[2]]

#make group results into single dataframe
final.df = rbind(pro[[2]],
                euk[[2]])

#round our results
final.df$SumsOfSqs = round(final.df$SumsOfSqs, 3)
final.df$F.Model = round(final.df$F.Model, 3)
final.df$R2 = round(final.df$R2, 3)

write.csv(final.df, "../results/PERMANOVA.csv")
