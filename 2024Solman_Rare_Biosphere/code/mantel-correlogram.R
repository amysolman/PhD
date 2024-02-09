rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan) #wascores
library(ggplot2)
library(picante)
library(ecodist) #for distance() function
library(parallel)
library(svglite)
library(dplyr) #for %>%
library(scales) #for percentages

# Step 2: Read in data
pro <- readRDS("../results/16S-phylo-object-rarefied-var-trans.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied-var-trans.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) > 0, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) > 0, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.int <- prune_samples(sample_sums(euk.ant.int)>0, euk.ant.int)
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")

#phylo = pro.ant

mantel_correlogram_func <- function(phylo){
  
  #PROPORTIONAL TRANSFORMATION
  ps.rel <- transform_sample_counts(phylo, function(x) x/sum(x))
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(ps.rel))
  
  #list of variables we're interested in 
  keeps1 <- c("WaterDepth", "SedimentDepth", "TotalDepth", "pH", "DOC", "Cl", "SO4", "Mg","Ca", "Area")
  
  keeps2 <-c("DistanceToSea", "Elevation", "HCO3")
  
  samp_data1 <- samp_data[ , (names(samp_data) %in% keeps1)]
  samp_data2 <- samp_data[ , (names(samp_data) %in% keeps2)]
  
  #remove columns with more than 50% missing variables
  samp_data_trim1 <- samp_data1[ lapply( samp_data1, function(x) sum(is.na(x)) / length(x) ) < 0.5 ]
  samp_data_trim2 <- samp_data2[ lapply( samp_data2, function(x) sum(is.na(x)) / length(x) ) < 0.7 ]
  
  #Only keep complete cases
  meta1 <- samp_data_trim1[complete.cases(samp_data_trim1),]
  meta2 <- samp_data_trim2[complete.cases(samp_data_trim2),]
  
  #MANTEL CORRELOGRAM ONE
  ps.smol1 = prune_samples(rownames(meta1), ps.rel)
  # sample_df1 <- data.frame(sample_data(trim1))
  meta.trim1 <- meta1
  
  #remove ASVs with 0 relative abundance
  trim.filter1 <- filter_taxa(ps.smol1, function(x) sum(x) > 0, TRUE)
  
  #Get phylogenetic distances
  tree1 <- phy_tree(trim.filter1) #pull out phylogenetic tree
  tree.dist1 <- cophenetic(tree1)
  asvs1 <- tree1$tip.label
  phylo.dists1 <- tree.dist1[asvs1, asvs1]
  phylo.dists1[upper.tri(phylo.dists1, diag=TRUE)] = NA
  
  #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  #turn the dataframe into a matrix
  niches1 <- as.matrix(meta.trim1) 
  
  asv.table1 <- data.frame(t(otu_table(trim.filter1)), check.names = FALSE)
  asv.niche1 <- wascores(niches1, asv.table1) #this will find the weighted mean environmental parameter for each asv
  asv.niche.df1 <- data.frame(asv.niche1, check.names = FALSE)
  
  #generate euclidean distance matrix for each ASV using combined environmental parameters
  dist.out1 = as.matrix(dist(asv.niche.df1), labels=TRUE)
  
  #calculate mantel correlogram
  corlg1 <- mantel.correlog(dist.out1, phylo.dists1, r.type="spearman")
  
  #MANTEL CORRELOGRAM TWO
  ps.smol2 = prune_samples(rownames(meta2), ps.rel)
  # sample_df1 <- data.frame(sample_data(trim1))
  meta.trim2 <- meta2
  
  #remove ASVs with 0 relative abundance
  trim.filter2 <- filter_taxa(ps.smol2, function(x) sum(x) > 0, TRUE)
  
  #Get phylogenetic distances
  tree2 <- phy_tree(trim.filter2) #pull out phylogenetic tree
  tree.dist2 <- cophenetic(tree2)
  asvs2 <- tree2$tip.label
  phylo.dists2 <- tree.dist2[asvs2, asvs2]
  phylo.dists2[upper.tri(phylo.dists2, diag=TRUE)] = NA
  
  #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  #turn the dataframe into a matrix
  niches2 <- as.matrix(meta.trim2) 
  
  asv.table2 <- data.frame(t(otu_table(trim.filter2)), check.names = FALSE)
  asv.niche2 <- wascores(niches2, asv.table2) #this will find the weighted mean environmental parameter for each asv
  asv.niche.df2 <- data.frame(asv.niche2, check.names = FALSE)
  
  #generate euclidean distance matrix for each ASV using combined environmental parameters
  dist.out2 = as.matrix(dist(asv.niche.df2), labels=TRUE)
  
  #calculate mantel correlogram
  corlg2 <- mantel.correlog(dist.out2, phylo.dists2, r.type="spearman")
  
  # Prep data
  crlg1 <- data.frame(corlg1$mantel.res)
  crlg2 <- data.frame(corlg2$mantel.res)
  
  crlg1 <- crlg1 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  
  crlg2 <- crlg2 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  
  crlg1$D.cl = rownames(crlg1)
  crlg1$Analysis = "1"
  crlg2$D.cl = rownames(crlg2)
  crlg2$Analysis = "2"
  
  final.crlg = rbind(crlg1, crlg2)
  
  
  return(final.crlg)
}

#not going to run the function because we did this on the HPC
# pro.ant.mc <- mantel_correlogram_func(pro.ant)
# pro.ant.abun.mc <- mantel_correlogram_func(pro.ant.abun)
# pro.ant.int.mc <- mantel_correlogram_func(pro.ant.int)
# pro.ant.rare.mc <- mantel_correlogram_func(pro.ant.rare)
# 
# euk.ant.mc <- mantel_correlogram_func(euk.ant)
# euk.ant.abun.mc <- mantel_correlogram_func(euk.ant.abun)
# euk.ant.int.mc <- mantel_correlogram_func(euk.ant.int)
# euk.ant.rare.mc <- mantel_correlogram_func(euk.ant.rare)
# 
# #combine into single dataframe
# pro.ant.mc$Group = "Prokaryote"
# pro.ant.mc$Subcommunity = "Full"
# pro.ant.abun.mc$Group = "Prokaryote"
# pro.ant.abun.mc$Subcommunity = "Abundant"
# pro.ant.int.mc$Group = "Prokaryote"
# pro.ant.int.mc$Subcommunity = "Intermediate"
# pro.ant.rare.mc$Group = "Prokaryote"
# pro.ant.rare.mc$Subcommunity = "Rare"
# euk.ant.mc$Group = "Eukaryote"
# euk.ant.mc$Subcommunity = "Full"
# euk.ant.abun.mc$Group = "Eukaryote"
# euk.ant.abun.mc$Subcommunity = "Abundant"
# euk.ant.int.mc$Group = "Eukaryote"
# euk.ant.int.mc$Subcommunity = "Intermediate"
# euk.ant.rare.mc$Group = "Eukaryote"
# euk.ant.rare.mc$Subcommunity = "Rare"
# 
# plot.df = rbind(pro.ant.mc, pro.ant.abun.mc, pro.ant.int.mc, pro.ant.rare.mc, euk.ant.mc, euk.ant.abun.mc, euk.ant.int.mc, euk.ant.rare.mc)
# 
# #save the results dataframe
# write.csv(plot.df, "../results/mantel-correlogram-results.csv")

#Read in the correlogran results from the HPC
plot.df = read.csv("../results/mantel-correlogram-results.csv")

plot.df$Analysis = as.factor(plot.df$Analysis)

#order our groups
df2plot = plot.df[plot.df$Subcommunity != "Full",]

df2plot$Group = factor(df2plot$Group, levels=c("Prokaryote", "Eukaryote"))
df2plot$Subcommunity = factor(df2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

p = ggplot(data=df2plot, aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis)) +
  geom_point(data=df2plot[df2plot$sig=="significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black", size=3, shape=16)+
  geom_point(data=df2plot[df2plot$sig=="non-significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black",size=3, shape=1)+
  geom_line(size=1)+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x = "Phylogenetic distance class", y="Mantel correlation")+
  # ylim(-0.15, 0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=7),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=10),
        legend.position = "bottom",
        legend.text = element_text(size=10))+
  facet_grid(rows=vars(Group), cols=vars(Subcommunity), scales = "free_x")

print(p)

pdf("../results/mantel-correlogram.pdf", width=10, height=5)
print(p)
dev.off()

#full community only
df2plot = plot.df[plot.df$Subcommunity == "Full",]

df2plot$Group = factor(df2plot$Group, levels=c("Prokaryote", "Eukaryote"))
df2plot$Subcommunity = factor(df2plot$Subcommunity, levels = c("Full"))

p2 = ggplot(data=df2plot, aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis)) +
  geom_point(data=df2plot[df2plot$sig=="significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black", size=3, shape=16)+
  geom_point(data=df2plot[df2plot$sig=="non-significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black",size=3, shape=1)+
  geom_line(size=1)+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x = "Phylogenetic distance class", y="Mantel correlation")+
  # ylim(-0.15, 0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=7),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=10),
        legend.position = "bottom",
        legend.text = element_text(size=10))+
  facet_grid(rows=vars(Group), cols=vars(Subcommunity), scales = "free_x")

print(p2)

pdf("../results/mantel-correlogram-full.pdf", width=10, height=5)
print(p2)
dev.off()