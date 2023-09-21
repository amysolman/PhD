#1. Clear workspace and load packages


rm(list=ls())
graphics.off()

#install packages
library(readr)
library(seqinr)
library(tidyverse)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)

#install BiocManager packages
#install.packages("BiocManager", repos = "https://cloud.r-project.org")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("decontam")
#library(decontam)
#BiocManager::install("phyloseq")
library(phyloseq)


#2. Import data


#16S

#Count Table
pro_count_table <- read_tsv(file="../16S/table.tsv", skip=1)

#Specify the first column in our count table is the row names (ASV IDs)
pro_count_table <- column_to_rownames(pro_count_table, var = colnames(pro_count_table)[1])

#Rest of the data
pro_tax <- read_tsv(file="../16S/taxonomy.tsv")
pro_fasta <- read.fasta(file="../16S/dna-sequences.fasta")
pro_tree <- read_tree("../16S/tree.nwk")

#18S

#Count Table
euk_count_table <- read_tsv(file="../18S/table.tsv", skip=1)

#Specify the first column in our count table is the row names (ASV IDs)
euk_count_table <- column_to_rownames(euk_count_table, var = colnames(euk_count_table)[1])

#Rest of the data
euk_tax <- read_tsv(file="../18S/taxonomy.tsv")
euk_fasta <- read.fasta(file="../18S/dna-sequences.fasta")
euk_tree <- read_tree("../18S/tree.nwk")

#Metadata
metadata <- read.csv(file="../metadata.csv", sep=",")

#3. Wrangle metadata
#Sort metadata and count dataframes so they have the same sample names in the same order

#16S
#remove unnecessary data from the metadata file
pro_metadata = metadata[metadata$SampleID %in% names(pro_count_table),]
#make sure row has the right names
rownames(pro_metadata) = pro_metadata$SampleID

#18S
#remove unnecessary data from the metadata file
euk_metadata = metadata[metadata$SampleID %in% names(euk_count_table),]
#make sure row has the right names
rownames(euk_metadata) = euk_metadata$SampleID

#4. Sort taxonomy table so our output from QIIME classification step is useful for analysis

#16S
pro_taxa_tab <- pro_tax %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#18S
euk_taxa_tab <- euk_tax %>%
  mutate(euk_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=";$", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=" ", replacement="")) %>%
  separate(euk_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')


#5. Output number of sequences per sample in raw reads

#16S
pro_samp.counts = as.data.frame(colSums(pro_count_table))
pro_samp.counts$Sample = rownames(pro_samp.counts)
names(pro_samp.counts) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_samp.counts, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

#18S
euk_samp.counts = as.data.frame(colSums(euk_count_table))
euk_samp.counts$Sample = rownames(euk_samp.counts)
names(euk_samp.counts) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_samp.counts, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p


#6. Raw data basic stats: total reads per dataset, number of unique ASVs

#16S
sum(pro_samp.counts$Reads) #reads
nrow(pro_count_table) #ASVs

#18S
sum(euk_samp.counts$Reads) #reads
nrow(euk_count_table) #ASVs



#7. Convert data into phyloseq objects

#16S
pro_ASV = otu_table(as.matrix(pro_count_table), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(pro_metadata, row.names = rownames(pro_metadata)))
pro_TREE = pro_tree

#make into phyloseq object
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE)

#18S
euk_ASV = otu_table(as.matrix(euk_count_table), taxa_are_rows = TRUE)
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(euk_metadata, row.names = rownames(euk_metadata)))
euk_TREE = euk_tree

#make into phyloseq object
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE)


#8. Remove ambiguous annotations and singletons

#for 16S data we remove those assigned to chloroplast at the order level 
#and mitochondria at the family level because bacteria don't have them! 
#they just have DNA like prokaryotes but they're from the organelles of eukaryotes
#remove vertebrata (animals) and embryophyta (plants) from 18S data 

#16S
pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukayote") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) > 1, pro_ps.ambig) 

#18S dataset with micrometazoan's removed
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata") & !Class %in% c("Embryophyta"))
euk_ps.ambig <- prune_taxa(taxa_sums(euk_ps.ambig) > 1, euk_ps.ambig) 


#11. Re-root tree. Use this function > (https://john-quensen.com/r/unifrac-and-tree-roots/) < picks longest branch

# First define the function from link above to find furthest outgroup
pick_new_outgroup <- function(tree.unrooted){
  #tablif parts of the tree that we need
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  #Take out the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}


#16S
pro_tree_new <- phy_tree(pro_ps.ambig)
pro_new_outgroup <- pick_new_outgroup(pro_tree_new)
# Re-root tree
pro_new_tree_root <- ape::root(pro_tree_new, outgroup=pro_new_outgroup, resolve.root=TRUE)
# Convert to dichotomy tree
pro_new_tree_dich <- ape::multi2di(pro_new_tree_root)
phy_tree(pro_ps.ambig) <- pro_new_tree_dich

#18S
euk_tree_new <- phy_tree(euk_ps.ambig)
euk_new_outgroup <- pick_new_outgroup(euk_tree_new)
# Re-root tree
euk_new_tree_root <- ape::root(euk_tree_new, outgroup=euk_new_outgroup, resolve.root=TRUE)
# Convert to dichotomy tree
euk_new_tree_dich <- ape::multi2di(euk_new_tree_root)
phy_tree(euk_ps.ambig) <- euk_new_tree_dich

#12. Filtered data basic stats: total reads per dataset, number of unique ASVs

#16S
sum(colSums(data.frame(otu_table(pro_ps.ambig)))) #reads
nrow(data.frame(otu_table(pro_ps.ambig))) #ASVs

#18S
sum(colSums(data.frame(otu_table(euk_ps.ambig)))) #reads
nrow(data.frame(otu_table(euk_ps.ambig))) #ASVs

#Normalise by relative abundance
pro.rel = transform_sample_counts(pro_ps.ambig, function(x) x/sum(x))
colSums(data.frame(otu_table(pro.rel)))

euk.rel = transform_sample_counts(euk_ps.ambig, function(x) x/sum(x))
colSums(data.frame(otu_table(euk.rel)))

#Normalise by rarefaction
pro.rar = rarefy_even_depth(pro_ps.ambig, replace=TRUE)
colSums(data.frame(otu_table(pro.rar)))
pro.rar

euk.rar = rarefy_even_depth(euk_ps.ambig, replace=TRUE)
colSums(data.frame(otu_table(euk.rar)))
euk.rar


#13. Save phyloseq objects

#16S
saveRDS(pro_ps.ambig, "../data/16S-phylo-object.rds")
saveRDS(pro.rel, "../data/16S-phylo-object-rel.rds")
saveRDS(pro.rar, "../data/16S-phylo-object-rar.rds")

#18S
saveRDS(euk_ps.ambig, "../data/18S-phylo-object.rds")
saveRDS(euk.rel, "../data/18S-phylo-object-rel.rds")
saveRDS(euk.rar, "../data/18S-phylo-object-rar.rds")
