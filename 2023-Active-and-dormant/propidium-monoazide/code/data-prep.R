# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Sort metadata
# 4. Sort taxonomy table
# 5. Output number of sequences per sample in raw reads
# 6. Raw data basic stats: total reads per dataset, number of unique ASVs
# 7. Convert data into phyloseq objects
# 8. Remove ambiguous annotations and singletons
# 9. Re-root tree 
# 10. Normalise data
# 11. Save phyloseq objects


#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)

#2. Import data
# #16S

#Count Table
pro_count_table <- read_tsv(file="../qiime-work-raw/16S/export/table.tsv", skip=1)
#Specify the first column in our count table is the row names (ASV IDs)
pro_count_table <- column_to_rownames(pro_count_table, var = colnames(pro_count_table)[1])
pro_sample_names = as.numeric(names(pro_count_table))

#Rest of the data
pro_tax <- read_tsv(file="../qiime-work-raw/16S/export/taxonomy.tsv")
pro_fasta <- read.fasta(file="../qiime-work-raw/16S/export/dna-sequences.fasta")
pro_tree <- read_tree("../qiime-work-raw/16S/export/tree.nwk")

#18S

#Count Table
euk_count_table <- read_tsv(file="../qiime-work-raw/18S/export/table.tsv", skip=1)
#Specify the first column in our count table is the row names (ASV IDs)
euk_count_table <- column_to_rownames(euk_count_table, var = colnames(euk_count_table)[1])
#sample_names = as.numeric(names(count_table))

#Rest of the data
euk_tax <- read_tsv(file="../qiime-work-raw/18S/export/taxonomy.tsv")
euk_fasta <- read.fasta(file="../qiime-work-raw/18S/export/dna-sequences.fasta")
euk_tree <- read_tree("../qiime-work-raw/18S/export/tree.nwk")

#Metadata
metadata <- read.csv(file="../data/metadata.csv", sep=",")


#3. Sort metadata

#Sort metadata and count dataframes so they have the same sample names in the same order
pro_metadata = metadata[metadata$Gene == "16S",]
# #make sure row has the right names
rownames(pro_metadata) = pro_metadata$R_Names

#remove unneccesary data from the metadata file
euk_metadata = metadata[metadata$Gene == "18S",]
# #make sure row has the right names
rownames(euk_metadata) = euk_metadata$R_Names


#4. Sort taxonomy table

#Now to modify the taxonomic table for my data
pro_taxa_tab <- pro_tax %>%
  #mutate(pro_tax=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#Now to modify the taxonomic table for my data
euk_taxa_tab <- euk_tax %>%
  #mutate(euk_tax=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=";$", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=" ", replacement="")) %>%
  separate(euk_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')


#5. Output number of sequences per sample in raw reads

pro_samp.counts = as.data.frame(colSums(pro_count_table))
pro_samp.counts$Sample = rownames(pro_samp.counts)
names(pro_samp.counts) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_samp.counts, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

pdf("../results/prokaryote-count-histogram.pdf", width=20, height=10)
print(pro_p)
dev.off()

#5. Output number of sequences per sample in raw reads

euk_samp.counts = as.data.frame(colSums(euk_count_table))
euk_samp.counts$Sample = rownames(euk_samp.counts)
names(euk_samp.counts) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_samp.counts, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p

pdf("../results/eukaryote-count-histogram.pdf", width=20, height=10)
print(euk_p)
dev.off()

#6. Raw data basic stats: total reads per dataset, number of unique ASVs

#how many 16S reads
sum(pro_samp.counts$Reads)
#how many 18S reads
nrow(pro_count_table)

#how many 18S reads
sum(euk_samp.counts$Reads)
#how many ASVs
nrow(euk_count_table)

#7. Convert data into phyloseq objects

pro_ASV = otu_table(as.matrix(pro_count_table), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(pro_metadata, row.names = rownames(pro_metadata)))
pro_TREE = pro_tree

#make into phyloseq object
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE)
pro_ps

euk_ASV = otu_table(as.matrix(euk_count_table), taxa_are_rows = TRUE)
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(euk_metadata, row.names = rownames(euk_metadata)))
euk_TREE = euk_tree

#make into phyloseq object
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE)
euk_ps

#8. Remove ambiguous annotations and singletons

pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukayote") & !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) >= 1, pro_ps.ambig) 

#compare before and after
pro_ps
pro_ps.ambig


#eukaryote dataset with micrometazoan's removed.
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Bacteria", "Archaea") & !is.na(Phylum) & !Phylum %in% c("Vertebrata", "Arthropoda") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
euk_ps.ambig <- prune_taxa(taxa_sums(euk_ps.ambig) >= 1, euk_ps.ambig) 

#compare before and after
euk_ps
euk_ps.ambig

#What was removed from our original phyloseq objects?

#kingdoms of the 16S dataset
table(pro_taxa_tab$Kingdom) #mostly bacteria

#orders of the 16S dataset
sort(table(pro_taxa_tab$Order), decreasing = TRUE) #quite a few chloroplasts

#kingdoms of the original 18S dataset
table(euk_taxa_tab$Kingdom) #mostly eukaryotes

#9. Re-root tree. Use this function > (https://john-quensen.com/r/unifrac-and-tree-roots/) < picks longest branch

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

# Run on my_phyloseq tree
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

#10. Normalise data

sort(sample_sums(pro_ps.ambig), decreasing = FALSE)
sort(sample_sums(euk_ps.ambig), decreasing = FALSE)

#as the minimum number of reads per environmental sample is 2300+ I'm not going to remove low abundance reads
#I'm going to use proportional transformation to normalise all data
pro.rel = transform_sample_counts(pro_ps.ambig, function(x) x/sum(x))
#sanity check
colSums(data.frame(otu_table(pro.rel)))
euk.rel = transform_sample_counts(euk_ps.ambig, function(x) x/sum(x))
#sanity check
colSums(data.frame(otu_table(euk.rel)))


#alternative method: Rarefy data
pro.rare = rarefy_even_depth(pro_ps.ambig, sample.size = 2402)
euk.rare = rarefy_even_depth(euk_ps.ambig, sample.size = 2266)

#11. Save phyloseq objects
saveRDS(pro_ps.ambig, "../results/16S-phylo-object.rds")
saveRDS(euk_ps.ambig, "../results/18S-phylo-object.rds")
saveRDS(pro.rare, "../results/16S-phylo-object-norm-rare.rds")
saveRDS(euk.rare, "../results/18S-phylo-object-norm-rare.rds")
saveRDS(pro.rel, "../results/16S-phylo-object-norm-rel.rds")
saveRDS(euk.rel, "../results/18S-phylo-object-norm-rel.rds")
