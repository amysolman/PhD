#1. Clear workspace and load packages


rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
#BiocManager::install("decontam")
library(decontam)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)


#2. Import data


#16S

#Count Table
pro_count_table <- read_tsv(file="../data/16S/export/table.tsv", skip=1)

#Specify the first column in our count table is the row names (ASV IDs)
pro_count_table <- column_to_rownames(pro_count_table, var = colnames(pro_count_table)[1])
pro_sample_names = as.numeric(names(pro_count_table))

#Rest of the data
pro_tax <- read_tsv(file="../data/16S/export/taxonomy.tsv")
pro_fasta <- read.fasta(file="../data/16S/export/dna-sequences.fasta")
pro_tree <- read_tree("../data/16S/export/tree.nwk")

#18S

#Count Table
euk_count_table <- read_tsv(file="../data/18S/export-11/table.tsv", skip=1)

#Specify the first column in our count table is the row names (ASV IDs)
euk_count_table <- column_to_rownames(euk_count_table, var = colnames(euk_count_table)[1])
#sample_names = as.numeric(names(count_table))

#Rest of the data
euk_tax <- read_tsv(file="../data/18S/export-11/taxonomy.tsv")
euk_fasta <- read.fasta(file="../data/18S/export-11/dna-sequences.fasta")
euk_tree <- read_tree("../data/18S/export-11/tree.nwk")

#Metadata
metadata <- read.csv(file="../data/metadata.csv", sep=",")

#3. Wrangle metadata
#Sort metadata and count dataframes so they have the same sample names in the same order

#16S
#remove unnecessary data from the metadata file
pro_metadata = metadata[! names(metadata) %in% c("SampleID_18S")]
#rename columns
names(pro_metadata)[names(pro_metadata) == 'SampleID_16S'] <- 'SampleID'
#make sure row has the right names
rownames(pro_metadata) = pro_metadata$SampleID

#18S
#remove unnecessary data from the metadata file
euk_metadata = metadata[! names(metadata) %in% c("SampleID_16S")]
#rename columns
names(euk_metadata)[names(euk_metadata) == 'SampleID_18S'] <- 'SampleID'
#make sure row has the right names
rownames(euk_metadata) = euk_metadata$SampleID

#4. Sort taxonomy table so our output from QIIME classification step is useful for analysis

#16S
pro_taxa_tab <- pro_tax %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#18S
euk_taxa_tab <- euk_tax %>%
  mutate(euk_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=";$", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=" ", replacement="")) %>%
  separate(euk_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
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


#7. Identify contaminants + print as table

#16S

# convert count table to dataframe
count_table_df <- as.data.frame(pro_count_table)
#what is the sample ID of the Negative Control sample
pro_neg.samp = pro_metadata[pro_metadata$Name == "NegCtrlBlank16",]
# A logical vector showing which entries are 'TRUE' for controls and 'FALSE' for genuine samples
pro_decontam <- startsWith(names(count_table_df), rownames(pro_neg.samp))
# Run isContaminant, while transforming the matrix with t()
contam_df <- isContaminant(t(count_table_df), neg=pro_decontam)
# What do our results look like?
table(contam_df$contaminant) 
# Make a vector holding contaminant IDs
pro_contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
# Check the taxonomy of these contaminants
pro_contam_tab <- pro_taxa_tab[rownames(pro_taxa_tab) %in% pro_contam_asvs, ]
write.csv(pro_contam_tab, "../results/16S-contaminant-asvs.csv")


#18S

# convert count table to dataframe
count_table_df <- as.data.frame(euk_count_table)
#what is the sample ID of the Negative Control sample
euk_neg.samp = euk_metadata[euk_metadata$Name == "NegCtrlBlank16",]
# A logical vector showing which entries are 'TRUE' for controls and 'FALSE' for genuine samples
euk_decontam <- startsWith(names(count_table_df), rownames(euk_neg.samp))
# Run isContaminant, while transforming the matrix with t()
contam_df <- isContaminant(t(count_table_df), neg=euk_decontam)
# What do our results look like?
table(contam_df$contaminant) 
# Make a vector holding contaminant IDs
euk_contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
# Check the taxonomy of these contaminants
euk_contam_tab <- euk_taxa_tab[rownames(euk_taxa_tab) %in% euk_contam_asvs, ]
write.csv(euk_contam_tab, "../results/18S-contaminant-asvs.csv")

#8. Remove contaminants from data objects (I didn't remove anything here)

#16S
#COUNT TABLE
# x <- rownames(pro_count_table) %in% pro_contam_asvs #which rows are contaminants?
# pro_count_table <- pro_count_table[!x,] #remove those rows
# 
# #what is the index of "TRUE" values in decontam?
# true_index <- NULL
# for (i in 1:length(pro_decontam)){
#   if (pro_decontam[[i]] == "TRUE"){
#     true_index <- c(true_index, i)
#   }
# }
# pro_count_table <- pro_count_table[,-c(true_index)]  #remove those samples

#METADATA
row_names_df_to_remove<-c(rownames(pro_neg.samp))
pro_metadata = pro_metadata[!(row.names(pro_metadata) %in% row_names_df_to_remove),]

###TAXONOMY FILE
#pro_taxonomy <- pro_taxa_tab[ ! rownames(pro_taxa_tab) %in% pro_contam_asvs, ]
# 
###FASTA FILE
#Remove contaminants from dna sequences file
#pro_fasta = pro_fasta[pro_contam_asvs] <- NULL

###PHYLOGENETIC TREE
# This is a phyloseq object so we can use a phyloseq command called prune_taxa
#pro_tree <- prune_taxa(pro_count_table$`#OTU ID`, pro_tree)


#18S
#COUNT TABLE
# x <- rownames(euk_count_table) %in% euk_contam_asvs #which rows are contaminants?
# euk_count_table <- euk_count_table[!x,] #remove those rows
# 
# #what is the index of "TRUE" values in decontam?
# true_index <- NULL
# for (i in 1:length(euk_decontam)){
#   if (euk_decontam[[i]] == "TRUE"){
#     true_index <- c(true_index, i)
#   }
# }
# euk_count_table <- euk_count_table[,-c(true_index)]  #remove those samples

#METADATA
row_names_df_to_remove<-c(rownames(euk_neg.samp))
euk_metadata = euk_metadata[!(row.names(euk_metadata) %in% row_names_df_to_remove),]

###TAXONOMY FILE
#euk_taxonomy <- euk_taxa_tab[ ! rownames(euk_taxa_tab) %in% euk_contam_asvs, ]
# 
###FASTA FILE
#Remove contaminants from dna sequences file
#euk_fasta = euk_fasta[euk_contam_asvs] <- NULL

###PHYLOGENETIC TREE
# This is a phyloseq object so we can use a phyloseq command called prune_taxa
#euk_tree <- prune_taxa(euk_count_table$`#OTU ID`, euk_tree)


#9. Convert data into phyloseq objects

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


#10. Remove ambiguous annotations and singletons

#16S
pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukayote") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) > 1, pro_ps.ambig) 

#18S dataset with micrometazoan's removed
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
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


#13. Save phyloseq objects

#16S
saveRDS(pro_ps.ambig, "../results/16S-phylo-object.rds")

#18S
saveRDS(euk_ps.ambig, "../results/18S-phylo-object.rds")
