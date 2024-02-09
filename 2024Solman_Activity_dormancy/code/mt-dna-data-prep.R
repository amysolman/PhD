
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
#install.packages("kableExtra")
#library(kableExtra)
library(RColorBrewer) #for plotting colours
library(tidyr) #wide to long format

#I'm going to import my previous amplicon phyloseq objects so I can compare the taxonomy
# pro_test <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds")
# euk_test <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds")

#2. Import data

metadata <- read.csv(file="../data/mt-metadata.csv", sep=",") #Metadata

#18S rRNA Gene Amplicon Data
euk_count_tab <- read_tsv(file="../data/qiime2-mt-dna/18S/table.tsv", skip=1) #Count Table
euk_count_tab <- column_to_rownames(euk_count_tab, var = colnames(euk_count_tab)[1]) #first col in our count tab is ASV IDs
euk_tax <- read_tsv(file="../data/qiime2-mt-dna/18S/taxonomy.tsv") #Taxonomy table
euk_fasta <- read.fasta(file="../data/qiime2-mt-dna/18S/dna-sequences.fasta") #Sequences
euk_tree <- read_tree("../data/qiime2-mt-dna/18S/tree.nwk") #Phylogenetic tree
names(euk_count_tab) = metadata$SampleID

#16S rRNA Gene Amplicon Data
pro_count_tab <- read_tsv(file="../data/qiime2-mt-dna/16S/table.tsv", skip=1) #Count Table
pro_count_tab <- column_to_rownames(pro_count_tab, var = colnames(pro_count_tab)[1]) #first col in our count tab is ASV IDs
pro_tax <- read_tsv(file="../data/qiime2-mt-dna/16S/taxonomy.tsv") #Taxonomy table
pro_fasta <- read.fasta(file="../data/qiime2-mt-dna/16S/dna-sequences.fasta") #Sequences
pro_tree <- read_tree("../data/qiime2-mt-dna/16S/tree.nwk") #Phylogenetic tree
names(pro_count_tab) = metadata$SampleID

#3. Sort metadata

#18S rRNA Gene Data
euk_meta = metadata[! names(metadata) %in% c("SampleID16S")] #remove unneccesary data from the metadata file
rownames(euk_meta) = euk_meta$SampleID #make sure row has the right names

#16S rRNA Gene Data
pro_meta = metadata[! names(metadata) %in% c("SampleID18S")] #remove unneccesary data from the metadata file
rownames(pro_meta) = pro_meta$SampleID #make sure row has the right names

#4. Sort taxonomy table

#Now to modify the taxonomic table for my data
euk_taxa_tab <- euk_tax %>%
  #mutate(euk_tax=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=";$", replacement="")) %>%
  mutate(euk_tax=str_replace_all(string=euk_tax, pattern=" ", replacement="")) %>%
  separate(euk_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

pro_taxa_tab <- pro_tax %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#let's look at the taxonomy of my other phylo objects (from PMA amplicon seq)
#pro.tax.test = data.frame(tax_table(pro_test))

#5. Output number of sequences per sample in raw reads

pro_count_tab_df = as.data.frame(colSums(pro_count_tab)) #make into df
pro_count_tab_df$Sample = rownames(pro_count_tab_df) #make sample names a column
names(pro_count_tab_df) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

pdf("../results/mt-dna-prokaryote-count-histogram.pdf", width=20, height=10)
print(pro_p)
dev.off()

#5. Output number of sequences per sample in raw reads

euk_count_tab_df = as.data.frame(colSums(euk_count_tab)) #make into df
euk_count_tab_df$Sample = rownames(euk_count_tab_df) #make sample names a column
names(euk_count_tab_df) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p

pdf("../results/mt-dna-eukaryote-count-histogram.pdf", width=20, height=10)
print(euk_p)
dev.off()

#6. Raw data basic stats: total reads per dataset, number of unique ASVs

euk_tot = sum(euk_count_tab_df$Reads) #total reads in the dataset
euk_asv = nrow(euk_count_tab) #number of unique ASVs

pro_tot = sum(pro_count_tab_df$Reads) #total reads in the dataset
pro_asv = nrow(pro_count_tab) #number of unique ASVs


#9. Convert data into phyloseq objects

#18S
euk_ASV = otu_table(as.matrix(euk_count_tab), taxa_are_rows = TRUE) 
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(euk_meta, row.names = rownames(euk_meta)))
euk_TREE = euk_tree
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE) #make into phyloseq object
euk_ps

#16S
pro_ASV = otu_table(as.matrix(pro_count_tab), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(pro_meta, row.names = rownames(pro_meta)))
pro_TREE = pro_tree
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE) #make into phyloseq object
pro_ps

#10. Remove ambiguous annotations and singletons
#do we have chloroplast/mitochondria in our 16S data?
# sort(get_taxa_unique(pro_ps, "Order")) == "Chloroplast"

pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukaryota") & !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) >= 1, pro_ps.ambig)

#eukaryote dataset with micrometazoan's removed.
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !is.na(Phylum) & !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa", "Arthropoda") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
euk_ps.ambig <- prune_taxa(taxa_sums(euk_ps.ambig) >= 1, euk_ps.ambig) 

#total eukaryote dataset
euk_ps.tot <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !is.na(Phylum) & !Phylum %in% c("Vertebrata", "Arthropoda") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
euk_ps.tot <- prune_taxa(taxa_sums(euk_ps.tot) >= 1, euk_ps.tot) 

# euk.check = data.frame(otu_table(euk_ps.ambig))

# sort(sample_sums(pro_ps.ambig))
# sort(sample_sums(euk_ps.ambig))

#make a micrometazoa phyloseq object
ps.mm <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
ps.mm

#Micrometazoa community 
plot_bar(ps.mm, "Class", fill="Genus", facet_grid=~Habitat2)

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

# Run on my_phyloseq tree
#16S
pro_tree_new <- phy_tree(pro_ps.ambig)
pro_new_outgroup <- pick_new_outgroup(pro_tree_new)
pro_new_tree_root <- ape::root(pro_tree_new, outgroup=pro_new_outgroup, resolve.root=TRUE) # Re-root tree
pro_new_tree_dich <- ape::multi2di(pro_new_tree_root) # Convert to dichotomy tree
phy_tree(pro_ps.ambig) <- pro_new_tree_dich

#18S without micro
euk_tree_new <- phy_tree(euk_ps.ambig)
euk_new_outgroup <- pick_new_outgroup(euk_tree_new)
euk_new_tree_root <- ape::root(euk_tree_new, outgroup=euk_new_outgroup, resolve.root=TRUE) # Re-root tree
euk_new_tree_dich <- ape::multi2di(euk_new_tree_root) # Convert to dichotomy tree
phy_tree(euk_ps.ambig) <- euk_new_tree_dich

#18S with micro
euk_tree_new.tot <- phy_tree(euk_ps.tot)
euk_new_outgroup.tot <- pick_new_outgroup(euk_tree_new.tot)
euk_new_tree_root.tot <- ape::root(euk_tree_new.tot, outgroup=euk_new_outgroup.tot, resolve.root=TRUE) # Re-root tree
euk_new_tree_dich.tot <- ape::multi2di(euk_new_tree_root.tot) # Convert to dichotomy tree
phy_tree(euk_ps.tot) <- euk_new_tree_dich.tot

#12. Filtered data basic stats: total reads per dataset, number of unique ASVs

filt_euk_tot = sum(colSums(data.frame(otu_table(euk_ps.ambig))))
filt_euk_asv = nrow(data.frame(otu_table(euk_ps.ambig)))

filt_euk_tot.tot = sum(colSums(data.frame(otu_table(euk_ps.tot))))
filt_euk_asv.tot = nrow(data.frame(otu_table(euk_ps.tot)))

filt_pro_tot = sum(colSums(data.frame(otu_table(pro_ps.ambig))))
filt_pro_asv = nrow(data.frame(otu_table(pro_ps.ambig)))

#13. Save phyloseq objects
saveRDS(pro_ps.ambig, "../results/mt-dna-16S-phylo-object.rds")
saveRDS(euk_ps.ambig, "../results/mt-dna-18S-phylo-object-micro-remove.rds")
saveRDS(euk_ps.tot, "../results/mt-dna-18S-phylo-object-micro-keep.rds")
saveRDS(ps.mm, "../results/mt-dna-18S-phylo-object-micro-only.rds")