#1. Clear workspace and load packages

rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
library(decontam)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)
library(RColorBrewer) #for plotting colours
library(tidyr) #wide to long format
library(microbiome)

#2. Import data

metadata <- read.csv(file="../data/metadata.csv", sep=",") #Metadata

#18S rRNA Gene Amplicon Data
euk_count_tab <- read_tsv(file="../data/18S/table.tsv", skip=1) #Count Table
euk_count_tab <- column_to_rownames(euk_count_tab, var = colnames(euk_count_tab)[1]) #first col in our count tab is ASV IDs
euk_tax <- read_tsv(file="../data/18S/taxonomy.tsv") #Taxonomy table
euk_fasta <- read.fasta(file="../data/18S/dna-sequences.fasta") #Sequences
euk_tree <- read_tree("../data/18S/tree.nwk") #Phylogenetic tree
names(euk_count_tab) = metadata$SampleID #match sample names in count table and metadata

#16S rRNA Gene Amplicon Data
pro_count_tab <- read_tsv(file="../data/16S/table.tsv", skip=1) #Count Table
pro_count_tab <- column_to_rownames(pro_count_tab, var = colnames(pro_count_tab)[1]) #first col in our count tab is ASV IDs
pro_tax <- read_tsv(file="../data/16S/taxonomy.tsv") #Taxonomy table
pro_fasta <- read.fasta(file="../data/16S/dna-sequences.fasta") #Sequences
pro_tree <- read_tree("../data/16S/tree.nwk") #Phylogenetic tree
names(pro_count_tab) = metadata$SampleID #match sample names in count table and metadata


#3.Remove sample S21.76 - this sample we contaminated during library preparation

#from metadata
metadata = metadata[!metadata$SampleID == "S21.76",]

#16S 
#asv table
pro_count_tab = pro_count_tab[,!names(pro_count_tab) == "S21.76"]
#remove any ASVs that now have zero counts
pro_count_tab = pro_count_tab[rowSums(pro_count_tab) > 0,]
#tab table
pro_tax = pro_tax[pro_tax$`Feature ID` %in% rownames(pro_count_tab),]
#tree
pro_tree <- prune_taxa(pro_count_tab$`#OTU ID`, pro_tree)
#fasta
pro_fasta = pro_fasta[rownames(pro_count_tab)]


#18S
#asv table
euk_count_tab = euk_count_tab[,!names(euk_count_tab) == "S21.76"]
#remove any ASVs that now have zero counts
euk_count_tab = euk_count_tab[rowSums(euk_count_tab) > 0,]
#tax table
euk_tax = euk_tax[euk_tax$`Feature ID` %in% rownames(euk_count_tab),]
#tree
euk_tree <- prune_taxa(euk_count_tab$`#OTU ID`, euk_tree)
#fasta
euk_fasta = euk_fasta[rownames(euk_count_tab)]

#4. Sort metadata

#18S rRNA Gene Data
euk_meta = metadata[! names(metadata) %in% c("SampleID_16S")] #remove unneccesary data from the metadata file
rownames(euk_meta) = euk_meta$SampleID #make sure row has the right names

#16S rRNA Gene Data
pro_meta = metadata[! names(metadata) %in% c("SampleID_18S")] #remove unneccesary data from the metadata file
rownames(pro_meta) = pro_meta$SampleID #make sure row has the right names


#5. Sort taxonomy table - we want our Qiime taxonomy output to be in a form appropriate for building a phyloseq object

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


#6. Output number of sequences per sample in raw reads

#16S
pro_count_tab_df = as.data.frame(colSums(pro_count_tab)) #make into df
pro_count_tab_df$Sample = rownames(pro_count_tab_df) #make sample names a column
names(pro_count_tab_df) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

#18S
euk_count_tab_df = as.data.frame(colSums(euk_count_tab)) #make into df
euk_count_tab_df$Sample = rownames(euk_count_tab_df) #make sample names a column
names(euk_count_tab_df) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p


#7. Raw data basic stats: total reads per dataset, number of unique ASVs

#18S
euk_tot = sum(euk_count_tab_df$Reads) #total reads in the datazset
euk_asv = nrow(euk_count_tab) #number of unique ASVs

pro_tot = sum(pro_count_tab_df$Reads) #total reads in the datazset
pro_asv = nrow(pro_count_tab) #number of unique ASVs


#8. Convert data into phyloseq objects

#16S
pro_ASV = otu_table(as.matrix(pro_count_tab), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(pro_meta, row.names = rownames(pro_meta)))
pro_TREE = pro_tree
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE) #make into phyloseq object

#18S
euk_ASV = otu_table(as.matrix(euk_count_tab), taxa_are_rows = TRUE) 
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(euk_meta, row.names = rownames(euk_meta)))
euk_TREE = euk_tree
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE) #make into phyloseq object


#9. Remove ambiguous annotations and singletons

#16S
pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukaryota") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) >= 10, pro_ps.ambig)

#eukaryote (18S) dataset with micrometazoan's removed.
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta"))
euk_ps.ambig <- prune_taxa(taxa_sums(euk_ps.ambig) >= 10, euk_ps.ambig) 

#total eukaryote (18S) dataset
euk_ps.tot <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata") & !Class %in% c("Embryophyta"))
euk_ps.tot <- prune_taxa(taxa_sums(euk_ps.tot) >= 10, euk_ps.tot) 

#make a micrometazoa phyloseq object
ps.mm <- subset_taxa(euk_ps, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta"))
ps.mm <- prune_taxa(taxa_sums(ps.mm) >= 10, ps.mm) 


#10. Re-root tree. Use this function > (https://john-quensen.com/r/unifrac-and-tree-roots/) < picks longest branch

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

#18S micrometazoans only
ps.mm_tree_new <- phy_tree(ps.mm)
ps.mm_outgroup <- pick_new_outgroup(ps.mm_tree_new)
ps.mm_root <- ape::root(ps.mm_tree_new, outgroup=ps.mm_outgroup, resolve.root=TRUE) # Re-root tree
ps.mm_dich <- ape::multi2di(ps.mm_root) # Convert to dichotomy tree
phy_tree(ps.mm) <- ps.mm_dich


#11. Filtered data basic stats: total reads per dataset, number of unique ASVs

#16S
filt_pro_tot = sum(colSums(data.frame(otu_table(pro_ps.ambig))))
filt_pro_asv = nrow(data.frame(otu_table(pro_ps.ambig)))

#18S without micrometazoans
filt_euk_tot = sum(colSums(data.frame(otu_table(euk_ps.ambig))))
filt_euk_asv = nrow(data.frame(otu_table(euk_ps.ambig)))

#all 18S
filt_euk_tot.tot = sum(colSums(data.frame(otu_table(euk_ps.tot))))
filt_euk_asv.tot = nrow(data.frame(otu_table(euk_ps.tot)))

#micrometazoans
filt_ps.mm_tot = sum(colSums(data.frame(otu_table(ps.mm))))
filt_ps.mm_asv = nrow(data.frame(otu_table(ps.mm)))


#12. Save phyloseq objects
saveRDS(pro_ps.ambig, "../results/16S-phylo-object.rds")
saveRDS(euk_ps.ambig, "../results/18S-phylo-object-micro-remove.rds")
saveRDS(euk_ps.tot, "../results/18S-phylo-object-micro-keep.rds")
saveRDS(ps.mm, "../results/18S-phylo-object-micro-only.rds")
