
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

#2. Import data

metadata <- read.csv(file="../data/metadata.csv", sep=",") #Metadata

#18S rRNA Gene Amplicon Data
euk_count_tab <- read_tsv(file="../data/18S/table.tsv", skip=1) #Count Table
euk_count_tab <- column_to_rownames(euk_count_tab, var = colnames(euk_count_tab)[1]) #first col in our count tab is ASV IDs
euk_tax <- read_tsv(file="../data/18S/taxonomy.tsv") #Taxonomy table
euk_fasta <- read.fasta(file="../data/18S/dna-sequences.fasta") #Sequences
euk_tree <- read_tree("../data/18S/tree.nwk") #Phylogenetic tree
names(euk_count_tab) = metadata$SampleID

#16S rRNA Gene Amplicon Data
pro_count_tab <- read_tsv(file="../data/16S/table.tsv", skip=1) #Count Table
pro_count_tab <- column_to_rownames(pro_count_tab, var = colnames(pro_count_tab)[1]) #first col in our count tab is ASV IDs
pro_tax <- read_tsv(file="../data/16S/taxonomy.tsv") #Taxonomy table
pro_fasta <- read.fasta(file="../data/16S/dna-sequences.fasta") #Sequences
pro_tree <- read_tree("../data/16S/tree.nwk") #Phylogenetic tree
names(pro_count_tab) = metadata$SampleID


#Check metadata for outliers
#Outliers have been removed as described in my metadata document

#chemical data to check
chem.keep = c("pH", "Conductivity_muS", "Temp")
meta.check = metadata[,names(metadata) %in% chem.keep]

for (i in 1:ncol(meta.check)){
  
  var = meta.check[,i]
  name = names(meta.check)[i]
  boxplot(var, main=name)
  hist(var, main=name, sub=paste0("Min: ", min(var, na.rm = TRUE), ", Max: ", max(var, na.rm=TRUE)), xlab=" ")
}


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
  separate(euk_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

pro_taxa_tab <- pro_tax %>%
  mutate(pro_tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=";$", replacement="")) %>%
  mutate(pro_tax=str_replace_all(string=pro_tax, pattern=" ", replacement="")) %>%
  separate(pro_tax, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon, -Confidence) %>%
  column_to_rownames(var='Feature ID')

#5. Output number of sequences per sample in raw reads

pro_count_tab_df = as.data.frame(colSums(pro_count_tab)) #make into df
pro_count_tab_df$Sample = rownames(pro_count_tab_df) #make sample names a column
names(pro_count_tab_df) = c("Reads", "Sample")

pro_p<-ggplot(data=pro_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pro_p

#5. Output number of sequences per sample in raw reads

euk_count_tab_df = as.data.frame(colSums(euk_count_tab)) #make into df
euk_count_tab_df$Sample = rownames(euk_count_tab_df) #make sample names a column
names(euk_count_tab_df) = c("Reads", "Sample")

euk_p<-ggplot(data=euk_count_tab_df, aes(x=Sample, y=Reads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

euk_p


#6. Raw data basic stats: total reads per dataset, number of unique ASVs

euk_tot = sum(euk_count_tab_df$Reads) #total reads in the datazset
euk_asv = nrow(euk_count_tab) #number of unique ASVs

pro_tot = sum(pro_count_tab_df$Reads) #total reads in the datazset
pro_asv = nrow(pro_count_tab) #number of unique ASVs

# pro_tot = sum(pro_samp.counts$Reads)
# pro_asv = nrow(pro_count_table)


# df = data.frame("Dataset"=c("16S", "18S"), "Total Reads" = c(pro_tot, euk_tot), "Unique ASVs" = c(pro_asv, euk_asv))
# names(df) = c("Dataset", "Total Reads", "Unique ASVs")


#9. Convert data into phyloseq objects

#18S
euk_ASV = otu_table(as.matrix(euk_count_tab), taxa_are_rows = TRUE) 
euk_TAX = tax_table(as.matrix(euk_taxa_tab))
euk_META = sample_data(data.frame(euk_meta, row.names = rownames(euk_meta)))
euk_TREE = euk_tree
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE) #make into phyloseq object

#16S
pro_ASV = otu_table(as.matrix(pro_count_tab), taxa_are_rows = TRUE)
pro_TAX = tax_table(as.matrix(pro_taxa_tab))
pro_META = sample_data(data.frame(pro_meta, row.names = rownames(pro_meta)))
pro_TREE = pro_tree
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE) #make into phyloseq object

#10. Remove ambiguous annotations and singletons
#do we have chloroplast/mitochondria in our 16S data?
# sort(get_taxa_unique(pro_ps, "Order")) == "Chloroplast"

pro_ps.ambig <- subset_taxa(pro_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro_ps.ambig <- prune_taxa(taxa_sums(pro_ps.ambig) >= 10, pro_ps.ambig)

#eukaryote dataset with micrometazoan's removed.
euk_ps.ambig <- subset_taxa(euk_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta"))
euk_ps.ambig <- prune_taxa(taxa_sums(euk_ps.ambig) >= 10, euk_ps.ambig) 

#total eukaryote dataset
euk_ps.tot <- subset_taxa(euk_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Bacteria", "Archaea") & !Phylum %in% c("Vertebrata") & !Class %in% c("Embryophyta"))
euk_ps.tot <- prune_taxa(taxa_sums(euk_ps.tot) >= 10, euk_ps.tot) 

# euk.check = data.frame(otu_table(euk_ps.ambig))

# sort(sample_sums(pro_ps.ambig))
# sort(sample_sums(euk_ps.ambig))

#make a micrometazoa phyloseq object
ps.mm <- subset_taxa(euk_ps, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Bacteria", "Archaea") & Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa") & !Class %in% c("Embryophyta"))
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


#What is in our blanks?
cont = c("Blank1", "Blank2", "Blank3", "Blank4", "Blank5", "Blank6", "Blank7", "Blank8")

#extract the count table with these samples
#get count and tax tables from our phyloseq obejcts
pro.count.filt = data.frame(otu_table(pro_ps.ambig))
euk.count.filt = data.frame(otu_table(euk_ps.ambig))
pro.tax.filt = data.frame(tax_table(pro_ps.ambig))
euk.tax.filt = data.frame(tax_table(euk_ps.ambig))

pro.count.blanks.only <- pro.count.filt[,names(pro.count.filt) %in% cont]
# names(pro.count.blanks.only) = c("LabBlank1", "LabBlank2", "LabBlank3", "FieldBlank1", "FieldBlank2", "FieldBlank3","FieldBlank4", "FieldBlank5")
euk.count.blanks.only <- euk.count.filt[,names(euk.count.filt) %in% cont]
# names(euk.count.blanks.only) = c("LabBlank1", "LabBlank2", "LabBlank3", "FieldBlank1", "FieldBlank2", "FieldBlank3","FieldBlank4", "FieldBlank5")
tab.sum.blank.reads = as.data.frame(rbind(colSums(pro.count.blanks.only), colSums(euk.count.blanks.only)))
tab.sum.blank.reads = data.frame(cbind(c("16S", "18S"), tab.sum.blank.reads))
names(tab.sum.blank.reads) = c("Amplicon", "Blank1", "Blank2", "Blank3", "Blank4", "Blank5", "Blank6", "Blank7", "Blank8")

#how many ASVs per blank?
pro.cont = pro.count.blanks.only
pro.cont[pro.cont >= 1] = 1
max(pro.cont)
colSums(pro.cont)
#how many ASVs in total 
x = pro.cont[rowSums(pro.cont) > 0,]
#what is the taxonomy of organisms in our control samples?
pro.cont.tax = pro.tax.filt[rownames(pro.tax.filt) %in% rownames(x),]
#merge at same taxonomic level
pro.cont.tax.m = pro.cont.tax %>%
  distinct()

euk.cont = euk.count.blanks.only
euk.cont[euk.cont >= 1] = 1
max(euk.cont)
colSums(euk.cont)
#how many ASVs in total
y = euk.cont[rowSums(euk.cont) > 0,]
#what is the taxonomy of organisms in our control samples?
euk.cont.tax = euk.tax.filt[rownames(euk.tax.filt) %in% rownames(y),]
euk.cont.tax.m = euk.cont.tax %>%
  distinct()

#combine dataframes to export for checking taxa
cont.tax.out = rbind(pro.cont.tax.m, euk.cont.tax.m)

write.csv(cont.tax.out, "../results/control-sample-ASV-tax.csv")

#combine with tax info
pro.blank.count.tax = cbind(pro.tax.filt, pro.count.blanks.only)
#to long format
pro_data_long <- gather(pro.blank.count.tax, Sample, Count, Blank1:Blank8, factor_key=TRUE)
pro_data_long$Amplicon = "16S"

#combine with tax info
euk.blank.count.tax = cbind(euk.tax.filt, euk.count.blanks.only)
#to long format
euk_data_long <- gather(euk.blank.count.tax, Sample, Count, Blank1:Blank8, factor_key=TRUE)
euk_data_long$Amplicon = "18S"

#bind for plotting
# data.2.plot = rbind(pro_data_long, euk_data_long)

#remove rows with zeros
pro_data_long.no.0 = pro_data_long[pro_data_long$Count != 0, ]
euk_data_long.no.0 = euk_data_long[euk_data_long$Count != 0, ]

pro.getPalette = colorRampPalette(brewer.pal(length(unique(pro_data_long.no.0$Genus)), "Paired"))
euk.getPalette = colorRampPalette(brewer.pal(length(unique(euk_data_long.no.0$Genus)), "Paired"))

plot.cols = c("cyan",
              "darkred",
              "green",
              "red",
              "darkgreen",
              "gold", 
              "deeppink",
              "gold4", 
              "lightblue", 
              "blue4", 
              "grey", 
              "firebrick",
              "grey30",
              "blueviolet",
              "yellow1",
              "chartreuse4",
              "cornsilk",
              "pink", 
              "aquamarine2",
              "deeppink4",
              "lemonchiffon3",
              "cornflowerblue",
              "orange", 
              "blue",
              "darkorange4", 
              "paleturquoise", 
              "turquoise4", 
              "tan",
              "darkorchid4",
              "tan4",
              "palegreen",
              "sienna3",
              "cyan4",
              "antiquewhite2",
              "hotpink",
              "darkolivegreen",
              "royalblue4",
              "tomato",
              "palegreen4",
              "goldenrod",
              "lightcyan",
              "mediumslateblue",
              "mediumseagreen",
              "red4",
              "steelblue3",
              "pink3",
              "springgreen",
              "darkorange",
              "ghostwhite",
              "grey53",
              "plum")


#16S Contaminants


#plot
ggplot(pro_data_long.no.0, aes(fill=Genus, y=Count, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  #facet_grid(~Amplicon)+
  scale_fill_manual(values = plot.cols)+
  guides(fill=guide_legend(ncol=3,byrow=FALSE))+
  theme(legend.posi="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

ggplot(euk_data_long.no.0, aes(fill=Genus, y=Count, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  #facet_grid(~Amplicon)+
  scale_fill_manual(values = plot.cols)+
  guides(fill=guide_legend(ncol=3,byrow=FALSE))+
  theme(legend.posi="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())


#7. Identify contaminants + print as table
cont = c("Blank1", "Blank2", "Blank3", "Blank4", "Blank5", "Blank6", "Blank7", "Blank8") #control sample IDs
# cont = c("Blank1", "Blank2", "Blank3") #control sample IDs
euk_decontam = names(euk.count.filt) %in% cont #get local vector of control samples
contam_df <- isContaminant(t(euk_count_tab), neg=euk_decontam) # Run isContaminant, while transforming the matrix with t()
table(contam_df$contaminant) # What do our results look like?
euk_contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ]) # Make a vector holding contaminant IDs
euk_contam_tab <- euk.tax.filt[rownames(euk.tax.filt) %in% euk_contam_asvs, ] # Check the taxonomy of these contaminants
write.csv(euk_contam_tab, "../results/18S-contaminant-asvs.csv")

pro_decontam = names(pro.count.filt) %in% cont #get local vector of control samples
pro_contam_df <- isContaminant(t(pro.count.filt), neg=pro_decontam) # Run isContaminant, while transforming the matrix with t()
table(pro_contam_df$contaminant) # What do our results look like?
pro_contam_asvs <- row.names(pro_contam_df[pro_contam_df$contaminant == TRUE, ]) # Make a vector holding contaminant IDs
pro_contam_tab <- pro.tax.filt[rownames(pro.tax.filt) %in% pro_contam_asvs, ] # Check the taxonomy of these contaminants
write.csv(pro_contam_tab, "../results/16S-contaminant-asvs.csv")

#13. Save phyloseq objects
saveRDS(pro_ps.ambig, "../results/16S-phylo-object.rds")
saveRDS(euk_ps.ambig, "../results/18S-phylo-object-micro-remove.rds")
saveRDS(euk_ps.tot, "../results/18S-phylo-object-micro-keep.rds")
saveRDS(ps.mm, "../results/18S-phylo-object-micro-only.rds")