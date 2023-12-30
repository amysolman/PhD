#clear workspace and load package
rm(list=ls())

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
library(microbiome)
library(cowplot)
#library(fantaxtic)

#load data
metadata <- read.csv(file="../data/metadata.csv", sep=",") #Metadata
ps.pro = readRDS("../results/16S-phylo-object.rds")
ps.euk.mm.rm = readRDS("../results/18S-phylo-object-micro-remove.rds")
ps.euk.mm = readRDS("../results/18S-phylo-object-micro-keep.rds")
ps.mm = readRDS("../results/18S-phylo-object-micro-only.rds")

#Dealing with Contaminants

# We can start with the preposition that contamination in our samples is minimal 
#and those ASVs within our negative controls are mainly carry over from our TRUE samples 
#with a small degree of contamination from sample processing reagents/equipment. 
# We can explore this by looking at the presence of bands in gels of our PCR products. 
#After amplifying my samples/blanks I had no visible bands in the blanks. Tapestation of 
#blanks after indexing PCR showed some contamination - this is why these samples were sequenced. 


#Methods

# 1) Do nothing. If your gels are clear of bands/Qubit doesn't register DNA concentration then the impact of contamination is likely negligible (https://www.tandfonline.com/doi/full/10.1657/AAAR0015-062).
# 2) Remove all ASVs that are found in blanks from true samples.
# 3) Look at taxonomy of ASVs in blanks and remove human-associated taxa only.
# 4) Remove samples with amplification levels below controls (https://ami-journals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.14366?casa_token=NS2g2Mg8PaUAAAAA%3A4gDRC6RPOXMyB0Wj16e2ynfggh0Coz1-oKN0KgFuJl-FFDKUtkKRe3b8UIS2ugZeCFAdxDL_ZTisgMg)
# 5) USE PREVELANCE METHODS
# 5.1) Use decontam to identify contaminants and remove them.
# I have done this using the prevalence based method. In this method the presence/absence of each ASV in the TRUE samples is compared to the presence/abundance in negative controls to identify contaminants. Those ASVs with greater prevalence in control samples than TRUE samples are considered contaminants with a probability threshold of p < 0.1. I think this means non-contaminants are those with with prevalence in TRUE samples 10x higher than in negative control samples. For example, prevalence of 1 in negative control samples and 20 in TRUE samples would be a non-contaminant. Prevalence of 3 in negative control samples and 20 in TRUE samples would be a contaminant.
# 5.2) For each ASV in the negative controls, calculate mean sum of reads in negative controls and divide by the mean sum of reads across control and TRUE samples. Exclude ASVs with >5% abundance in negative controls. (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.13344) For example:
# 50/500 = 0.1 (=10%) so this ASV would be excluded.
# 10/500 = 0.02 (=2%) so this ASV would be retained.
# 5.3) Find ASVs in blanks that represent more than 0.1% of reads in all blanks and remove from true samples. Remove ASVs from TRUE samples with less reads than the total number of reads in that sample that come from contaminant ASVs (https://journals.asm.org/doi/full/10.1128/AEM.01253-17).
# 5.4) Find ASVs in blanks that represent more than 0.05% of all reads in blanks and remove from true samples (https://tc.copernicus.org/articles/12/3653/2018/tc-12-3653-2018.pdf).
# 5.5) Remove ASVs that make up >=1% of sequences in negative controls AND >=1% of sequences in TRUE samples (https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/jzo.12832?casa_token=XFXCN_1ZenIAAAAA%3ApmvIns4b3wVSw08f961saH_sZ7nZt4Ggkj1ZA79cQZBejcTp-6zRo_yMtZIhdSfcslUhTM512HjYRu8).
# 
# My Method
# 
# 1. Compare library sizes of TRUE and control samples. Control samples should have low read depth (see controls-lib-size.R)
# 2. Make pofile bar plots to get a general idea of the taxa in our controls (see controls-profiles.R)
# 3. Find ASVs that are in blanks. Of those that show high read abundance in TRUE samples (>0.1% total reads) 
# check that they are expected in our environment and not human associated. 
# 4. If they are expected then we will assume these are carry over from TRUE samples. 
# If they are human associated they will be manually removed from the dataset (see controls-remove.R).
# 5. For each remaining ASV in the negative controls, calculate mean sum of reads in negative controls and divide by the mean sum of reads across control and TRUE samples. 
# Exclude ASVs with >1% abundance in negative controls (see controls-remove.R).
# 6. Compare datasets before and after removing ASVs in steps 4 and 5 using NMDS plots and ANOSIM (see controls-NMDS.R).

#Remove contaminants from individual samples.

#subset phyloseq objects 
#ps <- ps.pro

remove.tax <- function(ps){
  
  genus.rm = c("Paracoccus", "Pantoea", "Brevundimonas", "Corynebacterium", "Acinetobacter", "Neisseria", "Rothia", "Streptococcus", "Chryseobacterium", "Actinomyces", "Enterococcus", "Itersonilia", "Neoascochyta", "Pleospora", "Malassezia", "Cladosporium", "Penicillium", "Saccharomyces", "Pleurotus", "Haplotaxida")
  fam.rm = c("Peronosporomycetes", "Aspergillaceae")
  
  otu <- otu_table(ps)
  sam <- sample_data(ps)
  tax <- data.frame(tax_table(ps))
  
  #identify ASVs that correspond to the groups we want to remove
  bad_taxa <- c(rownames(tax[tax$Genus %in% genus.rm,]), rownames(tax[tax$Family %in% fam.rm,]))
  
  #remove "contaminant" ASVs from the phyloseq object
  my_subset <- subset(otu_table(ps), ! rownames(otu_table(ps)) %in% bad_taxa)
  new.ps <- merge_phyloseq(my_subset, tax_table(ps), sample_data(ps), phy_tree(ps))
  
  #finally remove blanks
  #new.ps = subset_samples(new.ps, Habitat != "Control")
  
  #remove any samples or ASVs with zero counts
  new.ps.t = prune_taxa(taxa_sums(new.ps) > 0, new.ps) 
  new.ps.s = prune_samples(sample_sums(new.ps.t)>0, new.ps.t)
  
  return(new.ps.s)
  
}

#run the function
ps.pro2 <- remove.tax(ps.pro)
ps.euk.mm2 <- remove.tax(ps.euk.mm)
ps.euk.mm.rm2 <- remove.tax(ps.euk.mm.rm)
ps.mm2 = remove.tax(ps.mm)

# For each remaining ASV in the negative controls, 
#calculate mean sum of reads in negative controls and divide by the mean sum of reads across control and TRUE samples. 
#Exclude ASVs with >1% abundance in negative controls.

# ps = ps.pro
# cont = c("Blank.Snow1", "Blank.Snow2")
# sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID

#identify contaminents to be removed
remove_contaminants <- function(ps, cont, sam){
  
  #subset out phyloseq objects by habitat and control samples
  ps.hab <- prune_samples(sam, ps)
  ps.cont = subset_samples(ps, SampleID %in% cont)
  
  #remove zero abundant ASVs
  ps.hab.rm = prune_taxa(taxa_sums(ps.hab) > 0, ps.hab)
  ps.cont.rm = prune_taxa(taxa_sums(ps.cont) > 0, ps.cont)
  
  #how many ASVs in our negative control samples?
  x = data.frame(otu_table(ps.cont.rm)) #get negative control count table
  #for each ASV in the negative controls find the mean number of reads
  neg.means = data.frame(rowMeans(x))
  neg.means$id = rownames(neg.means)
  names(neg.means) = c("NegMean", "id")
  
  #find the mean sum of reads of those ASVs in all samples
  y = data.frame(otu_table(ps.hab.rm))
  y.sub = y[rownames(y) %in% rownames(x),]
  #combine dataframes
  x.y = cbind(y.sub, x[rownames(x) %in% rownames(y.sub),])
  full.means = data.frame(rowMeans(x.y))
  full.means$id = rownames(full.means)
  names(full.means) = c("FullMean", "id")
  
  #join the results
  full.neg.means = inner_join(full.means, neg.means, by="id")
  
  #get our statistic
  full.neg.means$Stat = full.neg.means$NegMean/full.neg.means$FullMean
  
  #ASVs to remove
  asv.rm = full.neg.means[full.neg.means$Stat > 0.01,]$id
  
  return(asv.rm)
  
}

#ps = ps.pro2
#remove those contaminants
remove_all_contaminants <- function(ps){
  
  #get out phyloseq object data
  otu <- otu_table(ps)
  sam <- sample_data(ps)
  
  #Snow samples
  cont = c("Blank.Snow1", "Blank.Snow2")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat == "Snow")))$SampleID)
  
  #remove from snow samples
  samples1 <- sample_names(ps)[sam$Habitat == "Snow"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples1] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples1]
  }
  
  cont=c("Blank.Ice1")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat == "Spring Ice")))$SampleID)
  
  #remove from spring ice samples
  samples2 <- sample_names(ps)[sam$Habitat == "Spring Ice"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples2] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples2]
  }
  
  cont=c("Blank.E", "Blank.Ice2")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat == "Summer Ice")))$SampleID)
  
  #remove from summer ice samples
  samples3 <- sample_names(ps)[sam$Habitat == "Summer Ice"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples3] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples3]
  }
  
  cont=c("Blank.NC1", "Blank.NC2", "Blank.NC3", "Blank.NC4", "Blank.NC5", "Blank.Cryo")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat == "Cryoconite")))$SampleID)
  
  #remove from cryoconite samples
  samples4 <- sample_names(ps)[sam$Habitat == "Cryoconite"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples4] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples4]
  }
  
  cont=c("Blank.1", "Blank.2", "Blank.3", "Blank.A", "Blank.B", "Blank.C", "Blank.D")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat %in% c("Snow", "Spring Ice"))))$SampleID)
  
  #remove from snow and spring ice samples
  samples5 <- sample_names(ps)[sam$Habitat == "Snow" | sam$Habitat == "Spring Ice"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples5] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples5]
  }
  
  cont=c("Blank.4", "Blank.5", "Blank.6")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat %in% c("Cryoconite", "Summer Ice"))))$SampleID)
  
  #remove from cryoconite and summer ice samples
  samples6 <- sample_names(ps)[sam$Habitat == "Cryoconite" | sam$Habitat == "Summer Ice"]
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples6] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples6]
  }
  
  cont=c("Blank.PCR.NC")
  bad_taxa <- remove_contaminants(ps, cont=cont, sam=data.frame(sample_data(subset_samples(ps, Habitat %in% c("Snow", "Spring Ice", "Cryoconite", "Summer Ice"))))$SampleID)
  
  #remove from all samples
  samples7 <- sample_names(ps)
  # Now, set the bad_taxa in the snow_samples to 0. Note, need to check if taxa are rows!
  taxa_are_rows(ps)
  #> [1] TRUE
  if (length(bad_taxa) > 1){
    otu[bad_taxa, samples7] <- 0
    otu_table(ps) <- otu
    # Now they're 0
    otu_table(ps)[bad_taxa, samples7]
  }
  
  #remove control samples
  ps.new = subset_samples(ps, !SampleID %in% c("Blank.Snow1", "Blank.Snow2", "Blank.Ice1", "Blank.NC1", "Blank.NC2", "Blank.NC3", "Blank.NC4", "Blank.NC5", "Blank.Cryo", "Blank.E", "Blank.Ice2", "Blank.1", "Blank.2", "Blank.3", "Blank.A", "Blank.B", "Blank.C", "Blank.D", "Blank.4", "Blank.5", "Blank.6", "Blank.PCR.NC", "Blank.15", "Blank.19", "Blank.3a", "Blank.41"))
  
  #remove taxa with zero counts
  ps.new.rm <- prune_taxa(taxa_sums(ps.new) >= 1, ps.new) 
  
  #remove samples with zero counts
  ps.new.rm <- prune_samples(sample_sums(ps.new.rm)>=1, ps.new.rm)
  
  return(ps.new.rm)
  
}

ps.pro3 = remove_all_contaminants(ps.pro2)
ps.euk.mm3 = remove_all_contaminants(ps.euk.mm2)
ps.euk.mm.rm3 = remove_all_contaminants(ps.euk.mm.rm2)
ps.mm3 = remove_all_contaminants(ps.mm2)

#Export phyloseq objects with contaminants removed.
saveRDS(ps.pro3, "../results/16S-ps-decontam.rds")
saveRDS(ps.euk.mm3, "../results/18S-all-ps-decontam.rds")
saveRDS(ps.euk.mm.rm3, "../results/18S-no-mm-ps-decontam.rds")
saveRDS(ps.mm3, "../results/18S-mm-only-ps-decontam.rds")