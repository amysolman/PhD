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

#with "contaminants" removed
ps.pro.d = readRDS("../results/16S-ps-decontam.rds")
ps.euk.mm.rm.d = readRDS("../results/18S-no-mm-ps-decontam.rds")
ps.euk.mm.d = readRDS("../results/18S-all-ps-decontam.rds")
ps.mm.d = readRDS("../results/18S-mm-only-ps-decontam.rds")

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
# check that they are expected in our environment and not human associated (see controls-high-abundance.R).
# 4. If they are expected then we will assume these are carry over from TRUE samples. 
# If they are human associated they will be manually removed from the dataset (see controls-remove.R)
# 5. Compare datasets before and after removing ASVs in step  using NMDS plots and ANOSIM (see controls-NMDS.R).

#####################################################################################################################################

#merge old and new phyloseq objects
merge_phylo <- function(ps1, ps2){
  
  #get count data
  asv.tab1 = data.frame(otu_table(ps1))
  sum(colSums(asv.tab1))
  summarize_phyloseq(ps1)[3]
  asv.tab2 = data.frame(otu_table(ps2))
  sum(colSums(asv.tab2))
  
  #get metadata
  meta1 = data.frame(sample_data(ps1))
  meta2 = data.frame(sample_data(ps2))
  
  #change sample names
  meta1$SampleID = paste0(meta1$SampleID, "_Before")
  meta2$SampleID = paste0(meta2$SampleID, "_After")
  
  #change sample names
  names(asv.tab1) = meta1$SampleID
  names(asv.tab2) = meta2$SampleID
  
  #add dataset variable to metadata
  meta1$Data = "Before"
  meta2$Data = "After"
  
  #merge asv tables
  asv.tab1$ID = rownames(asv.tab1)
  asv.tab2$ID = rownames(asv.tab2)
  
  asv.tab.full = full_join(asv.tab1, asv.tab2, by = "ID")
  asv.tab.full[is.na(asv.tab.full)] <- 0 #change NAs to 0
  
  #give ASV ID to rownames
  rownames(asv.tab.full) = asv.tab.full$ID
  #remove ID column
  asv.tab.full = asv.tab.full[,-which(names(asv.tab.full) %in% c("ID"))]
  
  colSums(asv.tab.full)
  
  #merge metadata
  meta.full = data.frame(rbind(meta1, meta2))
  rownames(meta.full) = meta.full$SampleID
  
  #make phyloseq object
  new.ps = phyloseq(otu_table(asv.tab.full, taxa_are_rows = TRUE), sample_data(meta.full))
  
  #remove taxa with zero counts
  ps.new.rm <- prune_taxa(taxa_sums(new.ps) >= 1, new.ps) 
  #remove samples with zero counts
  ps.new.rm <- prune_samples(sample_sums(ps.new.rm)>=1, ps.new.rm)
  
  sample_sums(ps.new.rm)
  
  
  return(ps.new.rm)
  
}


#run the functions

#NMDS plot
#get phyloseq with controls removed only
ps.pro.sub = subset_samples(ps.pro, Habitat != "Control")
#remove ASVs with zero counts
ps.pro.sub = prune_taxa(taxa_sums(ps.pro.sub) > 0, ps.pro.sub)
#combine together phyloseq with contam removed and phyloseq with contam
full.pro = merge_phylo(ps.pro.sub, ps.pro.d)
#nmds
set.seed(333)
pro.nmds = phyloseq::ordinate(full.pro, method="NMDS")
p1 = phyloseq::plot_ordination(full.pro, pro.nmds, type="samples", color="Data")+
  theme_bw()+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))
#view plot
print(p1)
#test for significant differences
comm = as.matrix(t(as.data.frame(otu_table(full.pro))))
ano.pro = anosim(comm, data.frame(sample_data(full.pro))$Data, distance = "bray", permutations = 9999)
ano.pro
#finish by saving the original phyloseq object with blanks removed
saveRDS(ps.pro.sub, "../results/16S-ps-controls-removed.rds")

#FULL EUKARYOTE DATASET
#NMDS plot 
#get phyloseq with controls removed only
ps.euk.sub.mm = subset_samples(ps.euk.mm, Habitat != "Control")
#remove ASVs with zero counts
ps.euk.sub.mm = prune_taxa(taxa_sums(ps.euk.sub.mm) > 0, ps.euk.sub.mm)
#combine together phyloseq with contam removed and phyloseq with contam
full.euk.mm = merge_phylo(ps.euk.sub.mm, ps.euk.mm.d)
#nmds
euk.nmds.mm = phyloseq::ordinate(full.euk.mm, method="NMDS")
#plot
p2 = phyloseq::plot_ordination(full.euk.mm, euk.nmds.mm, type="samples", color="Data")+
  theme_bw()+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))
print(p2)

#test for significant differences
comm = as.matrix(t(as.data.frame(otu_table(full.euk.mm))))
ano.euk.mm = anosim(comm, data.frame(sample_data(full.euk.mm))$Data, distance = "bray", permutations = 9999)
ano.euk.mm

#finish by saving the original phyloseq object with blanks removed
saveRDS(ps.euk.sub.mm, "../results/18S-ps-mm-keep-controls-removed.rds")

#MICRO EUKARYOTES ONLY
#NMDS plot 
#get phyloseq with controls removed only
ps.euk.sub.mm.rm = subset_samples(ps.euk.mm.rm, Habitat != "Control")
#remove ASVs with zero counts
ps.euk.sub.mm.rm = prune_taxa(taxa_sums(ps.euk.sub.mm.rm) > 0, ps.euk.sub.mm.rm)
#combine together phyloseq with contam removed and phyloseq with contam
full.euk.mm.rm = merge_phylo(ps.euk.sub.mm.rm, ps.euk.mm.rm.d)
#nmds
euk.nmds.mm.rm = phyloseq::ordinate(full.euk.mm.rm, method="NMDS")
#plot
p3 = phyloseq::plot_ordination(full.euk.mm.rm, euk.nmds.mm.rm, type="samples", color="Data")+
  theme_bw()+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))
print(p3)

#test for significant differences
comm = as.matrix(t(as.data.frame(otu_table(full.euk.mm.rm))))
ano.euk.mm.rm = anosim(comm, data.frame(sample_data(full.euk.mm.rm))$Data, distance = "bray", permutations = 9999)
ano.euk.mm.rm

#finish by saving the original phyloseq object with blanks removed
saveRDS(ps.euk.sub.mm.rm, "../results/18S-ps-no-mm-controls-removed.rds")

#MICROMETAZOANS ONLY

#NMDS plot 
#get phyloseq with controls removed only
ps.mm.sub = subset_samples(ps.mm, Habitat != "Control")
#remove ASVs with zero counts
ps.mm.sub = prune_taxa(taxa_sums(ps.mm.sub) > 0, ps.mm.sub)
#combine together phyloseq with contam removed and phyloseq with contam
full.mm = merge_phylo(ps.mm.sub, ps.mm.d)

#nmds
nmds.mm = phyloseq::ordinate(full.mm, method="NMDS")

#plot
p4 = phyloseq::plot_ordination(full.mm, nmds.mm, type="samples", color="Data")+
  theme_bw()+
  theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))
print(p4)

#test for significant differences
comm = as.matrix(t(as.data.frame(otu_table(full.mm))))
ano.mm = anosim(comm, data.frame(sample_data(full.mm))$Data, distance = "bray", permutations = 9999)
ano.mm

#finish by saving the original phyloseq object with blanks removed
saveRDS(ps.mm.sub, "../results/18S-ps-mm-only-controls-removed.rds")


#Make into one image
legend <- get_legend(
  p1 + 
    #guides(color = guide_legend(nrow = 1, override.aes = list(size = 10))) +
    theme(legend.position = "bottom",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"),
          legend.text = element_text(size=10), 
          legend.title = element_blank())+
    guides(colour = guide_legend(override.aes = list(size = 6, alpha=.8)))+
    scale_color_discrete(breaks=c('Before', 'After'))
)

p = plot_grid(p1 + theme(legend.position = "none"),
              p2 + theme(legend.position = "none"),
              rel_widths = c(1,1),
              labels=c("Prokayrote", "Eukaryote"))
p

final.p = plot_grid(p, legend, ncol = 1, rel_heights = c(1, .1))
final.p

pdf("../results/NMDS-before-after-contaminant-removal.pdf", width=10, height=6)
print(final.p)
dev.off()

###################################################################################################
###################################################################################################

#Report NMDS Stress and ANOVA Results
sink("../results/contaminants-nmds-stress-and-anova.txt", type="output")
writeLines("===============================================================
NMDS STRESS AND ANOVA RESULTS
===============================================================")
writeLines("NMDS stress values below 0.05 indicate a good fit, 0.05-0.1 are a fair fit, 0.1-0.2 is a poor fit, 0.2-0.3 is a failed fit.") 

writeLines("NMDS stress for prokaryote dataset with Bray-Curtis Dissimilarities:")
round(pro.nmds$stress,3)
writeLines("ANOVA results for differences between datasets before and after removing contaminants in the prokaryote dataset:")
ano.pro

writeLines("NMDS stress for total eukaryote dataset with Bray-Curtis Dissimilarities:")
round(euk.nmds.mm$stress,3)
writeLines("ANOVA results for differences between datasets before and after removing contaminants in the total eukaryote dataset:")
ano.euk.mm

writeLines("NMDS stress for micro-eukaryote dataset with Bray-Curtis Dissimilarities:")
round(euk.nmds.mm.rm$stress,3)
writeLines("ANOVA results for differences between datasets before and after removing contaminants in the micro-eukaryote dataset:")
ano.euk.mm.rm

writeLines("NMDS stress for micro-fauna dataset with Bray-Curtis Dissimilarities:")
round(nmds.mm$stress,3)
writeLines("ANOVA results for differences between datasets before and after removing contaminants in the micro-fauna dataset:")
ano.mm

sink()
###################################################################################################
###################################################################################################
