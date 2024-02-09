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
# 3. Find ASVs that are in blanks. Of those that show high read abundance in TRUE samples (>0.5% total reads) check that they are expected in our environment and not human associated. 
#If they are expected then we will assume these are carry over from TRUE samples. 
#If they are human associated they will be manually removed from the dataset (see controls-high-abundance.R).
# 4. For each remaining ASV in the negative controls, calculate mean sum of reads in negative controls and divide by the mean sum of reads across control and TRUE samples. 
#Exclude ASVs with >5% abundance in negative controls (see controls-remove.R).
# 5. Compare datasets before and after removing ASVs in step 3 using NMDS plots and ANOSIM (see controls-NMDS.R).

# ps = ps.pro

compare_neg_lib <- function(ps){
  
  df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  
  df$Habitat = factor(df$Habitat,
                      levels=c("Control", "Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  df$LibrarySize <- sample_sums(ps)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  p = ggplot(data=df, aes(x=Index, y=LibrarySize, color=Habitat)) + geom_point()+
    theme_bw()+
    theme(plot.margin = margin(1,0.5,0.5,0.5, "cm"))+
    ylim(0, 150000)
  
  return(p)
  
}

pro.p = compare_neg_lib(ps.pro)
pro.p

euk.p = compare_neg_lib(ps.euk.mm)
euk.p

#add legend 

legend <- get_legend(
  pro.p + 
    #guides(color = guide_legend(nrow = 1, override.aes = list(size = 10))) +
    theme(legend.position = "bottom",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"),
          legend.text = element_text(size=10), 
          legend.title = element_blank())+
    guides(colour = guide_legend(override.aes = list(size = 6, alpha=.8)))
)

p = plot_grid(pro.p + theme(legend.position = "none"),
              euk.p + theme(legend.position = "none", axis.title.y = element_blank()),
              rel_widths = c(1,1),
              labels=c("Prokayrote", "Eukaryote"))

final.p = plot_grid(p, legend, ncol = 1, rel_heights = c(1, .1))
final.p

pdf("../results/read-depth.pdf", width=10, height=6)
print(final.p)
dev.off()
