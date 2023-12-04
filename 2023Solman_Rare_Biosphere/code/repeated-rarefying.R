# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Generate rarefaction curves to identify library size to rarefy to
# 4. Multiple iteration of rarefying libraries
# 5. Average rarefied tables
# 6. Generate new phyloseq object
# 7. Check sampling depth
# 8. Compare number of samples and species richness of average rarefied tables to raw data
# 9. Export data

# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
#BiocManager::install("escamero/mirlyn")
library(mirlyn)
library(ggplot2)

#number of repeats
num_rep = 500

# 2. Import data
pro <- readRDS("../results/16S-phylo-object.rds") 
euk <- readRDS("../results/18S-phylo-object.rds")

# 3. Generate rarefaction curves to identify library size to rarefy to
Rarefy_whole_rep_example <- rarefy_whole_rep(pro, rep = 2)

#prokaryotes
#Visualization of rarefaction curve
P1 <- mirlyn::rarecurve(Rarefy_whole_rep_example, sample = "Sample") 

P1 = P1 + 
  theme(legend.position="none") +
  scale_x_continuous(labels = scales::comma)

pdf("../results/16S-rarefaction-curve.pdf")
print(P1)
dev.off()

#eukaryotes
#lets look at the number of reads per sample
Rarefy_whole_rep_example <- rarefy_whole_rep(euk, rep = 2)

#Visualization of rarefaction curve
P2 <- mirlyn::rarecurve(Rarefy_whole_rep_example, sample = "Sample") 

P2 = P2 + theme(legend.position="none") +
  scale_x_continuous(labels = scales::comma)

pdf("../results/18S-rarefaction-curve.pdf")
print(P2)
dev.off()

#look at the library sizes
sort(sample_sums(pro))
sort(sample_sums(euk))

# 4. Multiple iteration of rarefying libraries
pro_lib = 5730
euk_lib = 2673 

# Repeatedly rarefies to specified library size n times
pro_mirl_object <- mirl(pro, libsize = pro_lib, rep = num_rep, set.seed = 666) 
euk_mirl_object <- mirl(euk, libsize = euk_lib, rep = num_rep, set.seed = 666) 

#Single rarefaction for comparison
pro_mirl_object_1 <- mirl(pro, libsize = pro_lib, rep = 1, set.seed = 666) 
euk_mirl_object_1 <- mirl(euk, libsize = euk_lib, rep = 1, set.seed = 666)

# 5. Average rarefied tables

average_tab <- function(mirl_object){
  
  count_tab = data.frame(otu_table(mirl_object[[1]]), check.names = FALSE)
  for (i in 2:length(mirl_object)){
    #extract the count table
    count_tab2 = data.frame(otu_table(mirl_object[[i]]), check.names = FALSE)
    #sum tables together
    count_tab = count_tab + count_tab2
  }
  final.tab = count_tab/num_rep
  final.tab2 = ceiling(final.tab)
  
  return(final.tab2)
}

pro.av.tab = average_tab(pro_mirl_object)
euk.av.tab = average_tab(euk_mirl_object)

# 6. Generate new phyloseq objects

pro_ASV = otu_table(as.matrix(pro.av.tab), taxa_are_rows = TRUE)
pro_TAX = tax_table(pro)
pro_META = sample_data(pro)
pro_TREE = phy_tree(pro)

#make into phyloseq object
pro_ps <- phyloseq(pro_ASV, pro_TAX, pro_META, pro_TREE)

euk_ASV = otu_table(as.matrix(euk.av.tab), taxa_are_rows = TRUE)
euk_TAX = tax_table(euk)
euk_META = sample_data(euk)
euk_TREE = phy_tree(euk)

#make into phyloseq object
euk_ps <- phyloseq(euk_ASV, euk_TAX, euk_META, euk_TREE)


# 7. Check sampling depth
pro.df = data.table::data.table(as(sample_data(pro_ps), "data.frame"),
                                TotalReads = sample_sums(pro_ps), keep.rownames = TRUE)

pro_p<-ggplot(data=pro.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")


pdf("../results/16S-sample-depth.pdf")
print(pro_p)
dev.off()

euk.df = data.table::data.table(as(sample_data(euk_ps), "data.frame"),
                                TotalReads = sample_sums(euk_ps), keep.rownames = TRUE)

euk_p<-ggplot(data=euk.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

pdf("../results/18S-sample-depth.pdf")
print(euk_p)
dev.off()

# 8. Compare number of samples and species richness of average rarefied tables to raw data

#remove taxa or samples with zero counts
final.pro = filter_taxa(pro_ps, function(x) sum(x) > 0, TRUE)
final.euk = filter_taxa(euk_ps, function(x) sum(x) > 0, TRUE)
final.pro = prune_samples(sample_sums(final.pro)>0, final.pro)
final.euk = prune_samples(sample_sums(final.euk)>0, final.euk)

#And from the one-shot rarefied data
final.pro1 = filter_taxa(pro_mirl_object_1[[1]], function(x) sum(x) > 0, TRUE)
final.euk1 = filter_taxa(euk_mirl_object_1[[1]], function(x) sum(x) > 0, TRUE)
final.pro1 = prune_samples(sample_sums(final.pro1)>0, final.pro1)
final.euk1 = prune_samples(sample_sums(final.euk1)>0, final.euk1)

#output results to text document
sink("../results/repeated-rarefaction.txt", type="output")
writeLines("===============================================================
REPEATED RAREFACTION
===============================================================")
writeLines("Prokaryotes rarefied to:")
pro_lib
writeLines("Eukaryotes rarefied to:")
euk_lib
writeLines("Number of prokaryote ASVs before:")
ntaxa(pro)
writeLines("Number of prokaryote ASVs after repeated rarefaction:")
ntaxa(final.pro)
writeLines("Number of prokaryote ASVs after one-shot rarefaction:")
ntaxa(final.pro1)
writeLines("Number of prokaryote samples before:")
nsamples(pro)
writeLines("Number of prokaryote samples after repeated rarefaction:")
nsamples(final.pro)
writeLines("Number of prokaryote samples after one-shot rarefaction:")
nsamples(final.pro1)
writeLines("Number of eukaryote ASVs before:")
ntaxa(euk)
writeLines("Number of eukaryote ASVs after repeated rarefaction:")
ntaxa(final.euk)
writeLines("Number of eukaryote ASVs after one-shot rarefaction:")
ntaxa(final.euk1)
writeLines("Number of eukaryote samples before:")
nsamples(euk)
writeLines("Number of eukaryote samples after repeated rarefaction:")
nsamples(final.euk)
writeLines("Number of eukaryote samples after one-shot rarefaction:")
nsamples(final.euk1)
sink()


# 9. Export data
saveRDS(final.pro, "../results/16S-phylo-object-rarefied.rds")
saveRDS(final.euk, "../results/18S-phylo-object-rarefied.rds")
