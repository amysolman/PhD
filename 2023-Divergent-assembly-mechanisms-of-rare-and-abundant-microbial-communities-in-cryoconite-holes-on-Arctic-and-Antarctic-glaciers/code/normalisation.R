
# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
#BiocManager::install("escamero/mirlyn")
library(mirlyn)
library(ggplot2)

#number of repeats
num_rep = 500 #took about 15 mins at 200 reps, 35 mins at 500 reps


## Method

#Repeated rarefying was used to normalise the data. 
#Rareying is a commonly used noralisation technique that equalises library sizes 
#by randomly subsampling observed counts. One-shot rarefying provides only one 
#snapshot of the community at a smaller library size. This omits a subset of 
#observed sequences leading to the loss of ASVs and samples. Repeated rarefying and 
#calculating an average rarefied table can preseve a greater number of ASVs and sequences. 
#Rarefying was carried out `r num_rep` times. The average count of each ASV was calculated 
#with floats rounded up to generate the final average rarefied table. 
#Rarefaction curves and repeated rarefying was carried out used the mirlyn R package. 

## Results


# 2. Import data
pro <- readRDS("../results/16S-phylo-object.rds") 
euk <- readRDS("../results/18S-phylo-object.rds")

# 3. Generate rarefaction curves to identify library size to rarefy to
sort(sample_sums(pro))

Rarefy_whole_rep_example <- rarefy_whole_rep(pro, rep = 2)

#Visualization of rarefaction curve
P1 <- mirlyn::rarecurve(Rarefy_whole_rep_example, sample = "Sample") 

P1 = P1 + 
  theme(legend.position="none") +
  scale_x_continuous(labels = scales::comma)
P1 


sort(sample_sums(euk))

#lets look at the number of reads per sample
Rarefy_whole_rep_example <- rarefy_whole_rep(euk, rep = 2)

#Visualization of rarefaction curve
P2 <- mirlyn::rarecurve(Rarefy_whole_rep_example, sample = "Sample") 

P2 = P2 + theme(legend.position="none") +
  scale_x_continuous(labels = scales::comma)
P2

# 4. Multiple iteration of rarefying libraries
pro_lib = 5730 #or 6000 or 5534 or 6268 (so far 6000 and 6268 works best)
euk_lib = 2673 #or 3000 or 2556 or 2905 (so far 3000 and 2905 works best)

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

pro_p

euk.df = data.table::data.table(as(sample_data(euk_ps), "data.frame"),
                                TotalReads = sample_sums(euk_ps), keep.rownames = TRUE)

euk_p<-ggplot(data=euk.df, aes(x=as.factor(SampleID), y=TotalReads)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("SampleID")+
  ylab("Library Depth")

euk_p


# 7. Compare rarefied phyloseq objects to original phyloseq objects

#remove taxa or samples with zero counts
final.pro = filter_taxa(pro_ps, function(x) sum(x) > 0, TRUE)
final.euk = filter_taxa(euk_ps, function(x) sum(x) > 0, TRUE)
final.pro = prune_samples(sample_sums(final.pro)>=1, final.pro)
final.euk = prune_samples(sample_sums(final.euk)>=1, final.euk)

#And from the one-shot rarefied data
final.pro1 = filter_taxa(pro_mirl_object_1[[1]], function(x) sum(x) > 0, TRUE)
final.euk1 = filter_taxa(euk_mirl_object_1[[1]], function(x) sum(x) > 0, TRUE)
final.pro1 = prune_samples(sample_sums(final.pro1)>=1, final.pro1)
final.euk1 = prune_samples(sample_sums(final.euk1)>=1, final.euk1)

#Look at 16S data
pro #original
final.pro #repeated rarefied
final.pro1 #one-shot rarefied

#Look at 18S data
euk #original
final.euk #repeated rarefied
final.euk1 #one-shot rarefied

# 9. Export data

saveRDS(final.pro, "../results/16S-phylo-object-rarefied.rds")
saveRDS(final.euk, "../results/18S-phylo-object-rarefied.rds")