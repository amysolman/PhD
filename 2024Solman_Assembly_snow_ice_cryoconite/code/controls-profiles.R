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
library(fantaxtic)

#load data
metadata <- read.csv(file="../data/metadata.csv", sep=",") #Metadata
ps.pro = readRDS("../results/16S-phylo-object.rds")
ps.euk.mm.rm = readRDS("../results/18S-phylo-object-micro-remove.rds")
ps.euk.mm = readRDS("../results/18S-phylo-object-micro-keep.rds")
ps.mm = readRDS("../results/18S-phylo-object-micro-only.rds")

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

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

#####################################################################################################################################

#GET TAXONOMY OF ASVs IN CONTROL SAMPLES

#What is in our blanks?
cont = metadata$SampleID[1:26]

#extract the count table with these samples
#get count and tax tables from our phyloseq obejcts
pro.count.filt = data.frame(otu_table(ps.pro))
euk.count.filt = data.frame(otu_table(ps.euk.mm.rm))
pro.tax.filt = data.frame(tax_table(ps.pro))
euk.tax.filt = data.frame(tax_table(ps.euk.mm.rm))

#only keep the control samples
pro.count.blanks.only <- pro.count.filt[,names(pro.count.filt) %in% cont]
euk.count.blanks.only <- euk.count.filt[,names(euk.count.filt) %in% cont]

#how many reads per blank?
tab.sum.blank.reads = as.data.frame(rbind(colSums(pro.count.blanks.only), colSums(euk.count.blanks.only)))
tab.sum.blank.reads = data.frame(cbind(c("16S", "18S"), tab.sum.blank.reads))
names(tab.sum.blank.reads) = c("Amplicon", cont) 

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

###Table with dataset, sample, number of reads, number of asvs
pro.samp = names(pro.cont)
pro.asvs = colSums(pro.cont)
pro.conts = colSums(pro.count.blanks.only)
pro.cont.df = data.frame(Sample=pro.samp, ASVs16S=pro.asvs, ReadCounts16S=pro.conts)

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

###Table with dataset, sample, number of reads, number of asvs
euk.samp = names(euk.cont)
euk.asvs = colSums(euk.cont)
euk.conts = colSums(euk.count.blanks.only)
euk.cont.df = data.frame(Sample=euk.samp, ASVs18S=euk.asvs, ReadCounts18S=euk.conts)

#combine tables
pro.euk.join = full_join(pro.cont.df, euk.cont.df, by = c("Sample"))

#combine dataframes to export for checking taxa
cont.tax.out = rbind(pro.cont.tax.m, euk.cont.tax.m)

write.csv(cont.tax.out, "../results/control-sample-ASV-tax.csv")

#####################################################################################################################################

###COMMUNITY BAR PLOTS

#16S Contaminants
ps = subset_samples(ps.pro, Habitat == "Control")
#remove samples of ASVs with no reads
ps = prune_samples(sample_sums(ps)>0, ps)
ps = filter_taxa(ps, function(x) sum(x) > 0, TRUE)

#which are our most abundant phyla?
top.phy <- fantaxtic::top_taxa(ps, tax_level = "Phylum", n_taxa = 15)
#which phyla are we keeping?
phy.keep = unique(top.phy[[2]]$Phylum)
phy.keep

#Modify tax names so it looks better when we plot them
x = data.frame(tax_table(ps)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))

#replace anything that says NA with Genus Unknown
y = x
y$Phylum <- gsub("P.NA", "(Genus Unknown)", y$Phylum)
y$Class <- gsub("P.NA C.NA", "(Genus Unknown)", y$Class)
y$Class <- gsub("C.NA", "(Genus Unknown)", y$Class)
y$Order <- gsub("P.NA C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("O.NA", "(Genus Unknown)", y$Order)
y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("F.NA", "(Genus Unknown)", y$Family)
y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)

tax_table(ps) = as.matrix(y)

#get our data
data <- ps %>%
  #subset_samples(Habitat == "Control") %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum - we'll use this dataframe to assign colours
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

#remove anything not assigned at the phylum level
#out.col = pro_data_long.no.0[!is.na(pro_data_long.no.0$Phylum),]
  
for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = pro.cols[pro.cols$Phylum == unique(mini.df$Phylum),]$Colour1
  col2 = pro.cols[pro.cols$Phylum == unique(mini.df$Phylum),]$Colour2
  
  #get colour function for that phylum
  cols.fun = colorRampPalette(c(col1, col2))
  
  #get colours for each class
  save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#add phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Genus2 = paste0(data2plot$Phylum, ": ", data2plot$Genus)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Genus2[i] %in% out.col$Genus){
    #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Genus2[i] %in% out.col$Genus){ 
    #else they are just given the phyla name with the genus name
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    
  }
}

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "genus" with Other
    data2plot$Genus[i] <- "Other"
  }
}

#add black to genera classified as "Other"
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))

#add colours to plotting dataframe
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2                           

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Genus as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Genus)), ]
levs = gen.df$Genus
levs = levs[levs != "Other"]
df.plot$Genus <- factor(df.plot$Genus, levels = out.col2$Genus)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Genus)

df.plot2 = 
  df.plot %>%
  group_by(Genus, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p1 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  #facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol = 5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p1

pdf("../results/prokaryote-control-profiles.pdf", width=35, height=10)
print(p1)
dev.off()


#18S Contaminants

ps = subset_samples(ps.euk.mm, Habitat == "Control")
#remove samples of ASVs with no reads
ps = prune_samples(sample_sums(ps)>0, ps)
ps = filter_taxa(ps, function(x) sum(x) > 0, TRUE)

#which are our most abundant phyla?
top.phy <- fantaxtic::top_taxa(ps, tax_level = "Phylum", n_taxa = 15)
#which phyla are we keeping?
phy.keep = unique(top.phy[[2]]$Phylum)
phy.keep

#Modify tax names so it looks better when we plot them
x = data.frame(tax_table(ps)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))

#replace anything that says NA with Genus Unknown
y = x
y$Phylum <- gsub("P.NA", "(Genus Unknown)", y$Phylum)
y$Class <- gsub("P.NA C.NA", "(Genus Unknown)", y$Class)
y$Class <- gsub("C.NA", "(Genus Unknown)", y$Class)
y$Order <- gsub("P.NA C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("O.NA", "(Genus Unknown)", y$Order)
y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("F.NA", "(Genus Unknown)", y$Family)
y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)

tax_table(ps) = as.matrix(y)

#get our data
data <- ps %>%
  #subset_samples(Habitat == "Control") %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum - we'll use this dataframe to assign colours
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

#remove anything not assigned at the phylum level
#out.col = pro_data_long.no.0[!is.na(pro_data_long.no.0$Phylum),]

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = euk.cols[euk.cols$Phylum == unique(mini.df$Phylum),]$Colour1
  col2 = euk.cols[euk.cols$Phylum == unique(mini.df$Phylum),]$Colour2
  
  #get colour function for that phylum
  cols.fun = colorRampPalette(c(col1, col2))
  
  #get colours for each class
  save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#add phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Genus2 = paste0(data2plot$Phylum, ": ", data2plot$Genus)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Genus2[i] %in% out.col$Genus){
    #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Genus2[i] %in% out.col$Genus){ 
    #else they are just given the phyla name with the genus name
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    
  }
}

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "genus" with Other
    data2plot$Genus[i] <- "Other"
  }
}

#add black to genera classified as "Other"
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))

#add colours to plotting dataframe
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2                           

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Genus as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Genus)), ]
levs = gen.df$Genus
levs = levs[levs != "Other"]
df.plot$Genus <- factor(df.plot$Genus, levels = out.col2$Genus)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Genus)

df.plot2 = 
  df.plot %>%
  group_by(Genus, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p2 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  #facet_grid(cols = vars(SubPole), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol = 5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=18),
        strip.text.x = element_text(size = 20))

p2

pdf("../results/eukaryote-control-profiles.pdf", width=35, height=10)
print(p2)
dev.off()