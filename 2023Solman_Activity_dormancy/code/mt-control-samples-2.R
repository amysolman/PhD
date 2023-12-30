
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
pro <- readRDS("../results/mt-dna-16S-phylo-object-rarefied.rds")
euk <- readRDS("../results/mt-dna-18S-phylo-object-micro-keep.rds")

#What is in our blanks?
cont = c("Blank1", "Blank2", "Blank3", "Blank4", "Blank5", "Blank6", "Blank7", "Blank8")

#extract the count table with these samples
#get count and tax tables from our phyloseq obejcts
pro.count.filt = data.frame(otu_table(pro))
euk.count.filt = data.frame(otu_table(euk))
pro.tax.filt = data.frame(tax_table(pro))
euk.tax.filt = data.frame(tax_table(euk))

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

write.csv(cont.tax.out, "../results/mt-dna-control-sample-ASV-tax.csv")

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
write.csv(euk_contam_tab, "../results/mt-dna-18S-contaminant-asvs.csv")

pro_decontam = names(pro.count.filt) %in% cont #get local vector of control samples
pro_contam_df <- isContaminant(t(pro.count.filt), neg=pro_decontam) # Run isContaminant, while transforming the matrix with t()
table(pro_contam_df$contaminant) # What do our results look like?
pro_contam_asvs <- row.names(pro_contam_df[pro_contam_df$contaminant == TRUE, ]) # Make a vector holding contaminant IDs
pro_contam_tab <- pro.tax.filt[rownames(pro.tax.filt) %in% pro_contam_asvs, ] # Check the taxonomy of these contaminants
write.csv(pro_contam_tab, "../results/mt-dna-16S-contaminant-asvs.csv")