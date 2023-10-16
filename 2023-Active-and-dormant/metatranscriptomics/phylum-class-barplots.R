# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
#library(vegan)
#source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
library(dplyr)
library(cowplot)
library(reshape) 
library(funrar) #for make relative
library(stringr) #mutate function
library(gridExtra) #for exporting as pdf
library(scales)
# devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic")
library(tidyverse)
library(patchwork)
library(ggpubr)
library(plyr)
library(zoo) 

# 2. Import data

#prokaryotes
pro.DNA <- readRDS("../results/16S-phylo-object.rds") 
pro.RNA <- readRDS("../results/mt-pro-phylo-object.rds") 

#eukaryotes without micrometazoans
euk.DNA <- readRDS("../results/18S-phylo-object-micro-keep.rds") 
euk.RNA <- readRDS("../results/mt-euk-phylo-object.rds") 


#16S Dataset

#combine the prokaryote DNA and RNA phyloseq objects into one

#get otu tables
tab1 = data.frame(otu_table(pro.DNA))
tab2 = data.frame(otu_table(pro.RNA))

#rename the samples
names(tab1) = paste0(names(tab1), ".DNA")
names(tab2) = paste0(names(tab2), ".RNA")

#combine the otu tables
total.tab = rbind.fill(tab1, tab2)
rownames(total.tab) = c(rownames(tab1), rownames(tab2)) #give the ASVs/NTUs the right names
total.tab[is.na(total.tab)] <- 0 #make NA into zero

#get taxa tables
tax.tab1 = data.frame(tax_table(pro.DNA))
tax.tab2 = data.frame(tax_table(pro.RNA))

#combine the taxa tables (excluding the columns they don't share)
total.tax = rbind(tax.tab1[, !(names(tax.tab1) %in% 'Species')], 
                  tax.tab2[, !(names(tax.tab2) %in% c('Major_clade', 'Kingdom'))])

#remove underscores in taxa names to make them the same between DNA and RNA datasets
total.tax$Class <-gsub("_"," ",as.character(total.tax$Class))
total.tax = as.matrix(total.tax)

#get metadata
meta1 = data.frame(sample_data(pro.DNA))
meta2 = data.frame(sample_data(pro.RNA))

#add new sample names
rownames(meta1) = paste0(rownames(meta1), ".DNA")
rownames(meta2) = paste0(rownames(meta2), ".RNA")

#add data type column
meta1$Data = "DNA"
meta2$Data = "RNA"

#remove random RNA metadata column
meta2 = subset(meta2, select = -c(SampleID18S) )

#combine metadata
total.meta = rbind(meta1, meta2)

#combine into one phyloseq object
pro = phyloseq(otu_table(total.tab, taxa_are_rows = TRUE), tax_table(total.tax), sample_data(total.meta))

#remove control samples
ps.pro2 = subset_samples(pro, Habitat != "Control")

#remove taxa with zero counts
ps.pro2 = filter_taxa(ps.pro2, function(x) sum(x) > 0, TRUE)

#modify tax names
df = data.frame(tax_table(ps.pro2))
for (i in 1:nrow(df)){
  data = unlist(as.vector(df[i,]))
  x = na.locf(data)
  y = names(data[is.na(data)]) #which columns had NA?
  for (j in 1:length(y)){
    col = y[j]
    x[col] = paste0(x[col], " (",  y[j], " Unknown)")
  }
  df[i,] = x
  
}
tax_table(ps.pro2) = as.matrix(df)

#find the 10 most abundant Phyla
top.phy <- top_taxa(ps.pro2, 
                    tax_level = "Phylum", 
                    n_taxa = 10,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#find the three most abundant Class within each phyla
#get our data
data <-
  ps.pro2 %>%
  tax_glom("Class") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant Class by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Class) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#get colours for plotting
col.list = c("orange", "darkorange4",
             "blue", "darkblue",
             "darkolivegreen1", "darkolivegreen",
             "tan", "tan4", 
             "turquoise", "turquoise4",
             "red", "darkred",
             "green", "darkgreen",
             "violetred1", "violetred4",
             "yellow", "yellow4", 
             "pink", "deeppink4", 
             "purple","darkorchid4", 
             "coral", "coral4", 
             "aquamarine", "aquamarine4",
             "brown1", "brown4",
             "antiquewhite", "antiquewhite4",
             "chocolate", "chocolate4",
             "deepskyblue", "deepskyblue4",
             "darkgoldenrod1", "darkgoldenrod4",
             "cyan", "cyan4")

#sanity check
length(unique(data$Phylum))

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
num4cols = 1
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #get colour function for that phylum
  cols = colorRampPalette(c(col.list[num4cols], col.list[num4cols+1]))
  
  #get colours for each Class
  save.cols = c(save.cols, cols(nrow(mini.df)+1))
  
  #add nums for moving along our colour list
  num4cols = num4cols + 2
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Class=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#out phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))

#remove unwanted genera and replace with other
for (i in 1:nrow(data2plot)){
  if(!data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": Other")
  } else if (data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Class[i])
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))
sort(unique(data2plot$Class))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
    data2plot$Class[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))
sort(unique(data2plot$Class))

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Class="Other", col="#000000"))
#check we have the same values
#length(intersect(out.col2$Class, data2plot$Class)) == length(unique(data2plot$Class)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat2, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Data <- factor(df$Data, levels = c("DNA", "RNA"))

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Class as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Class)), ]
levs = gen.df$Class
levs = levs[levs != "Other"]
df.plot$Class <- factor(df.plot$Class, levels = out.col2$Class)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Class)

df.plot2 = 
  df.plot %>%
  group_by(Class, Habitat, Data, SampleID) %>%
  dplyr::summarise(across(c(Abundance), sum))

pro.p = ggplot(df.plot2, aes(fill=Class, y=Abundance, x=SampleID)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), rows = vars(Data), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y="Relative abundance of 16S rRNA genes (DNA) and gene transcripts (RNA)", x = "Sample")+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=4, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

pro.p


pdf("../results/prokaryote-community-plot.pdf", width=15, height=10)
print(pro.p)
dev.off()

###18S Data Set

#combine the eukaryote DNA and RNA phyloseq objects into one

#remove microfauna from the DNA data
mm.rm = subset_taxa(euk.DNA, !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa", "Arthropoda"))

#get otu tables
tab1 = data.frame(otu_table(mm.rm))
tab2 = data.frame(otu_table(euk.RNA))

#rename the samples
names(tab1) = paste0(names(tab1), ".DNA")
names(tab2) = paste0(names(tab2), ".RNA")

#combine the otu tables
total.tab = rbind.fill(tab1, tab2)
rownames(total.tab) = c(rownames(tab1), rownames(tab2)) #give the ASVs/NTUs the right names
total.tab[is.na(total.tab)] <- 0 #make NA into zero

#get taxa tables
tax.tab1 = data.frame(tax_table(euk.DNA))
tax.tab2 = data.frame(tax_table(euk.RNA))

#combine the taxa tables
total.tax = rbind(tax.tab1[, !(names(tax.tab1) %in% 'Species')], 
                  tax.tab2[, !(names(tax.tab2) %in% c('Major_clade', 'Kingdom'))])

#remove underscores in taxa names to make them the same between DNA and RNA datasets
total.tax$Class <-gsub("_"," ",as.character(total.tax$Class))
total.tax = as.matrix(total.tax)

#get metadata
meta1 = data.frame(sample_data(euk.DNA))
meta2 = data.frame(sample_data(euk.RNA))

#add new sample names
rownames(meta1) = paste0(rownames(meta1), ".DNA")
rownames(meta2) = paste0(rownames(meta2), ".RNA")

#add data type column
meta1$Data = "DNA"
meta2$Data = "RNA"

#combine metadata
total.meta = rbind(meta1, meta2)

#combine into one phyloseq object
euk = phyloseq(otu_table(total.tab, taxa_are_rows = TRUE), tax_table(total.tax), sample_data(total.meta))

#remove control samples
ps.euk2 = subset_samples(euk, Habitat != "Control")

#remove taxa with zero counts
ps.euk2 = filter_taxa(ps.euk2, function(x) sum(x) > 0, TRUE)

#modify tax names
df = data.frame(tax_table(ps.euk2))
for (i in 1:nrow(df)){
  data = unlist(as.vector(df[i,]))
  x = na.locf(data)
  y = names(data[is.na(data)]) #which columns had NA?
  for (j in 1:length(y)){
    col = y[j]
    x[col] = paste0(x[col], " (",  y[j], " Unknown)")
  }
  df[i,] = x
  
}

tax_table(ps.euk2) = as.matrix(df)

#find the 10 most abundant Phyla
top.phy <- top_taxa(ps.euk2, 
                    tax_level = "Phylum", 
                    n_taxa = 10,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#find the three most abundant Class within each phyla
#get our data
data <-
  ps.euk2 %>%
  tax_glom("Class") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant Class by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Class) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#get colours for plotting
col.list = c("orange", "darkorange4",
             "blue", "darkblue",
             "darkolivegreen1", "darkolivegreen",
             "tan", "tan4", 
             "turquoise", "turquoise4",
             "red", "darkred",
             "green", "darkgreen",
             "violetred1", "violetred4",
             "yellow", "yellow4", 
             "pink", "deeppink4", 
             "purple","darkorchid4", 
             "coral", "coral4", 
             "aquamarine", "aquamarine4",
             "brown1", "brown4",
             "antiquewhite", "antiquewhite4",
             "chocolate", "chocolate4",
             "deepskyblue", "deepskyblue4",
             "darkgoldenrod1", "darkgoldenrod4",
             "cyan", "cyan4")

#sanity check
length(unique(data$Phylum))

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
num4cols = 1
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #get colour function for that phylum
  cols = colorRampPalette(c(col.list[num4cols], col.list[num4cols+1]))
  
  #get colours for each Class
  save.cols = c(save.cols, cols(nrow(mini.df)+1))
  
  #add nums for moving along our colour list
  num4cols = num4cols + 2
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Class=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#out phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))

#remove unwanted genera and replace with other
for (i in 1:nrow(data2plot)){
  if(!data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": Other")
  } else if (data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Class[i])
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))
sort(unique(data2plot$Class))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
    data2plot$Class[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum))
length(unique(data2plot$Class))
sort(unique(data2plot$Class))

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Class="Other", col="#000000"))
#check we have the same values
#length(intersect(out.col2$Class, data2plot$Class)) == length(unique(data2plot$Class)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat2, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Data <- factor(df$Data, levels = c("DNA", "RNA"))

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Class as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Class)), ]
levs = gen.df$Class
levs = levs[levs != "Other"]
df.plot$Class <- factor(df.plot$Class, levels = out.col2$Class)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Class)

df.plot2 = 
  df.plot %>%
  group_by(Class, Habitat, Data, SampleID) %>%
  dplyr::summarise(across(c(Abundance), sum))

euk.p = ggplot(df.plot2, aes(fill=Class, y=Abundance, x=SampleID)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), rows = vars(Data), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y="Relative abundance of 18S rRNA genes (DNA) and gene transcripts (RNA)", x = "Sample")+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=4, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

euk.p


pdf("../results/eukaryote-community-plot.pdf", width=15, height=10)
print(euk.p)
dev.off()


#merge 16S and 18S results into one plot

p3 = pro.p / euk.p + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p3


pdf("../results/all-community-plot.pdf", width=18, height=20)
print(p3)
dev.off()


