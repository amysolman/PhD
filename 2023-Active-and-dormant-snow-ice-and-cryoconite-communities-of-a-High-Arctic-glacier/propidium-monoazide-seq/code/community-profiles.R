
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

#prokaryotes
ps.pro <- readRDS("../results/16S-phylo-object-norm-rare.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-phylo-object-norm-rare.rds") 



#####16S RELATIVE ABUNDANCE BAR CHART##########

#remove control samples
ps.pro2 = subset_samples(ps.pro, Habitat != "Control")
#remove taxa with zero countrs
ps.pro2 = filter_taxa(ps.pro2, function(x) sum(x) > 0, TRUE)

#modify tax names
x = data.frame(tax_table(ps.pro2)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Kingdom, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))
tax_table(ps.pro2) = as.matrix(x)

#find the 10 most abundant Phyla
top.phy <- top_taxa(ps.pro2, 
                    tax_level = "Phylum", 
                    n_taxa = 10,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#find the three most abundant genus within each phyla
#get our data
data <-
  ps.pro2 %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#check our data is correct with relative abundances
sum(data$Abundance) #should equal 96
sum(data[data$Sample == "251a",]$Abundance) #should equal 1

#get 3 most abundant genus by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#get colours for plotting
col.list = c("orange", "darkorange4",
             "blue", "darkblue",
             "tan", "tan4", 
             "turquoise", "turquoise4",
             "red", "darkred",
             "antiquewhite", "antiquewhite4",
             "green", "darkgreen",
             "darkgoldenrod1", "darkgoldenrod4",
             "violetred1", "violetred4", 
             "cyan", "cyan4",
             "yellow", "yellow4", 
             "pink", "deeppink4", 
             "purple","darkorchid4", 
             "coral", "coral4", 
             "aquamarine", "aquamarine4",
             "brown1", "brown4",
             "darkolivegreen1", "darkolivegreen",
             "chocolate", "chocolate4",
             "deepskyblue", "deepskyblue4")

#sanity check
length(unique(data$Phylum)) #40

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
num4cols = 1
#insert_rows = c(4,8,12,16,20,24,28,32,36,40,42)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #get colour function for that phylum
  cols = colorRampPalette(c(col.list[num4cols], col.list[num4cols+1]))
  
  #get colours for each genus
  save.cols = c(save.cols, cols(nrow(mini.df)+1))
  
  #add nums for moving along our colour list
  num4cols = num4cols + 2
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#out phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
length(unique(data2plot$Phylum)) #40
length(unique(data2plot$Genus)) #795

#remove unwanted genera and replace with other
for (i in 1:nrow(data2plot)){
  if(!data2plot$Genus[i] %in% out$Genus){
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
  } else if (data2plot$Genus[i] %in% out$Genus){
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
  }
}

#sanity check
length(unique(data2plot$Phylum)) #40
length(unique(data2plot$Genus)) #81
sort(unique(data2plot$Genus))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
    data2plot$Genus[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum)) #40
length(unique(data2plot$Genus)) #57
sort(unique(data2plot$Genus))

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Genus, data2plot$Genus)) == length(unique(data2plot$Genus)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Treatment <- factor(df$Treatment, levels = c("No-PMA", "PMA"))
df$HabTreat <- paste0(df$Habitat, " (", df$Treatment, ")")
df$HabTreat <- factor(df$HabTreat, levels = c("Snow (No-PMA)", "Snow (PMA)", "Spring Ice (No-PMA)", 
                                              "Spring Ice (PMA)", "Summer Ice (No-PMA)", "Summer Ice (PMA)",
                                              "Cryoconite (No-PMA)", "Cryoconite (PMA)"))

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
  group_by(Genus, HabTreat, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p1 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(HabTreat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p1

pdf("../results/prokaryote-community-plot-rare.pdf", width=15, height=10)
print(p1)
dev.off()

#save colours
write.csv(gen.df, "../results/prokaryote-colours-for-plotting-rare.csv")

#####18S RELATIVE ABUNDANCE BAR CHART##########

#remove control samples
ps.euk2 = subset_samples(ps.euk, Habitat != "Control")
#remove taxa with zero countrs
ps.euk2 = filter_taxa(ps.euk2, function(x) sum(x) > 0, TRUE)

#colours for plotting

#modify tax names
x = data.frame(tax_table(ps.euk2)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Kingdom, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))
tax_table(ps.euk2) = as.matrix(x)

#find the 10 most abundant Phyla
top.phy <- top_taxa(ps.euk2, 
                    tax_level = "Phylum", 
                    n_taxa = 10,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#find the three most abundant genus within each phyla
#get our data
data <-
  ps.euk2 %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#check our data is correct with relative abundances
sum(data$Abundance) #should equal 96
sum(data[data$Sample == "251a",]$Abundance) #should equal 1

#get 3 most abundant genus by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)


#get colours for plotting
col.list = c("antiquewhite", "antiquewhite4",
             "orange", "darkorange4",
             "blue", "darkblue", 
             "violetred1", "violetred4", 
             "green", "darkgreen",
             "red", "darkred",
             "purple","darkorchid4", 
             "turquoise", "turquoise4",
             "yellow", "yellow4", 
             "pink", "deeppink4", 
             "chocolate", "chocolate4",
             "deepskyblue", "deepskyblue4",
             "brown1", "brown4",
             "darkolivegreen1", "darkolivegreen",
             "tan", "tan4", 
             "coral", "coral4", 
             "aquamarine", "aquamarine4",
             "darkgoldenrod1", "darkgoldenrod4",
             "cyan", "cyan4")


#sanity check
length(unique(data$Phylum)) #36

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
num4cols = 1
#insert_rows = c(4,8,12,16,20,24,28,32,36,40,42)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #get colour function for that phylum
  cols = colorRampPalette(c(col.list[num4cols], col.list[num4cols+1]))
  
  #get colours for each genus
  save.cols = c(save.cols, cols(nrow(mini.df)+1))
  
  #add nums for moving along our colour list
  num4cols = num4cols + 2
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#out phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
length(unique(data2plot$Phylum)) #36
length(unique(data2plot$Genus)) #253

#remove unwanted genera and replace with other
for (i in 1:nrow(data2plot)){
  if(!data2plot$Genus[i] %in% out$Genus){
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
  } else if (data2plot$Genus[i] %in% out$Genus){
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
  }
}

#sanity check
length(unique(data2plot$Phylum)) #36
length(unique(data2plot$Genus)) #68
sort(unique(data2plot$Genus))

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
    data2plot$Genus[i] <- "Other"
  }
}

#sanity check
length(unique(data2plot$Phylum)) #36
length(unique(data2plot$Genus)) #47
sort(unique(data2plot$Genus))

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Genus, data2plot$Genus)) == length(unique(data2plot$Genus)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Treatment <- factor(df$Treatment, levels = c("No-PMA", "PMA"))
df$HabTreat <- paste0(df$Habitat, " (", df$Treatment, ")")
df$HabTreat <- factor(df$HabTreat, levels = c("Snow (No-PMA)", "Snow (PMA)", "Spring Ice (No-PMA)", 
                                              "Spring Ice (PMA)", "Summer Ice (No-PMA)", "Summer Ice (PMA)",
                                              "Cryoconite (No-PMA)", "Cryoconite (PMA)"))

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
  group_by(Genus, HabTreat, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p2 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(HabTreat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p2

pdf("../results/eukaryote-community-plot.pdf", width=15, height=10)
print(p2)
dev.off()

#save colours
write.csv(gen.df, "../results/eukaryote-colours-for-plotting-rare.csv")


###PUT THE TWO PLOTS TOGETHER

p3 = p1 / p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p3

pdf("../results/all-community-plot-rare.pdf", width=18, height=20)
print(p3)
dev.off()

#Community structure numbers

#How many ASVs
ps.pro2
ps.euk2

#16S Snow No-PMA
sub = subset_samples(ps.pro2, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Snow PMA
sub = subset_samples(ps.pro2, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Spring Ice No-PMA
sub = subset_samples(ps.pro2, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Spring Ice PMA
sub = subset_samples(ps.pro2, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Summer Ice No-PMA
sub = subset_samples(ps.pro2, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Summer Ice PMA
sub = subset_samples(ps.pro2, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Cryoconite No-PMA
sub = subset_samples(ps.pro2, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S Cryoconite PMA
sub = subset_samples(ps.pro2, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Snow No-PMA
sub = subset_samples(ps.euk2, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Snow PMA
sub = subset_samples(ps.euk2, Habitat == "Snow")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Spring Ice No-PMA
sub = subset_samples(ps.euk2, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Spring Ice PMA
sub = subset_samples(ps.euk2, Habitat == "Spring Ice")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Summer Ice No-PMA
sub = subset_samples(ps.euk2, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Summer Ice PMA
sub = subset_samples(ps.euk2, Habitat == "Summer Ice")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Cryoconite No-PMA
sub = subset_samples(ps.euk2, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "No-PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S Cryoconite PMA
sub = subset_samples(ps.euk2, Habitat == "Cryoconite")
sub = subset_samples(sub, Treatment == "PMA")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#any archeae and if so which communities are they in? = subset_taxa(p.ar.a, Kingdom == "Archaea")
arch.sub = subset_taxa(ps.pro2, Kingdom == "Archaea")
#remove samples with zero counts
arch.sub1 = prune_samples(sample_sums(arch.sub) >0, arch.sub)
#remove taxa with zero counts
arch.sub2 = filter_taxa(arch.sub1, function(x) sum(x) > 0, TRUE)

x = data.frame(tax_table(arch.sub2))
sort(table(x$Family), decreasing = TRUE)
arch.counts = data.frame(otu_table(arch.sub2))
sort(colSums(arch.counts), decreasing = TRUE)


####Alpha Diversity

#16S
pro.alpha = estimate_richness(ps.pro2)
pro.alpha$Sample = sample_names(ps.pro2)
pro.alpha$Treatment = data.frame(sample_data(ps.pro2))$Treatment
pro.alpha$Habitat = data.frame(sample_data(ps.pro2))$Habitat
pro.alpha$Habitat = factor(pro.alpha$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

p1 <- ggboxplot(pro.alpha, x = "Treatment", y = "Observed",
                color = "Treatment", palette = "jco",
                add = "jitter")+
  facet_grid(cols = vars(Habitat))+
  theme(legend.position = "none")
p1
#  Add p-value
p1 = p1  + stat_compare_means()
p1


#18S
euk.alpha = estimate_richness(ps.euk2)
euk.alpha$Sample = sample_names(ps.euk2)
euk.alpha$Treatment = data.frame(sample_data(ps.euk2))$Treatment
euk.alpha$Habitat = data.frame(sample_data(ps.euk2))$Habitat
euk.alpha$Habitat = factor(euk.alpha$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

p2 <- ggboxplot(euk.alpha, x = "Treatment", y = "Observed",
                color = "Treatment", palette = "jco",
                add = "jitter")+
  facet_grid(cols = vars(Habitat))+
  theme(legend.position = "none")
p2
#  Add p-value
p2 = p2  + stat_compare_means()
p2

#Combine plots

p3 = p1 / p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p3

pdf("../results/alpha-comparisons-observed.pdf", width=10, height=10)
print(p3)
dev.off()


#16S
pro.alpha = estimate_richness(ps.pro2)
pro.alpha$Sample = sample_names(ps.pro2)
pro.alpha$Treatment = data.frame(sample_data(ps.pro2))$Treatment
pro.alpha$Habitat = data.frame(sample_data(ps.pro2))$Habitat
pro.alpha$Habitat = factor(pro.alpha$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

p1 <- ggboxplot(pro.alpha, x = "Treatment", y = "Shannon",
                color = "Treatment", palette = "jco",
                add = "jitter")+
  facet_grid(cols = vars(Habitat))+
  theme(legend.position = "none")
p1
#  Add p-value
p1 = p1  + stat_compare_means()
p1


#18S
euk.alpha = estimate_richness(ps.euk2)
euk.alpha$Sample = sample_names(ps.euk2)
euk.alpha$Treatment = data.frame(sample_data(ps.euk2))$Treatment
euk.alpha$Habitat = data.frame(sample_data(ps.euk2))$Habitat
euk.alpha$Habitat = factor(euk.alpha$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

p2 <- ggboxplot(euk.alpha, x = "Treatment", y = "Shannon",
                color = "Treatment", palette = "jco",
                add = "jitter")+
  facet_grid(cols = vars(Habitat))+
  theme(legend.position = "none")
p2
#  Add p-value
p2 = p2  + stat_compare_means()
p2

#Combine plots

p3 = p1 / p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p3

pdf("../results/alpha-comparisons-shannon.pdf", width=10, height=10)
print(p3)
dev.off()


