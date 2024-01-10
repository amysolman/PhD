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
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
#eukaryotes without micrometazoans
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#####16S RELATIVE ABUNDANCE BAR CHART##########

#STEP ONE: find the 10 most abundant phyla in each habitat
top.phy <- top_taxa(pro, 
                    tax_level = "Phylum", 
                    n_taxa = 8,
                    grouping = c("Habitat"))

top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
pro.m = pro
x = data.frame(tax_table(pro.m)) %>% 
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

tax_table(pro.m) = as.matrix(y)

#find the three most abundant genus within each phyla
#get our data
data <-
  pro.m %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = pro.cols[pro.cols$Phylum == mini.df$Phylum,]$Colour1
  col2 = pro.cols[pro.cols$Phylum == mini.df$Phylum,]$Colour2
  
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

#sanity check
#number of all phyla in our dataframe
length(unique(data2plot$Phylum))

#number of all genera in our dataframe
length(unique(data2plot$Genus))

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

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Genus, data2plot$Genus)) == length(unique(data2plot$Genus)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Treatment <- factor(df$Treatment, levels = c("tDNA", "iDNA"))
df$HabTreat <- paste0(df$Habitat, " (", df$Treatment, ")")
df$HabTreat <- factor(df$HabTreat, levels = c("Snow (tDNA)", "Snow (iDNA)", "Spring Ice (tDNA)", 
                                              "Spring Ice (iDNA)", "Summer Ice (tDNA)", "Summer Ice (iDNA)",
                                              "Cryoconite (tDNA)", "Cryoconite (iDNA)"))

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
  group_by(Genus, Habitat, Treatment, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))

#make axis title
my_y_title <- expression(paste(italic("tDNA"), " Relative Abundance"))

p1a = ggplot(df.plot2[df.plot2$Treatment == "tDNA",], aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y=my_y_title)+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=3, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p1a

#change axis text
df.plot2$Sample <- gsub('.{1}$', '', df.plot2$Sample)

#make axis title
my_y_title <- expression(paste(italic("iDNA"), " Relative Abundance"))

p1b = ggplot(df.plot2[df.plot2$Treatment == "iDNA",], aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y=my_y_title)+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=3, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p1b

#combine plots
p1a2 = p1a+theme(legend.position = "none", 
                 axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())
p1b2 = p1b+theme(strip.text.x = element_blank())
pro.p = p1a2 / p1b2

pdf("../results/pma-prokaryote-community-plot.pdf", width=15, height=10)
print(pro.p)
dev.off()


#####18S RELATIVE ABUNDANCE BAR CHART##########


#STEP ONE: find the 10 most abundant phyla in each habitat
top.phy <- top_taxa(euk, 
                    tax_level = "Phylum", 
                    n_taxa = 8,
                    grouping = c("Habitat"))

top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
euk.m = euk
x = data.frame(tax_table(euk.m)) %>% 
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

tax_table(euk.m) = as.matrix(y)

#find the three most abundant genus within each phyla
#get our data
data <-
  euk.m %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 3)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = euk.cols[euk.cols$Phylum == mini.df$Phylum,]$Colour1
  col2 = euk.cols[euk.cols$Phylum == mini.df$Phylum,]$Colour2
  
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

#sanity check
#number of all phyla in our dataframe
length(unique(data2plot$Phylum))

#number of all genera in our dataframe
length(unique(data2plot$Genus))

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

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Genus, data2plot$Genus)) == length(unique(data2plot$Genus)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Treatment <- factor(df$Treatment, levels = c("tDNA", "iDNA"))
df$HabTreat <- paste0(df$Habitat, " (", df$Treatment, ")")
df$HabTreat <- factor(df$HabTreat, levels = c("Snow (tDNA)", "Snow (iDNA)", "Spring Ice (tDNA)", 
                                              "Spring Ice (iDNA)", "Summer Ice (tDNA)", "Summer Ice (iDNA)",
                                              "Cryoconite (tDNA)", "Cryoconite (iDNA)"))

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
  group_by(Genus, Habitat, Treatment, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))

#make axis title
my_y_title <- expression(paste(italic("tDNA"), " Relative Abundance"))

p2a = ggplot(df.plot2[df.plot2$Treatment == "tDNA",], aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y=my_y_title)+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=3, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p2a

#change axis text
df.plot2$Sample <- gsub('.{1}$', '', df.plot2$Sample)

#make axis title
my_y_title <- expression(paste(italic("iDNA"), " Relative Abundance"))

p2b = ggplot(df.plot2[df.plot2$Treatment == "iDNA",], aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  labs(y=my_y_title)+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=3, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p2b

#combine plots
p2a2 = p2a+theme(legend.position = "none", 
                 axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.ticks.x = element_blank())

p2b2 = p2b+theme(strip.text.x = element_blank())
euk.p = p2a2 / p2b2
euk.p

pdf("../results/pma-eukaryote-community-plot.pdf", width=15, height=10)
print(euk.p)
dev.off()

p4 = p1a2 / p1b2 / p2a2 / p2b2 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 25))
p4

pdf("../results/pma-all-community-plot.pdf", width=18, height=20)
print(p4)
dev.off()


p5 = (p1a2 + ggtitle("Prokaryote") + theme(plot.title = element_text(size=40),
                                           strip.text.x = element_text(size=20))) / (p1b2 + theme(legend.text = element_text(size=15)))/ (p2a2 + ggtitle("Eukaryote") + theme(plot.title = element_text(size=40),
                                                                                                                                strip.text.x = element_text(size=20))) / (p2b2 + theme(legend.text = element_text(size=15)))

pdf("../results/pma-all-community-plot2.pdf", width=18, height=25)
print(p5)
dev.off()
