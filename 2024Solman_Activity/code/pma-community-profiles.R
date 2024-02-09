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

#modify tax names
name_mod <- function(ps){
  
  ps.m = ps
  x = data.frame(tax_table(ps.m)) %>% 
    mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
    mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
    mutate(Class = ifelse(Class == "Cyanobacteriia", "Cyanophyceae", Class)) %>%
    mutate(Class = ifelse(Class == "uncultured", paste0(Phylum, " uncult."), Class)) %>%
    mutate(Class = ifelse(Class == "Incertae_Sedis", paste0(Phylum, " (incertae sedis)"), Class)) %>%
    mutate(Class = ifelse(Class == "Incertae Sedis", paste0(Phylum, " (incertae sedis)"), Class)) %>%
    mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
    mutate(Order = ifelse(Order == "Leptolyngbyales", "Synechococcales", Order)) %>%
    mutate(Order = ifelse(Order == "uncultured", paste0(Class, " uncult."), Order)) %>%
    mutate(Order = ifelse(Order == "Incertae_Sedis", paste0(Class, " (incertae sedis)"), Order)) %>%
    mutate(Order = ifelse(Order == "Incertae Sedis", paste0(Class, " (incertae sedis)"), Order)) %>%
    mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
    mutate(Family = ifelse(Family == "uncultured", paste0(Order, " uncult."), Family)) %>%
    mutate(Family = ifelse(Family == "Unknown_Family", paste0(Order, " (family unknown)"), Family)) %>%
    mutate(Family = ifelse(Family == "Unknown Family", paste0(Order, " (family unknown)"), Family)) %>%
    mutate(Family = ifelse(Family == "Incertae_Sedis", paste0(Order, " (incertae sedis)"), Family)) %>%
    mutate(Family = ifelse(Family == "Incertae Sedis", paste0(Order, " (incertae sedis)"), Family)) %>%
    mutate(Family = ifelse(Family == "Pleosporales", "Pleosporales (Genus Unknown)", Family)) %>%
    mutate(Family = ifelse(Family == "Rhizophydiales", "Rhizophydiales (Genus Unknown)", Family)) %>%
    mutate(Family = ifelse(Family == "Thecofilosea", "Thecofilosea (Genus Unknown)", Family)) %>%
    mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
    mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus)) %>%
    mutate(Genus = ifelse(Genus == "Unknown_Family", paste0(Order, " (family and genus unknown)"), Genus)) %>%
    mutate(Phylum = ifelse(Class == "Labyrinthulomycetes", "Bigyra", Phylum)) %>%
    mutate(Phylum = ifelse(Phylum == "Labyrinthulomycetes", "Bigyra", Phylum)) %>%
    mutate(Class = ifelse(Class == "Labyrinthulomycetes (Genus Unknown)", "Labyrinthulomycetes", Class)) %>%
    mutate(Class = ifelse(Order == "Burkholderiales", "Betaproteobacteria", Class)) %>%
    mutate(Phylum = ifelse(Phylum == "Actinobacteria", "Actinobacteriota", Phylum)) %>%
    mutate(Class = ifelse(Class == "Actinobacteria", "Actinomycetia", Class)) %>%
    mutate(Class = ifelse(Genus == "Ferruginibacter", "Chitinophagia", Class)) %>%
    mutate(Class = ifelse(Genus == "Hymenobacter", "Cytophagia", Class))%>%
    mutate(Class = ifelse(Genus == "Arcicella", "Cytophagia", Class))%>%
    mutate(Class = ifelse(Order == "Cytophagales", "Cytophagia", Class)) %>%
    mutate(Class = ifelse(Order == "Flavobacteriales", "Flavobacteriia", Class)) %>%
    mutate(Class = ifelse(Order == "Sphingobacteriales", "Sphingobacteriia", Class)) %>%
    mutate(Class = ifelse(Order == "Chitinophagales", "Chitinophagia", Class))%>%
    mutate(Order = ifelse(Genus == "Nostoc_PCC-73102", "Nostocales", Order)) %>%
    mutate(Order = ifelse(Family == "Phormidiaceae", "Oscillatoriales", Order)) %>%
    mutate(Order = ifelse(Family == "Nostocaceae", "Nostocales", Order))
  
  
  
  
  #replace anything that says NA with Genus Unknown
  y = x
  y$Phylum <- gsub("P.NA", "(Phylum Unknown)", y$Phylum)
  y$Class <- gsub("P.NA C.NA", "(Class Unknown)", y$Class)
  y$Class <- gsub("C.NA", "(Class Unknown)", y$Class)
  y$Order <- gsub("P.NA C.NA O.NA", "(Order Unknown)", y$Order)
  y$Order <- gsub("C.NA O.NA", "(Order Unknown)", y$Order)
  y$Order <- gsub("O.NA", "(Order Unknown)", y$Order)
  y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Family Unknown)", y$Family)
  y$Family <- gsub("C.NA O.NA F.NA", "(Family Unknown)", y$Family)
  y$Family <- gsub("O.NA F.NA", "(Family Unknown)", y$Family)
  y$Family <- gsub("F.NA", "(Family Unknown)", y$Family)
  y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
  y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
  y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
  y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
  y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)
  y$Class <- gsub("Labyrinthulomycetes (Genus Unknown)", "Labyrinthulomycetes", y$Class)
  y$Class <- ifelse(y$Class %in% c("Labyrinthulomycetes (Genus Unknown)"), "Labyrinthulomycetes", y$Class)
  y$Class <- ifelse(y$Class %in% c("Labyrinthulomycetes uncult."), "Labyrinthulomycetes", y$Class)
  y$Order <- ifelse(y$Order %in% c("Thecofilosea (Genus Unknown)"), "Thecofilosea", y$Order)
  
  
  tax_table(ps.m) = as.matrix(y)
  
  return(ps.m)
}

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
pro.m = name_mod(pro)

#find the three most abundant Class within each phyla
#get our data
data <-
  pro.m %>%
  tax_glom("Class") %>%
  psmelt() %>%
  as_tibble()

#get 2 most abundant Class by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Class) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 2)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())

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
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Class=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#out phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#remove unwanted genera and replace with other
for (i in 1:nrow(data2plot)){
  if(!data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": Other")
  } else if (data2plot$Class[i] %in% out$Class){
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Class[i])
  }
}


#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
    data2plot$Class[i] <- "Other"
  }
}


#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Class="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Class, data2plot$Class)) == length(unique(data2plot$Class)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")

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
  group_by(Class, Habitat, Treatment, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))



#make axis title
my_y_title <- expression(paste(italic("tDNA"), " Relative Abundance"))

p1a = ggplot(df.plot2[df.plot2$Treatment == "tDNA",], aes(fill=Class, y=Abundance, x=Sample)) + 
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

p1b = ggplot(df.plot2[df.plot2$Treatment == "iDNA",], aes(fill=Class, y=Abundance, x=Sample)) + 
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
euk.m = name_mod(euk)

#find the three most abundant Class within each phyla
#get our data
data <-
  euk.m %>%
  tax_glom("Class") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant Class by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Class) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 2)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())

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
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Class=paste0(unique(out.col$Phylum)[i], ": Other")))
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
length(unique(data2plot$Class))

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Class2 = paste0(data2plot$Phylum, ": ", data2plot$Class)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Class2[i] %in% out.col$Class){
    #if the Class' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Class2[i] %in% out.col$Class){ 
    #else they are just given the phyla name with the Class name
    data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Class[i])
    
  }
}

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "Class" with Other
    data2plot$Class[i] <- "Other"
  }
}

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Class="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Class, data2plot$Class)) == length(unique(data2plot$Class)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")

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
  group_by(Class, Habitat, Treatment, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))

df.plot2$Class <- gsub("Genus", "Class", df.plot2$Class)

#make axis title
my_y_title <- expression(paste(italic("tDNA"), " Relative Abundance"))

p2a = ggplot(df.plot2[df.plot2$Treatment == "tDNA",], aes(fill=Class, y=Abundance, x=Sample)) + 
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

p2b = ggplot(df.plot2[df.plot2$Treatment == "iDNA",], aes(fill=Class, y=Abundance, x=Sample)) + 
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


p5 = (p1a2 + ggtitle("A: Prokaryote") + theme(plot.title = element_text(size=40),
                                              strip.text = element_text(size=20),
                                              legend.text = element_text(size=16),
                                              axis.text = element_text(size=18),
                                              axis.title = element_text(size=18)))  / (p1b2 + theme(plot.title = element_text(size=40),
                                                                                                    strip.text = element_text(size=20),
                                                                                                    legend.text = element_text(size=16),
                                                                                                    axis.text = element_text(size=18),
                                                                                                    axis.title = element_text(size=18)))/ plot_spacer() / (p2a2 + ggtitle("B: Microbial Eukaryote") + theme(plot.title = element_text(size=40),
                                                                                                                                                                                                            strip.text = element_text(size=20),
                                                                                                                                                                                                            legend.text = element_text(size=16),
                                                                                                                                                                                                            axis.text = element_text(size=18),
                                                                                                                                                                                                            axis.title = element_text(size=18))) / (p2b2 + theme(plot.title = element_text(size=40),
                                                                                                                                                                                                                                                               strip.text = element_text(size=20),
                                                                                                                                                                                                                                                               legend.text = element_text(size=16),
                                                                                                                                                                                                                                                               axis.text = element_text(size=18),
                                                                                                                                                                                                                                                               axis.title = element_text(size=18))) + plot_layout(heights = c(1,1,0.1,1,1))

pdf("../results/pma-all-community-plot2.pdf", width=18, height=25)
print(p5)
dev.off()
