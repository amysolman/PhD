
# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
#library(vegan)
source("00-solman-functions.R")
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

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk.rm <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 
#18S WITH MICROFAUNA
ps.full <- readRDS("../results/18S-mm-included-ps-norm.rds")

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#STEP ONE: find the 8 phyla with the highest mean relative abundance for each subcommunity in each pole (ignoring NA phyla)
#find the 8 most abundant Phyla
top.phy <- top_taxa(ps.pro, 
                    tax_level = "Phylum", 
                    n_taxa = 8,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#STEP TWO: modify tax names so it looks better when we plot them
pro.m = ps.pro
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

#STEP THREE: find the three most abundant genus within each phyla by agglomerating data at the 
#genus level then finding the genera with the highest mean abundance in each phylum

#find the three most abundant genus within each phyla
#get our data
data <-
  pro.m %>%
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
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Location <- factor(df$Location, levels = c("Lower", "Upper"))
df$HabLoc <- paste0(df$Habitat, " (", df$Location, ")")
df$HabLoc <- factor(df$HabLoc, levels = c("Snow (Lower)", "Snow (Upper)", "Spring Ice (Lower)", "Spring Ice (Upper)",
                                          "Summer Ice (Lower)", "Summer Ice (Upper)", "Cryoconite (Lower)", "Cryoconite (Upper)"))

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
  group_by(Genus, HabLoc, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p1 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) +
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=5, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

p1

pdf("../results/prokaryote-community-plot.pdf", width=15, height=10)
print(p1)
dev.off()


######PROKARYOTE PLOT FOR PRESENTATION

#for creating space between rows
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

p1a = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1, key_glyph = "polygon3")+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=4, byrow=FALSE))+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"))

p1a

pdf("../results/prokaryote-community-plot-for-presentation.pdf", width=27, height=18)
print(p1a)
dev.off()


######SEND PLOT FOR PRESENTATION
df.plot2a = 
  df.plot %>%
  group_by(Genus, HabLoc) %>%
  dplyr::summarise(across(c(Abundance), sum))


p1b = ggplot(df.plot2a, aes(fill=Genus, y=Abundance, x=HabLoc)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  ylab("Relative Abundance")+ 
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol=5, byrow=FALSE))+
  theme(legend.position="bottom",
        #axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(),
        legend.text = element_text(size=25),
        strip.text.x = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=30),
        axis.title.y = element_text(size=30),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"))

p1b

pdf("../results/prokaryote-community-plot-for-presentation2.pdf", width=40, height=22)
print(p1b)
dev.off()


##########################plot for eukaryotes#########################

#STEP ONE: find the 8 phyla with the highest mean relative abundance for each subcommunity in each pole (ignoring NA phyla)
#find the 8 most abundant Phyla
top.phy <- top_taxa(ps.euk.rm, 
                    tax_level = "Phylum", 
                    n_taxa = 8,
                    grouping = "Habitat")
top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)

#STEP TWO: modify tax names so it looks better when we plot them
euk.m = ps.euk.rm
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

#STEP THREE: find the three most abundant genus within each phyla by agglomerating data at the 
#genus level then finding the genera with the highest mean abundance in each phylum

#find the three most abundant genus within each phyla
#get our data
data <-
  euk.m %>%
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
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Location <- factor(df$Location, levels = c("Lower", "Upper"))
df$HabLoc <- paste0(df$Habitat, " (", df$Location, ")")
df$HabLoc <- factor(df$HabLoc, levels = c("Snow (Lower)", "Snow (Upper)", "Spring Ice (Lower)", "Spring Ice (Upper)",
                                          "Summer Ice (Lower)", "Summer Ice (Upper)", "Cryoconite (Lower)", "Cryoconite (Upper)"))

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
  group_by(Genus, HabLoc, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))


p2 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
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


####MICROMETAZOAN PLOT

my_subset <- subset(otu_table(ps.full), rownames(otu_table(ps.full)) %in% rownames(otu_table(ps.mm)))
new_physeq <- merge_phyloseq(my_subset, tax_table(ps.full), sample_data(ps.full))


get_my_data <- function(phylo, habitat){
  
  #Separate habitat data
  ps <- prune_samples(habitat, phylo)
  
  #agglomerate taxa at the species level
  ps2 <- tax_glom(ps, taxrank = "Genus")
  
  #get presence, absence dataframe
  df <- data.frame(otu_table(ps2))
  df[df > 0] <- 1
  
  #add genus names
  rownames(df) == rownames(data.frame(tax_table(ps2))) #should be TRUE
  rownames(df) <- paste0(data.frame(tax_table(ps2))$Genus, " (", data.frame(tax_table(ps2))$Phylum, ")")
  
  #put into long format for plotting
  df2 <- cbind(rownames(df), df)
  names(df2) <- c("Genus", names(df2)[2:length(names(df2))])
  df3 <- gather(df2, Sample, Presence, names(df2)[2]:names(df2)[length(names(df2))], factor_key=TRUE)
  
  #add Habitat/Location variable
  data = data.frame(sample_data(ps2))
  data2add = data.frame(Sample=rownames(data), HabLoc=paste0(data$Habitat, " (", data$Location, ")"))
  df4 = merge(df3, data2add[, c("Sample", "HabLoc")], by="Sample")
  
  return(df4)
  
}

snow.df <- get_my_data(new_physeq, habitat = data.frame(sample_data(subset_samples(new_physeq, Habitat == "Snow")))$SampleID)
sp.ice.df <- get_my_data(new_physeq, habitat = data.frame(sample_data(subset_samples(new_physeq, Habitat == "Spring Ice")))$SampleID)
sum.ice.df <- get_my_data(new_physeq, habitat = data.frame(sample_data(subset_samples(new_physeq, Habitat == "Summer Ice")))$SampleID)
cry.df <- get_my_data(new_physeq, habitat = data.frame(sample_data(subset_samples(new_physeq, Habitat == "Cryoconite")))$SampleID)

#plot
snow.p <- ggplot(snow.df, aes(x = Sample, y = Genus, fill = Presence)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))
#theme(legend.position = "none")

sp.p <- ggplot(sp.ice.df, aes(x = Sample, y = Genus, fill = Presence)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))

sum.p <- ggplot(sum.ice.df, aes(x = Sample, y = Genus, fill = Presence)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))

cry.p <- ggplot(cry.df, aes(x = Sample, y = Genus, fill = Presence)) +
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))

plot_grid(snow.p + theme(legend.position = "none"), 
          sp.p + theme(axis.text.y=element_blank(),  #remove y axis labels
                       axis.ticks.y=element_blank(), axis.title.y=element_blank()), 
          sum.p + theme(legend.position = "none"),
          cry.p + theme(axis.text.y=element_blank(),  #remove y axis labels
                        axis.ticks.y=element_blank(), axis.title.y=element_blank()),
          labels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
#hjust = 0.1, vjust = -0.5,
#scale = 0.9)
#plot

snow.df$Genus <- factor(snow.df$Genus,
                        levels = c("Dorylaimida (Nematozoa)", "Monhysterida (Nematozoa)", "Echiniscoidea (Tardigrada)",  "Parachela (Tardigrada)", "Ploimida (Rotifera)", "Adinetida (Rotifera)"))

snow.p <-ggplot(snow.df, aes(x=Sample, y=Genus, fill = factor(Presence))) +
  geom_point(size = 20, shape = 22)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=30),
        axis.text.y = element_text(size=30),
        legend.title = element_blank(),
        plot.margin = unit(c(1,0.2,0.5,0), "cm"), strip.text = element_text(size=40))+
  scale_fill_discrete(labels=c('Absent', 'Present'))+
  scale_y_discrete(position = "right")
snow.p

sp.ice.df$Genus <- factor(sp.ice.df$Genus,
                          levels = c("Dorylaimida (Nematozoa)", "Monhysterida (Nematozoa)", "Echiniscoidea (Tardigrada)",  "Parachela (Tardigrada)", "Ploimida (Rotifera)", "Adinetida (Rotifera)"))

sp.p <-ggplot(sp.ice.df, aes(x=Sample, y=Genus, fill = factor(Presence))) +
  geom_point(size = 20, shape = 22)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=30),
        axis.text.y = element_text(size=30),
        legend.title = element_blank(),
        plot.margin = unit(c(1,0.2,0.5,0), "cm"), strip.text = element_text(size=40))+
  scale_fill_discrete(labels=c('Absent', 'Present'))+
  scale_y_discrete(position = "right")

sum.ice.df$Genus <- factor(sum.ice.df$Genus,
                           levels = c("Dorylaimida (Nematozoa)", "Monhysterida (Nematozoa)", "Echiniscoidea (Tardigrada)",  "Parachela (Tardigrada)", "Ploimida (Rotifera)", "Adinetida (Rotifera)"))

sum.p <-ggplot(sum.ice.df, aes(x=Sample, y=Genus, fill = factor(Presence))) +
  geom_point(size = 20, shape = 22)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=30),
        axis.text.y = element_text(size=30),
        legend.title = element_blank(),
        plot.margin = unit(c(1,0.2,0.5,0), "cm"), strip.text = element_text(size=40))+
  scale_fill_discrete(labels=c('Absent', 'Present'))+
  scale_y_discrete(position = "right")

cry.df$Genus <- factor(cry.df$Genus,
                       levels = c("Dorylaimida (Nematozoa)", "Monhysterida (Nematozoa)", "Echiniscoidea (Tardigrada)",  "Parachela (Tardigrada)", "Ploimida (Rotifera)", "Adinetida (Rotifera)"))

cry.p <-ggplot(cry.df, aes(x=Sample, y=Genus, fill = factor(Presence))) +
  geom_point(size = 20, shape = 22)+
  facet_grid(cols = vars(HabLoc), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size=30),
        axis.text.y = element_text(size=50),
        legend.title = element_blank(),
        plot.margin = unit(c(1,0.2,0.5,0), "cm"), strip.text = element_text(size=40))+
  scale_fill_discrete(labels=c('Absent', 'Present'))+
  scale_y_discrete(position = "right")

multi.p = plot_grid(snow.p + theme(legend.position = "none", #remove y axis labels
                                   axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title=element_blank()), 
                    sp.p + theme(axis.text.y=element_blank(),  #remove y axis labels
                                 axis.ticks.y=element_blank(), axis.title=element_blank(), legend.position = "none"), 
                    sum.p + theme(axis.text.y=element_blank(),  #remove y axis labels
                                  axis.ticks.y=element_blank(), axis.title=element_blank(), legend.position = "none"),
                    cry.p + theme(legend.position = "none", axis.title.x = element_blank()),
                    #labels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"),
                    ncol=4,
                    vjust = 1, hjust = c(-3, 0, 0, 0), rel_widths = c(.9, .9, .9, 1.4))
#hjust = 0.1, vjust = -0.5,
#scale = 0.9)

leg = get_legend(
  snow.p + 
    theme(legend.position = "bottom",
          legend.text = element_text(size=100), 
          legend.title = element_blank())+
    guides(colour = guide_legend(override.aes = list(size = 150, alpha=.8)))
)

# final.p = plot_grid(multi.p, leg, ncol = 1, rel_heights = c(1, .1))
# final.p

# pdf("../results/mm-presence-absence-graph.pdf", width = 20, height=3.5)
# print(final.p)
# dev.off()

######MERGE ALL PLOTS

title1 <- ggdraw() + draw_label("Prokaryote Community Profiles", size=90, hjust=1.72, fontface = "bold")
title2 <- ggdraw() + draw_label("Eukaryote Community Profiles", size=90, hjust=1.75, fontface = "bold")
title3 <- ggdraw() + draw_label("Microfauna Community Profiles", size=90, hjust=1.7, fontface = "bold")

final.all.coms = title1 / plot_spacer() / (p1 + theme(legend.text = element_text(size=60), strip.text = element_text(size=50), axis.text = element_text(size=30), axis.title = element_text(size=40)) + guides(fill=guide_legend(ncol = 3, byrow=FALSE))) / 
  plot_spacer() / title2 / plot_spacer() / (p2 + theme(legend.text = element_text(size=60), strip.text = element_text(size=50), axis.text = element_text(size=30), axis.title = element_text(size=40)) + guides(fill=guide_legend(ncol = 3, byrow=FALSE))) / 
  plot_spacer() / title3 / plot_spacer() / (multi.p  + theme(legend.text = element_text(size=60), axis.text = element_text(size=40), strip.text = element_text(size=50), axis.title = element_text(size=40))) / plot_spacer()  / leg + plot_layout(heights = c(0.1, 0.05, 1.4, 0.1, 0.1, 0.05, 1.4, 0.1, 0.1, 0.05, 0.8, 0.1, 0.2))
final.all.coms

pdf("../results/full-community-profiles.pdf", height=85, width=65)
print(final.all.coms)
dev.off()