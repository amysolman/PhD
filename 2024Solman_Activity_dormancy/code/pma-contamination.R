#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(tidyverse)
library(phyloseq)
library(readr)
library(seqinr)
library(vegan)
library(magrittr)
library(data.table)
library(ape) #appe::Ntip
library(dplyr)
library(fantaxtic)
library(patchwork)

#2. Import data
pro <- readRDS("../results/pma-16S-phylo-object.rds")
euk <- readRDS("../results/pma-18S-phylo-object.rds")

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#how many reads in control samples?
pro.bar = plot_bar(pro) + theme(axis.text = element_text(size=12),
                                plot.title = element_text(size=40),
                                axis.title = element_text(size=18)) + ggtitle("A: Prokaryote")
euk.bar = plot_bar(euk) + theme(axis.text = element_text(size=12),
                                plot.title = element_text(size=40),
                                axis.title = element_text(size=18)) + ggtitle("B: Microbial Eukaryote")


####################################################################################################################
#mean reads per control sample and environment sample
#output results to text document
sink("../results/pma-control-env-reads.txt", type="output")
writeLines("===============================================================
NUMBER OF READS IN CONTROL SAMPLES VS ENVIRONMENTAL SAMPLES
===============================================================")
writeLines("Mean number of reads in prokaryote control samples:")
mean(colSums(data.frame(otu_table(subset_samples(pro, Habitat == "Control")))))
writeLines("Mean number of reads in prokaryote environmental samples:")
mean(colSums(data.frame(otu_table(subset_samples(pro, Habitat != "Control")))))
writeLines("Mean number of reads in eukaryote control samples:")
mean(colSums(data.frame(otu_table(subset_samples(euk, Habitat == "Control")))))
writeLines("Mean number of reads in eukaryote environmental samples:")
mean(colSums(data.frame(otu_table(subset_samples(euk, Habitat != "Control")))))
sink()
###############################################################################################################

#taxonomy of ASVs in control samples?

#Prokaryotes

#subset to prokaryote control samples
pro.cont = subset_samples(pro, Habitat == "Control")

#STEP ONE: find the 10 most abundant phyla
top.phy <- top_taxa(pro.cont, 
                    tax_level = "Phylum", 
                    n_taxa = 10)

top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
pro.m = pro.cont
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
  guides(fill=guide_legend(ncol = 3, byrow=FALSE))+
  theme(plot.title = element_text(size=40),
        strip.text = element_text(size=25),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        #legend.key.size = unit(3,"line"),
        axis.text.y = element_text(size=18),
        axis.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size=18),
        legend.position="bottom")

p1

pdf("../results/pma-prokaryote-control-plot.pdf", width=35, height=10)
print(p1)
dev.off()


#Eukaryotes

#subset to prokaryote control samples
euk.cont = subset_samples(euk, Habitat == "Control")

x = data.frame(tax_table(euk.cont))

#STEP ONE: find the 10 most abundant phyla
top.phy <- top_taxa(euk.cont, 
                    tax_level = "Phylum", 
                    n_taxa = 10)

top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
euk.m = euk.cont

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
  guides(fill=guide_legend(ncol = 3, byrow=FALSE))+
  theme(plot.title = element_text(size=40),
          strip.text = element_text(size=25),
          legend.text = element_text(size=15),
          legend.title = element_blank(),
          #legend.key.size = unit(3,"line"),
          axis.text.y = element_text(size=18),
          axis.title = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size=18),
        legend.position="bottom")

p2

# pdf("../results/pma-eukaryote-control-plot.pdf", width=35, height=10)
# print(p2)
# dev.off()

#combine our plots
# full.p = (pro.bar + euk.bar) / (p1 + p2)

full.p = pro.bar / plot_spacer() / p1 / plot_spacer() / euk.bar / plot_spacer() /p2 + plot_layout(heights = c(1,0.1,1,0.2,1,0.1,1))

pdf("../results/pma-control-plots.pdf", width=18, height=25)
print(full.p)
dev.off()


#remove contaminants and control samples
#no contaminants were removed due to the low number of reads in control samples

#remove control samples and same new phyloseq objects
pro.no.cont = subset_samples(pro, Habitat != "Control")
euk.no.cont = subset_samples(euk, Habitat != "Control")

saveRDS(pro.no.cont, "../results/pma-16S-phylo-object-no-controls.rds")
saveRDS(euk.no.cont, "../results/pma-18S-phylo-object-no-controls.rds")
