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
pro.DNA <- readRDS("../results/mt-dna-16S-phylo-object-no-controls.rds") 
pro.RNA <- readRDS("../results/mt-rna-16S-phylo-object-no-controls.rds") 

#eukaryotes without micrometazoans
euk.DNA <- readRDS("../results/mt-dna-18S-phylo-object-no-controls.rds") 
euk.RNA <- readRDS("../results/mt-rna-18S-phylo-object-no-controls.rds") 

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#combine the DNA and RNA phyloseq objects into one

comb_phylo <- function(ps1, ps2){
  
  #get otu tables
  tab1 = data.frame(otu_table(ps1))
  tab2 = data.frame(otu_table(ps2))
  
  #rename the samples
  names(tab1) = paste0(names(tab1), ".DNA")
  names(tab2) = paste0(names(tab2), ".RNA")
  
  #combine the otu tables
  total.tab = rbind.fill(tab1, tab2)
  rownames(total.tab) = c(rownames(tab1), rownames(tab2)) #give the ASVs/NTUs the right names
  total.tab[is.na(total.tab)] <- 0 #make NA into zero
  
  #get taxa tables
  tax.tab1 = data.frame(tax_table(ps1))
  tax.tab2 = data.frame(tax_table(ps2))
  
  #combine the taxa tables (excluding the columns they don't share)
  total.tax = rbind(tax.tab1[, !(names(tax.tab1) %in% 'Species')], tax.tab2)
  
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
  ps = phyloseq(otu_table(total.tab, taxa_are_rows = TRUE), tax_table(total.tax), sample_data(total.meta))
  

  return(ps)
}

#run the function
pro.total <- comb_phylo(pro.DNA, pro.RNA)
euk.total <- comb_phylo(euk.DNA, euk.RNA)

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

#run the function
pro <- name_mod(pro.total)
euk <- name_mod(euk.total)

#plot function 

# ps = euk
# col.df = euk.cols
# gene = "18S"
# num_groups = 2
# n_col = 3

com_plot_function <- function(ps, col.df, gene, num_groups, n_col){
  
  #find the 10 most abundant Phyla
  top.phy <- top_taxa(ps, 
                      tax_level = "Phylum", 
                      n_taxa = 10,
                      grouping = "Habitat")
  top.phy2 = top.phy[[2]]
  length(unique(top.phy2$Phylum))
  phy.keep = unique(top.phy2$Phylum)
  
  #find the three most abundant Class within each phyla
  #get our data
  data <-
    ps %>%
    tax_glom("Class") %>%
    psmelt() %>%
    as_tibble()
  
  #get 3 most abundant Class by phylum
  out = data %>%
    filter(Phylum %in% phy.keep) %>%
    group_by(Phylum, Class) %>%
    dplyr::summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance)%>%
    top_n(n = num_groups)
  
  
  #sort the out dataframe and add colours
  out.col = out[with(out, order(Phylum)), ]
  out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
  save.cols <- vector()
  df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())
  
  for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
    
    #get our phyla data only
    mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
    
    #which colours are assigned to that phyla?
    col1 = col.df[col.df$Phylum == mini.df$Phylum,]$Colour1
    col2 = col.df[col.df$Phylum == mini.df$Phylum,]$Colour2
    
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
  #length(intersect(out.col2$Class, data2plot$Class)) == length(unique(data2plot$Class)) #should = TRUE
  data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")
  
  
  #replace Genus with Class
  #data2plot2$Class <- gsub('Genus', 'Class', data2plot2$Class)
  
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

  
  p1 = ggplot(df.plot2, aes(fill=Class, y=Abundance, x=SampleID)) + 
    geom_bar(position="fill", stat="identity", alpha=1)+
    facet_grid(cols = vars(Habitat), rows = vars(Data), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
    labs(y=paste0("Relative abundance of ", gene, " rRNA (RNA) and rRNA genes (DNA)"), x = "Sample")+
    theme_bw()+
    scale_fill_manual(values = col)+
    guides(fill=guide_legend(ncol=n_col, byrow=FALSE))+
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          legend.title=element_blank())
  print(p1)
  
  
  p2 = ggplot(df.plot2, aes(fill=Class, y=Abundance, x=Habitat)) + 
    geom_bar(position="fill", stat="identity", alpha=1)+
    facet_grid(cols = vars(Habitat), rows = vars(Data), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
    labs(y=paste0("Relative abundance of ", gene, " rRNA (RNA) and rRNA genes (DNA)"), x = "Sample")+
    theme_bw()+
    scale_fill_manual(values = col)+
    guides(fill=guide_legend(ncol=n_col, byrow=FALSE))+
    theme(legend.position="bottom",
          axis.text.x = element_blank(),
          legend.title=element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  print(p2)
  
  
  return(list(p1, p2))

  
}


#run function
pro.p = com_plot_function(pro, pro.cols, "16S", 2, 3)
euk.p = com_plot_function(euk, euk.cols, "18S", 2, 3)
euk.p[[1]]
euk.p[[2]]

#merge 16S and 18S results into one plot

final.p = (pro.p[[1]] +ggtitle("A: Prokaryote") + theme(plot.title = element_text(size=40),
                                                     strip.text = element_text(size=20),
                                                     legend.text = element_text(size=16),
                                                     axis.text = element_text(size=18),
                                                                  axis.title = element_text(size=18))) / plot_spacer() / (euk.p[[1]] + ggtitle("B: Microbial Eukaryote")+ theme(plot.title = element_text(size=40),
                                                                                                                                                  strip.text = element_text(size=20),
                                                                                                                                                  legend.text = element_text(size=16),
                                                                                                                                                  axis.text = element_text(size=18),
                                                                                                                                                  axis.title = element_text(size=18))) + plot_layout(heights = c(1,0.1,1))
final.p


pdf("../results/my-all-community-plot.pdf", width=18, height=25)
print(final.p)
dev.off()


