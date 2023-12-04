
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
#devtools::install_github("gmteunisse/fantaxtic")
library(fantaxtic)
library(tidyverse)
library(patchwork)
library(ggpubr)

#load data
ps.pro <- readRDS("../data/16S-phylo-object.rds") 
ps.euk <- readRDS("../data/18S-phylo-object.rds") 

#remove the 2 most abundant taxa from ps.euk
#which ASVs have the most counts?
tax = data.frame(tax_table(ps.euk))
x = data.frame(sort(taxa_sums(ps.euk), decreasing = FALSE))
names(x) = "value"
y = x %>% 
  top_n(2, value)
z = tax[rownames(tax) %in% rownames(y),]
allTaxa = taxa_names(ps.euk)
allTaxa <- allTaxa[!(allTaxa %in% rownames(z))]
ps.euk.rm = prune_taxa(allTaxa, ps.euk)

#RELATIVE ABUNDANCE BAR CHART FUNCTIONS

#data for testing
# ps = ps.euk
# controls = TRUE
# data_name = "18S"
# w = 20
# h = 10

# ps = ps.pro
# controls = FALSE
# data_name = "16S"
# w = 15
# h = 10

bar_chart_func <- function(ps, controls, data_name, w, h){
  
  
  # if (controls == TRUE){
  #   #remove control samples
  #   ps = subset_samples(ps, SampleType != "Control")
  # }

  #remove taxa with zero countrs
  ps2 = filter_taxa(ps, function(x) sum(x) > 0, TRUE)
  
  #modify tax names
  x = data.frame(tax_table(ps2)) %>% 
    mutate(Phylum = ifelse(is.na(Phylum), paste0(Kingdom, " P.",Phylum), Phylum)) %>%
    mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
    mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
    mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
    mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
    mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))
  tax_table(ps2) = as.matrix(x)
  
  #find the 10 most abundant Phyla
  top.phy <- top_taxa(ps2, 
                      tax_level = "Phylum", 
                      n_taxa = 10,
                      grouping = "SampleType")
  top.phy2 = top.phy[[2]]

  phy.keep = unique(top.phy2$Phylum)
  
  #find the three most abundant genus within each phyla
  #get our data
  data <-
    ps2 %>%
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
  
  
  #sort the out dataframe and add colours
  out.col = out[with(out, order(Phylum)), ]
  out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
  num4cols = 1
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
  
  
  #remove unwanted genera and replace with other
  for (i in 1:nrow(data2plot)){
    if(!data2plot$Genus[i] %in% out$Genus){
      data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    } else if (data2plot$Genus[i] %in% out$Genus){
      data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    }
  }
  
  #remove unwanted phyla
  for (i in 1:nrow(data2plot)){
    if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
      data2plot$Genus[i] <- "Other"
    }
  }
  
  #add those colours to the corresponding row in our main dataframe
  #add "Other" to out.col
  out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
  
  data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")
  
  #plot data
  df = data2plot2
  if (controls == TRUE){
  df$SampleType <- factor(df$SampleType, levels = c("Interface", "Meltwater", "Surface", "Control"))
  df$Site <- factor(df$Site, levels = c("A", "D", "G"))
  } else {
    df$SampleType <- factor(df$SampleType, levels = c("Interface", "Meltwater", "Surface"))
    df$Site <- factor(df$Site, levels = c("A", "D"))
  }
  df$Day <- factor(df$Day, levels = c("1", "2", "3"))
  
  
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
    group_by(Genus, SampleType, Sample) %>%
    dplyr::summarise(across(c(Abundance), sum))
  
  
  p1 = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
    geom_bar(position="fill", stat="identity", alpha=1)+
    facet_grid(cols = vars(SampleType), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
    ylab("Relative Abundance")+ 
    theme_bw()+
    scale_fill_manual(values = col)+
    guides(fill=guide_legend(ncol=5, byrow=FALSE))+
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5, size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=15),
          axis.title.y = element_text(size=15),
          legend.title=element_blank(),
          strip.text = element_text(
            size = 20))
  
  p1
  
  pdf(paste0("../results/", data_name, "-community-plot.pdf"), width=w, height=h)
  print(p1)
  dev.off()
  
  #save colours
  #write.csv(gen.df, paste0("../results/", data, "-colours-for-plotting-rare.csv"))
  
  return(p1)
  
}

p1 = bar_chart_func(ps = ps.pro, controls = FALSE, data_name = "16S", w = 15, h = 10)
p2 = bar_chart_func(ps = ps.euk, controls = TRUE, data_name = "18S", w = 20, h = 10)
p2a = bar_chart_func(ps = ps.euk.rm, controls = TRUE, data_name = "18S-rm", w = 20, h = 10)

###PUT THE TWO PLOTS TOGETHER

p3 = p1 / p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p3

pdf("../results/all-community-plot.pdf", width=18, height=20)
print(p3)
dev.off()

#and after removing those most abundant ASVs from the euk data
p4 = p1 / p2a + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p4

pdf("../results/all-community-plot-euk-rm.pdf", width=18, height=20)
print(p4)
dev.off()




#WHATS IN OUR CONTROLS?

euk.c = subset_samples(ps.euk, SampleType == "Control")
euk.c = filter_taxa(euk.c, function(x) sum(x) > 0, TRUE)
euk.c


#Community structure numbers

#How many ASVs
ps.pro

#number of ASVs after removing blanks
ps.euk.no.blank = subset_samples(ps.euk, SampleType != "Control")
ps.euk.no.blank = filter_taxa(ps.euk.no.blank, function(x) sum(x) > 0, TRUE)
ps.euk.no.blank

#16S interface
sub = subset_samples(ps.pro, SampleType == "Interface")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S meltwater
sub = subset_samples(ps.pro, SampleType == "Meltwater")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#16S surface
sub = subset_samples(ps.pro, SampleType == "Surface")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)


#18S interface
sub = subset_samples(ps.euk.no.blank, SampleType == "Interface")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S meltwater
sub = subset_samples(ps.euk.no.blank, SampleType == "Meltwater")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#18S surface
sub = subset_samples(ps.euk.no.blank, SampleType == "Surface")
glom <- tax_glom(sub, taxrank = 'Phylum')
t.glom = data.frame(otu_table(glom))
rownames(t.glom) = data.frame(tax_table(glom))$Phylum
sort(rowSums(t.glom))/sum(t.glom)

#exploring Camptobasidium
taxa.camp = data.frame(tax_table(ps.euk.no.blank))
taxa.camp.sub = taxa.camp[which(taxa.camp$Genus == "Camptobasidium"),]
euk.camp = prune_taxa(rownames(taxa.camp.sub), ps.euk)
euk.camp
count.camp = data.frame(otu_table(euk.camp))
sort(rowSums(count.camp))


#any archeae and if so which communities are they in? = subset_taxa(p.ar.a, Kingdom == "Archaea")

#NO ARCHAEA
# arch.sub = subset_taxa(ps.pro, Kingdom == "Archaea")
# #remove samples with zero counts
# arch.sub1 = prune_samples(sample_sums(arch.sub) >0, arch.sub)
# #remove taxa with zero counts
# arch.sub2 = filter_taxa(arch.sub1, function(x) sum(x) > 0, TRUE)
# 
# x = data.frame(tax_table(arch.sub2))
# sort(table(x$Family), decreasing = TRUE)
# arch.counts = data.frame(otu_table(arch.sub2))
# sort(colSums(arch.counts), decreasing = TRUE)


####Alpha Diversity - NOTE THESE ARE ONLY VALID IF USING RAREFIED DATA

#16S
pro.alpha = estimate_richness(ps.pro)
pro.alpha$Sample = sample_names(ps.pro)
pro.alpha$SampleType = data.frame(sample_data(ps.pro))$SampleType
pro.alpha$SampleType = factor(pro.alpha$SampleType, levels = c("Interface", "Meltwater", "Surface"))

p1 <- ggboxplot(pro.alpha, x = "SampleType", y = "Observed",
                color = "SampleType", palette = "jco",
                add = "jitter")+
  theme(legend.position = "none")
p1

#  Add p-value
p1 = p1  + stat_compare_means()
p1


#18S
euk.alpha = estimate_richness(ps.euk.no.blank)
euk.alpha$Sample = sample_names(ps.euk.no.blank)
euk.alpha$SampleType = data.frame(sample_data(ps.euk.no.blank))$SampleType
euk.alpha$SampleType = factor(euk.alpha$SampleType, levels = c("Interface", "Meltwater", "Surface"))

p2 <- ggboxplot(euk.alpha, x = "SampleType", y = "Observed",
                color = "SampleType", palette = "jco",
                add = "jitter")+
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
pro.alpha = estimate_richness(ps.pro)
pro.alpha$Sample = sample_names(ps.pro)
pro.alpha$SampleType = data.frame(sample_data(ps.pro))$SampleType
pro.alpha$SampleType = factor(pro.alpha$SampleType, levels = c("Interface", "Meltwater", "Surface"))

p1 <- ggboxplot(pro.alpha, x = "SampleType", y = "Shannon",
                color = "SampleType", palette = "jco",
                add = "jitter")+
  theme(legend.position = "none")
p1

#  Add p-value
p1 = p1  + stat_compare_means()
p1


#18S
euk.alpha = estimate_richness(ps.euk.no.blank)
euk.alpha$Sample = sample_names(ps.euk.no.blank)
euk.alpha$SampleType = data.frame(sample_data(ps.euk.no.blank))$SampleType
euk.alpha$SampleType = factor(euk.alpha$SampleType, levels = c("Interface", "Meltwater", "Surface"))

p2 <- ggboxplot(euk.alpha, x = "SampleType", y = "Shannon",
                color = "SampleType", palette = "jco",
                add = "jitter")+
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


