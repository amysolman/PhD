#1. clear workspace and load packages
rm(list=ls())
library(ggplot2)
library(phyloseq)
library(patchwork)
library(ggpubr)
library(glue)
library(tibble) # for add_column
library(reshape2)
library(cowplot)
library(betapart)
library(dplyr)
require("fantaxtic")
library(gtools)

# 2. Import data
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 

#1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.
#2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
#3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
#4. Plot total relative abundances from untreated samples and viability ratios.

#Step One: #1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.

#phylo = pro

viability_phylum <- function(phylo){
  
  #subset to tDNA samples only
  noPMA = subset_samples(phylo, Treatment == "tDNA")
  #remove taxa with zero countrs
  noPMA = filter_taxa(noPMA, function(x) sum(x) > 0, TRUE)
  #remove samples with zero counts
  noPMA = prune_samples(sample_sums(noPMA)>=1, noPMA)
  
  #get our data
  data <- noPMA %>%
    tax_glom("Phylum") %>%
    psmelt() %>%
    as_tibble()
  
  #group data by SampleID and ASV and 
  out = data %>%
    group_by(SampleID, Phylum) %>% #for each ASV in each sample
    dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance
  
  #from this get the relative abundance of each 
  out2 = out %>%
    group_by(SampleID) %>%
    dplyr::summarise(NoPMARelAbun = Abundance/sum(Abundance))
  
  #give ASV IDs
  out2$Phylum = out$Phylum
  
  #sanity check
  test = out2[out2$SampleID == "S21.25",]
  sum(test$NoPMARelAbun) #should equal 1
  
  #Step Two: #2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
  
  #subset to PMA samples only
  PMA = subset_samples(phylo, Treatment == "iDNA")
  #remove taxa with zero countrs
  PMA = filter_taxa(PMA, function(x) sum(x) > 0, TRUE)
  #remove samples with zero counts
  PMA = prune_samples(sample_sums(PMA)>=1, PMA)
  
  #get our data
  data <- PMA %>%
    tax_glom("Phylum") %>%
    psmelt() %>%
    as_tibble()
  
  #group data by SampleID and ASV and 
  out3 = data %>%
    group_by(SampleID, Phylum) %>% #for each ASV in each sample
    dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance
  
  #from this get the relative abundance of each 
  out4 = out3 %>%
    group_by(SampleID) %>%
    dplyr::summarise(PMARelAbun = Abundance/sum(Abundance))
  
  #give ASV IDs
  out4$Phylum = out3$Phylum
  
  #sanity check
  test = out4[out4$SampleID == "S21.25",]
  sum(test$PMARelAbun) #should equal 1
  
  
  #Step Three: #3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
  
  #combine dataframes
  df = full_join(out2, out4, by=c("SampleID", "Phylum"))
  #set NA values to zero
  df[is.na(df)] <- 0
  #get viability ratios
  df$viability = df$PMARelAbun / df$NoPMARelAbun
  
  #if an ASV has zero relative abundance in untreated samples omit them from the viability ratios
  df$viability[df$NoPMARelAbun == 0] <- NA
  
  #Step Four: #4. Plot viability ratios for the whole communities.
  
  #change viability ratios greater than 1 to 1
  df$viability[df$viability > 1] <- 1
  
  #change viability ratios NaN to NA
  df$viability[is.nan(df$viability)] <- NA
  
  #add Habitat data
  df$Habitat <- data$Habitat[match(df$SampleID, data$SampleID)]
  
  #set habitat as factor
  df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  #find the 15 most abundant Phyla
  top.phy <- top_taxa(noPMA, tax_level = "Phylum", 
                      n_taxa = 15)
  top15 = top.phy[[2]]
  top15$ASV = top15$taxid
  
  #remove any ASVs not in the top 15
  top15 = top15[top15$tax_rank < 16,]
  
  #only keep abundant ASVs
  df.chop = df[df$Phylum %in% top15$Phylum,]
  
  #give ASVs simple names
  # df.chop.id <- transform(df.chop, id=match(ASV, unique(ASV)))
  # df.chop.id = df.chop
  df.chop.id = full_join(df.chop, top15[,c("tax_rank", "Phylum")], by="Phylum")
  
  #give them a new name
  df.chop.id$name = paste0(df.chop.id$tax_rank, ": ", df.chop.id$Phylum)
  
  #make sure the names are in the right order
  df.chop.id$name <- factor(df.chop.id$name, levels=c(mixedsort(unique(df.chop.id$name))))
  
  #only keep abundant phyla
  #df.chop = df[df$Phylum %in% top.phy[[2]]$Phylum,]
  
  # plot as heat map
  p = ggplot(df.chop.id, aes(as.factor(name), SampleID, fill= viability)) + 
    geom_tile() +
    facet_grid(rows = vars(Habitat), 
               scales = "free", 
               space = "free", 
               labeller = label_wrap_gen(width=10),
               switch = "y")+
    scale_fill_gradient(low="grey", high="blue", na.value="white",
                        name="Viability Ratio iDNA/tDNA")+
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right")+
    theme_bw()+
    theme(legend.position="bottom", 
          legend.box = "horizontal",
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=0.05),
          plot.margin = margin(1,1,1,1, "cm"))+
    xlab("Phylum")+
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
  
  return(p)
  
}

pro.phy.plot <- viability_phylum(pro)
pro.phy.plot
euk.phy.plot <- viability_phylum(euk)
euk.phy.plot

#combine plots
# pro.phy.plot2 = pro.phy.plot+theme(legend.position = "none")
# phyla.plot = pro.phy.plot2/ euk.phy.plot + plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 20))

# phyla.plot = ((pro.phy.plot + ggtitle("Prokaryote") + theme(plot.title = element_text(size=35), axis.title.x = element_blank(), legend.position = "none")) 
#       / (euk.phy.plot + ggtitle("Eukaryote") + theme(plot.title = element_text(size=35), axis.title.x = element_blank())))

phyla.plot = ((pro.phy.plot + ggtitle("Prokaryote") + theme(plot.title = element_text(size=40),
                                                              strip.text = element_text(size=20),
                                                              legend.text = element_text(size=25),
                                                              legend.title = element_text(size=25),
                                                              legend.key.size = unit(3,"line"),
                                                              axis.text = element_text(size=18),
                                                              axis.title.y = element_text(size=18),
                                                              axis.title.x = element_blank(),
                                                              legend.position = "none")) 
              / (euk.phy.plot + ggtitle("Microbial Eukaryote") + theme(plot.title = element_text(size=40),
                                                                         strip.text = element_text(size=20),
                                                                         legend.text = element_text(size=25),
                                                                         legend.title = element_text(size=25),
                                                                         legend.key.size = unit(3,"line"),
                                                                         axis.text = element_text(size=18),
                                                                         axis.title.y = element_text(size=18),
                                                                         axis.title.x = element_blank())))

pdf("../results/pma-16S-18S-phyla-dna-ratio.pdf", height=25, width=18)
print(phyla.plot)
dev.off()