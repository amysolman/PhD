
rm(list=ls())
graphics.off()

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

# 2. Import data
ps.pro <- readRDS("../data/16S-phylo-object-rel.rds") 
ps.euk <- readRDS("../data/18S-phylo-object-rel.rds") 
#remove blanks
ps.euk = subset_samples(ps.euk, SampleType != "Control")
ps.euk = filter_taxa(ps.euk, function(x) sum(x) > 0, TRUE)

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


#PCOA function

# ps = ps.euk
# d = "bray" #options unifrac, wunifrac, bray

PCoA_function <- function(ps, d){
  
  #perform ordination
  ord <- ordinate(ps, method="PCoA", distance=d)
  
  #extract data for plotting
  positions <- ord$vectors[,1:2]
  colnames(positions) <- c("pcoa1", "pcoa2")
  
  #get percentages explained by the first 2 axis
  percent_explained <- c(100*sum(ord$values$Relative_eig[1]), 100*sum(ord$values$Relative_eig[2]))
  pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
  
  #create axis labels
  labels <- c(glue("PCoA Axis 1 ({pretty_pe[1]}%)"),
              glue("PCoA Axis 2 ({pretty_pe[2]}%)"))
  
  #give sample names as row names and add regions and poles
  x = data.frame(sample_data(ps))
  data2plot = positions %>%
    as_tibble(rownames = "samples") %>%
    add_column(data.frame(sample_data(ps))$SampleType, data.frame(sample_data(ps))$Site, data.frame(sample_data(ps))$Day)
  names(data2plot) = c("Sample", "pcoa1", "pcoa2",  "SampleType", "Site", "Day")
  
  #default legend order
  #data2plot$Habitat <- factor(data2plot$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  #plot the data
  p1 = ggplot(data=data2plot) +
    geom_point(data=data2plot, aes(x=pcoa1, y=pcoa2, fill=SampleType, shape=Site),
               alpha=.5, size=7)+
    labs(x=labels[1], y=labels[2])+
    theme_bw()+
    scale_shape_manual(values = c(21,23,24))+
    #scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
    guides(fill=guide_legend(override.aes=list(shape=21), order = 1),
           shape = guide_legend(order = 2))+
    theme(legend.position="bottom", legend.title = element_blank(),
          legend.text = element_text(size = 14))
  
  print(p1)
  
  return(p1)
  
}

#16S
pro.pcoa.uni <- PCoA_function(ps.pro, "unifrac")
pro.pcoa.wuni <- PCoA_function(ps.pro, "wunifrac")
pro.pcoa.bray <- PCoA_function(ps.pro, "bray")

#18S
euk.pcoa.uni <- PCoA_function(ps.euk, "unifrac")
euk.pcoa.wuni <- PCoA_function(ps.euk, "wunifrac")
euk.pcoa.bray <- PCoA_function(ps.euk, "bray")

#18S with top 2 most abundant ASVs removed
euk.rm.pcoa.uni <- PCoA_function(ps.euk.rm, "unifrac")
euk.rm.pcoa.wuni <- PCoA_function(ps.euk.rm, "wunifrac")
euk.rm.pcoa.bray <- PCoA_function(ps.euk.rm, "bray")

#Note PCoA was used here (instead of NMDS) because there wasn't sufficient data in
#the 16S dataset for NMDS analysis and we want to use the same method for 16S + 18S.


#put all the plots together
final.p = pro.pcoa.wuni / euk.pcoa.wuni + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))

pdf("../results/wunifrac-pcoa.pdf", height=8, width=10)
print(final.p)
dev.off()

#18S PCoA with Weighted Unifrac distances (top 2 ASVs removed)
pdf("../results/18S-rm-wunifrac-pcoa.pdf", height=8, width=10)
print(euk.rm.pcoa.wuni)
dev.off()
