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

####Alpha Diversity
alpha_function <- function(phylo){
  
  alpha = estimate_richness(phylo)
  alpha$Sample = sample_names(phylo)
  alpha$Treatment = data.frame(sample_data(phylo))$Treatment
  alpha$Habitat = data.frame(sample_data(phylo))$Habitat
  alpha$Habitat = factor(alpha$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  #change names of treatments
  alpha$Treatment[alpha$Treatment == "No-PMA"] <- "tDNA"
  alpha$Treatment[alpha$Treatment == "PMA"] <- "iDNA"
  
  p1 <- ggboxplot(alpha, x = "Treatment", y = "Observed",
                  color = "Treatment", palette = "jco",
                  add = "jitter")+
    facet_grid(cols = vars(Habitat))+
    theme(legend.position = "none")
  p1
  #  Add p-value
  p1 = p1  + stat_compare_means()
  p1
  
  
  
  p2 <- ggboxplot(alpha, x = "Treatment", y = "Shannon",
                  color = "Treatment", palette = "jco",
                  add = "jitter")+
    facet_grid(cols = vars(Habitat))+
    theme(legend.position = "none")
  p2
  #  Add p-value
  p2 = p2  + stat_compare_means()
  p2
  
  
  p3 <- ggboxplot(alpha, x = "Treatment", y = "Simpson",
                  color = "Treatment", palette = "jco",
                  add = "jitter")+
    facet_grid(cols = vars(Habitat))+
    theme(legend.position = "none")
  #  Add p-value
  p3 = p3  + stat_compare_means()
  p3
  
  
  p4 <- ggboxplot(alpha, x = "Treatment", y = "Chao1",
                  color = "Treatment", palette = "jco",
                  add = "jitter")+
    facet_grid(cols = vars(Habitat))+
    theme(legend.position = "none")
  #  Add p-value
  p4 = p4  + stat_compare_means()
  p4
  
  return(list(p1,p2,p3,p4))
  
}

#16S
pro.plots = alpha_function(pro)
pro.plots[[1]] #observed
pro.plots[[2]] #shannon
pro.plots[[3]] #simpson
pro.plots[[4]] #chao1

#18S
euk.plots = alpha_function(euk)
euk.plots[[1]] #observed
euk.plots[[2]] #shannon
euk.plots[[3]] #simpson
euk.plots[[4]] #chao1


#Combine plots

p = pro.plots[[1]] / pro.plots[[2]] / pro.plots[[3]] / pro.plots[[4]]
p

pro.p = plot_grid(pro.plots[[1]] + theme(plot.title = element_text(size=40),
                                         strip.text = element_text(size=25),
                                         legend.text = element_text(size=20),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(3,"line"),
                                         axis.text = element_text(size=18),
                                         axis.title = element_text(size=18)),
                  pro.plots[[2]] + theme(plot.title = element_text(size=40),
                                         strip.text = element_text(size=25),
                                         legend.text = element_text(size=20),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(3,"line"),
                                         axis.text = element_text(size=18),
                                         axis.title = element_text(size=18)),
                  pro.plots[[3]] + theme(plot.title = element_text(size=40),
                                         strip.text = element_text(size=25),
                                         legend.text = element_text(size=20),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(3,"line"),
                                         axis.text = element_text(size=18),
                                         axis.title = element_text(size=18)),
                  pro.plots[[4]] + theme(plot.title = element_text(size=40),
                                         strip.text = element_text(size=25),
                                         legend.text = element_text(size=20),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(3,"line"),
                                         axis.text = element_text(size=18),
                                         axis.title = element_text(size=18)),
                  ncol=1,
                  align=c("v"))
pro.p

#add title
# now add the title
title <- ggdraw() + draw_label(
    "A: Prokaryote",
    x = 0,
    hjust = 0,
    size=40
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 50)
  )

pro.p = plot_grid(
  title, pro.p,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

euk.p = plot_grid(euk.plots[[1]] + theme(plot.title = element_text(size=40),
                                       strip.text = element_text(size=25),
                                       legend.text = element_text(size=20),
                                       legend.title = element_blank(),
                                       legend.key.size = unit(3,"line"),
                                       axis.text = element_text(size=18),
                                       axis.title = element_text(size=18)),
                  euk.plots[[2]]+ theme(plot.title = element_text(size=40),
                                        strip.text = element_text(size=25),
                                        legend.text = element_text(size=20),
                                        legend.title = element_blank(),
                                        legend.key.size = unit(3,"line"),
                                        axis.text = element_text(size=18),
                                        axis.title = element_text(size=18)),
                  euk.plots[[3]]+ theme(plot.title = element_text(size=40),
                                        strip.text = element_text(size=25),
                                        legend.text = element_text(size=20),
                                        legend.title = element_blank(),
                                        legend.key.size = unit(3,"line"),
                                        axis.text = element_text(size=18),
                                        axis.title = element_text(size=18)),
                  euk.plots[[4]]+ theme(plot.title = element_text(size=40),
                                        strip.text = element_text(size=25),
                                        legend.text = element_text(size=20),
                                        legend.title = element_blank(),
                                        legend.key.size = unit(3,"line"),
                                        axis.text = element_text(size=18),
                                        axis.title = element_text(size=18)),
                  ncol=1,
                  align=c("v"))
euk.p

# now add the title
title <- ggdraw() + draw_label(
  "B: Microbial Eukaryote",
  x = 0,
  hjust = 0,
  size=40
) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 50),
  )

euk.p = plot_grid(
  title, euk.p,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)
euk.p

p = pro.p+ euk.p
p

pdf("../results/pma-alpha-comparisons.pdf", width=18, height=25)
print(p)
dev.off()