# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(ecodist) #for distance() function
library(parallel)
library(svglite)


#Load data 
man.pro = read.csv("../results/prokaryote-mantel-results.csv")
man.euk = read.csv("../results/eukaryote-mantel-results.csv")
man.mm = read.csv("../results/microfauna-mantel-results.csv")

plot.df = rbind(man.euk, man.mm, man.pro)

#make sure our data has the right labels
plot.df$Group[plot.df$Group == "Micrometazoan"] <- "Microfauna"
plot.df$Group[plot.df$Group == "Eukaryote"] <- "Microbial Eukaryote"

#make sure our data are in the right order
plot.df$Group = factor(plot.df$Group, 
                       levels = c("Prokaryote", "Microbial Eukaryote", "Microfauna"))
plot.df$Habitat = factor(plot.df$Habitat,
                         levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
p = ggplot(data=plot.df, aes(x=class.index, y=Mantel.cor)) +
  geom_point(data=plot.df[plot.df$sig=="significant",], aes(x=class.index, y=Mantel.cor), color = "black", size=3, shape=16)+
  geom_point(data=plot.df[plot.df$sig=="non-significant",], aes(x=class.index, y=Mantel.cor), color = "black",size=3, shape=1)+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x = "Phylogenetic distance class", y="Mantel correlation")+
  # ylim(-0.03, 0.2)+
  theme_bw()+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=13),
        legend.position = "bottom")+
  facet_wrap(~ Group + Habitat, scales="free", ncol=2)

p

pdf("../results/mantel-correlogram.pdf", height=15, width=10)
print(p)
dev.off()
