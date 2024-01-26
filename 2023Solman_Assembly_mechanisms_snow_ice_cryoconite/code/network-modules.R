rm(list=ls())
graphics.off()

library(tibble)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(phyloseq)
#install.packages("ggtern")
library(ggtern)
source("00-solman-functions.R")
library(RColorBrewer) #for our pie chart colours
# install.packages("randomcoloR")
library(randomcoloR)

#Import gephi data tables and add colour column for taxonomy on node df 
#and colour column for positive/negative correlations on edge df.


#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#Load graphics dataframes from Gephi step
#Snow
sn.node <- read.csv("../results/network-nodes-snow.csv")
sn.edge <- read.csv("../results/network-edges-snow.csv")

#Spring Ice
sp.node <- read.csv("../results/network-nodes-spring-ice.csv")
sp.edge <- read.csv("../results/network-edges-spring-ice.csv")

#Summer Ice
sm.node <- read.csv("../results/network-nodes-summer-ice.csv")
sm.edge <- read.csv("../results/network-edges-summer-ice.csv")

#cryoconite
cr.node <- read.csv("../results/network-nodes-summer-ice.csv")
cr.edge <- read.csv("../results/network-edges-summer-ice.csv")

# edge.df = sn.edge
# node.df = sn.node
# mod_classes = c(as.numeric(names(sort(table(sn.node$modularity_class), decreasing = TRUE)[1:6])))

module_analysis <- function(edge.df, node.df, mod_classes){
  
  #find out which ASVs are present in each modules and how many edges they have
  x = data.frame(table(edge.df$Source))
  x2 = data.frame(table(edge.df$Target))
  df = rbind(x, x2)
  df$Var1 = as.character(df$Var1)
  names(df) = c("Id", "Edges")
  
  #merge where ASV is the same
  df2 = df %>%
    dplyr::group_by(Id) %>%
    dplyr::summarise(Edges = sum(Edges))
  
  #get ASV tax info
  node.df$Id = as.character(node.df$Id)
  df3 = full_join(df2, node.df, by="Id")
  
  #subset to module 1
  mode.1 = df3[df3$modularity_class == mod_classes[1],]
  #which is the most connected ASV?
  hub.1 = subset(mode.1, Edges == max(Edges))
  hub.1$Module = 1
  
  #subset to module 2
  mode.2 = df3[df3$modularity_class == mod_classes[2],]
  #which is the most connected ASV?
  hub.2 = subset(mode.2, Edges == max(Edges))
  hub.2$Module = 2
  
  #subset to module 3
  mode.3 = df3[df3$modularity_class == mod_classes[3],]
  #which is the most connected ASV?
  hub.3 = subset(mode.3, Edges == max(Edges))
  hub.3$Module = 3
  
  #subset to module 4
  mode.4 = df3[df3$modularity_class == mod_classes[4],]
  #which is the most connected ASV?
  hub.4 = subset(mode.4, Edges == max(Edges))
  hub.4$Module = 4
  
  #subset to module 5
  mode.5 = df3[df3$modularity_class == mod_classes[5],]
  #which is the most connected ASV?
  hub.5 = subset(mode.5, Edges == max(Edges))
  hub.5$Module = 5
  
  #subset to module 6
  mode.6 = df3[df3$modularity_class == mod_classes[6],]
  #which is the most connected ASV?
  hub.6 = subset(mode.6, Edges == max(Edges))
  hub.6$Module = 6
  
  res = rbind(hub.1, hub.2, hub.3, hub.4, hub.5, hub.6)
  
  return(res)
  
}


#run function

#Snow
#get the numbers of the top 6 modules
mods.snow = c(as.numeric(names(sort(table(sn.node$modularity_class), decreasing = TRUE)[1:6])))
sn.mod <- module_analysis(sn.edge, sn.node, mod_classes = mods.snow)
sn.mod$Habitat = "Snow"

#Spring Ice
mods.sp = c(as.numeric(names(sort(table(sp.node$modularity_class), decreasing = TRUE)[1:6])))
sp.mod <- module_analysis(sp.edge, sp.node, mod_classes = mods.sp)
sp.mod$Habitat = "Spring Ice"

#Summer Ice
mods.sm = c(as.numeric(names(sort(table(sm.node$modularity_class), decreasing = TRUE)[1:6])))
sm.mod <- module_analysis(sm.edge, sm.node, mod_classes = mods.sm)
sm.mod$Habitat = "Summer Ice"

#Cryoconite
mods.cr = c(as.numeric(names(sort(table(cr.node$modularity_class), decreasing = TRUE)[1:6])))
cr.mod <- module_analysis(cr.edge, cr.node, mod_classes = mods.cr)
cr.mod$Habitat = "Cryoconite"

#full results
full.df = rbind(sn.mod, sp.mod, sm.mod, cr.mod)



