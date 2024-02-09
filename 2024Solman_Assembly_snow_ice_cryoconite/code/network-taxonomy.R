#Script breakdown
#Get taxonomy of nodes + report

#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(Hmisc)
library(igraph)
library(vegan) #rarecurve function
library(dplyr) #for summarising data 
library(ggpubr)
library(car) #levene test
library(plyr)
library(cowplot)
library(tibble) #rownames to column
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(microViz) #for tax_filter
source("00-solman-functions.R")
library(patchwork)

#Import networks
#positive and negative with negative values to absolute values for calculating topological features
g.sn1 = read.graph("../results/full-snow-network-for-calc-properties.gml", format="gml")
g.sp1 = read.graph("../results/full-springice-network-for-calc-properties.gml", format="gml")
g.sm1 = read.graph("../results/full-summerice-network-for-calc-properties.gml", format="gml")
g.cr1 = read.graph("../results/full-cryoconite-network-for-calc-properties.gml", format="gml")

#positive and negative connections
g.sn2 = read.graph("../results/positive-and-negative-snow-network.gml", format="gml")
g.sp2 = read.graph("../results/positive-and-negative-springice-network.gml", format="gml")
g.sm2 = read.graph("../results/positive-and-negative-summerice-network.gml", format="gml")
g.cr2 = read.graph("../results/positive-and-negative-cryoconite-network.gml", format="gml")

#positive connections only
g.sn3 = read.graph("../results/positive-snow-network.gml", format="gml")
g.sp3 = read.graph("../results/positive-springice-network.gml", format="gml")
g.sm3 = read.graph("../results/positive-summerice-network.gml", format="gml")
g.cr3 = read.graph("../results/positive-cryoconite-network.gml", format="gml")

#combine for analysis
g.sn <- list(g.sn1, g.sn2, g.sn3)
g.sp <- list(g.sp1, g.sp2, g.sp3)
g.sm <- list(g.sm1, g.sm2, g.sm3)
g.cr <- list(g.cr1, g.cr2, g.cr3)

#get combined phyloseq object
ps <- readRDS("../results/combined-phylo-for-network-analysis.rds")

#2. Get taxonomy of nodes + report

# ps = ps
# node_df = sn.node

network_tax <- function(ps, node_df){
  
  df <- node_df
  
  #get taxa info of ASVs
  taxonomy <- data.frame(tax_table(ps))
  network_taxa <- subset(taxonomy, rownames(taxonomy) %in% rownames(df))
  
  x = network_taxa %>% 
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
  
  network_taxa = y
  
  out.dom <- network_taxa %>%
    dplyr::group_by(Domain) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  #what percentage are micrometazoans
  out.met = network_taxa[network_taxa$Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa"),]
  met.perc = round((nrow(out.met)/nrow(network_taxa))*100, 2)
  
  out.phy <- network_taxa %>%
    dplyr::group_by(Phylum) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  out.cla <- network_taxa %>%
    dplyr::group_by(Class) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  out.ord <- network_taxa %>%
    dplyr::group_by(Order) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  out.fam <- network_taxa %>%
    dplyr::group_by(Family) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  out.gen <- network_taxa %>%
    dplyr::group_by(Genus) %>% 
    dplyr::summarise(n= n()) %>%
    mutate(perc=n/sum(n)*100)
  
  return(res.list = list(network_taxa, out.dom, out.phy, out.cla, out.ord, out.fam, out.gen, met.perc, out.met$Genus))
  
}

#DOMAIN BREAKDOWN
#Snow Domain
sn.net.tax = network_tax(ps, sn.node)
sn.df.dom1 = sn.net.tax[[2]]
sn.df.dom2 <- sn.df.dom1[order(sn.df.dom1$perc, decreasing = TRUE), ]
sn.df.dom2$Habitat = "Snow"

#Spring ice Domain
sp.net.tax = network_tax(ps, sp.node)
sp.df.dom1 = sp.net.tax[[2]]
sp.df.dom2 <- sp.df.dom1[order(sp.df.dom1$perc, decreasing = TRUE), ]
sp.df.dom2$Habitat = "Spring Ice"

#Summer ice Domain
sm.net.tax = network_tax(ps, sm.node)
sm.df.dom1 = sm.net.tax[[2]]
sm.df.dom2 <- sm.df.dom1[order(sm.df.dom1$perc, decreasing = TRUE), ]
sm.df.dom2$Habitat = "Summer Ice"

#cryoconite Domain
cr.net.tax = network_tax(ps, cr.node)
cr.df.dom1 = cr.net.tax[[2]]
cr.df.dom2 <- cr.df.dom1[order(cr.df.dom1$perc, decreasing = TRUE), ]
cr.df.dom2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.dom2, sp.df.dom2, sm.df.dom2, cr.df.dom2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-domains.csv")

#PHYLUM BREAKDOWN

#Snow phylum
sn.df.phy1 = sn.net.tax[[3]]
sn.df.phy2 <- sn.df.phy1[order(sn.df.phy1$perc, decreasing = TRUE), ]

#spring ice phylum
sp.df.phy1 = sp.net.tax[[3]]
sp.df.phy2 <- sp.df.phy1[order(sp.df.phy1$perc, decreasing = TRUE), ]

#summer ice phylum
sm.df.phy1 = sm.net.tax[[3]]
sm.df.phy2 <- sm.df.phy1[order(sm.df.phy1$perc, decreasing = TRUE), ]

#cryoconite phylum
cr.df.phy1 = cr.net.tax[[3]]
cr.df.phy2 <- cr.df.phy1[order(cr.df.phy1$perc, decreasing = TRUE), ]


#CLASS BREAKDOWN

#Snow class
sn.df.cla1 = sn.net.tax[[4]]
sn.df.cla2 <- sn.df.cla1[order(sn.df.cla1$perc, decreasing = TRUE), ]

#spring ice class
sp.df.cla1 = sp.net.tax[[4]]
sp.df.cla2 <- sp.df.cla1[order(sp.df.cla1$perc, decreasing = TRUE), ]

#summer ice class
sm.df.cla1 = sm.net.tax[[4]]
sm.df.cla2 <- sm.df.cla1[order(sm.df.cla1$perc, decreasing = TRUE), ]

#cryoconite class
cr.df.cla1 = cr.net.tax[[4]]
cr.df.cla2 <- cr.df.cla1[order(cr.df.cla1$perc, decreasing = TRUE), ]

#ORDER BREAKDOWN

#Snow order
sn.df.ord1 = sn.net.tax[[5]]
sn.df.ord2 <- sn.df.ord1[order(sn.df.ord1$perc, decreasing = TRUE), ]

#spring ice order
sp.df.ord1 = sp.net.tax[[5]]
sp.df.ord2 <- sp.df.ord1[order(sp.df.ord1$perc, decreasing = TRUE), ]

#summer ice order
sm.df.ord1 = sm.net.tax[[5]]
sm.df.ord2 <- sm.df.ord1[order(sm.df.ord1$perc, decreasing = TRUE), ]

#cryoconite order
cr.df.ord1 = cr.net.tax[[5]]
cr.df.ord2 <- cr.df.ord1[order(cr.df.ord1$perc, decreasing = TRUE), ]

#FAMILY BREAKDOWN

#Snow family
sn.df.fam1 = sn.net.tax[[6]]
sn.df.fam2 <- sn.df.fam1[order(sn.df.fam1$perc, decreasing = TRUE), ]

#spring ice family
sp.df.fam1 = sp.net.tax[[6]]
sp.df.fam2 <- sp.df.fam1[order(sp.df.fam1$perc, decreasing = TRUE), ]

#summer ice family
sm.df.fam1 = sm.net.tax[[6]]
sm.df.fam2 <- sm.df.fam1[order(sm.df.fam1$perc, decreasing = TRUE), ]

#cryoconite family
cr.df.fam1 = cr.net.tax[[6]]
cr.df.fam2 <- cr.df.fam1[order(cr.df.fam1$perc, decreasing = TRUE), ]

#GENUS BREAKDOWN

#Snow genus
sn.df.gen1 = sn.net.tax[[7]]
sn.df.gen2 <- sn.df.gen1[order(sn.df.gen1$perc, decreasing = TRUE), ]

#spring ice genus
sp.df.gen1 = sp.net.tax[[7]]
sp.df.gen2 <- sp.df.gen1[order(sp.df.gen1$perc, decreasing = TRUE), ]

#summer ice genus
sm.df.gen1 = sm.net.tax[[7]]
sm.df.gen2 <- sm.df.gen1[order(sm.df.gen1$perc, decreasing = TRUE), ]

#cryoconite genus
cr.df.gen1 = cr.net.tax[[7]]
cr.df.gen2 <- cr.df.gen1[order(cr.df.gen1$perc, decreasing = TRUE), ]