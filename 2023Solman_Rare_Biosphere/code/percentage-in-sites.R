# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

# library(ggpubr)
#install.packages("ggstatsplot")
# library(ggstatsplot)
source("00-solman-functions.R")
# library(ggcorrplot)
#library(fdrtool)
# library(psych)
# library(cowplot)
# #devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)
library(stringr) #for subsetting strings in lat long function
library(fossil) #for earth.dist function
library(vegan) #for community dissimilarity matrix
library(tibble) #for rownames to column
library(dplyr) #coalesce function + rows update function
library(reshape2) #for melt function
library(funrar)
library(phyloseq)
library(ggplot2)
#install.packages("ggpmisc")
library(ggpmisc)
library(picante) #for faith's pd
library(cowplot)

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) >= 1, TRUE)

#Arctic prokaryotes
pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds") 
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds") 

#Antarctic prokaryotes
pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun.rds") 
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int.rds") 
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare.rds") 

#Arctic eukaryotes
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds") 
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds") 

#Antarctic eukaryotes
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun.rds") 
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int.rds") 
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare.rds") 


#find the percentage of abundant, intermediate and rare ASVs within samples

#get a dataframe with each ASV, it's abundance classification, the number of samples it is found in, the percentage of total samples that represents

# ps = euk.ant
# ps.abun = euk.ant.abun
# ps.int = euk.ant.int
# ps.rar = euk.ant.rare
# pole = "Antarctic"
# group = "Eukaryote"

PercInSamples <- function(ps, ps.abun, ps.int, ps.rar, pole, group){
  
  com = data.frame(otu_table(ps)) #get count table
  com[com > 1] = 1 #change to presence/absence dataframe
  max(com) #make sure it worked - this should equal 1
  x = rowSums(com)
  df = data.frame(ASV=row.names(com), NumSites = x, PercSites=x/ncol(com)*100)
  abun = data.frame(ASV=row.names(data.frame(otu_table(ps.abun))), Subcommunity="Abundant")
  int = data.frame(ASV=row.names(data.frame(otu_table(ps.int))), Subcommunity="Intermediate")
  rar = data.frame(ASV=row.names(data.frame(otu_table(ps.rar))), Subcommunity="Rare")
  y = rbind(abun, int, rar)
  final.df = left_join(df, y)
  
  final.df = final.df[final.df$NumSites > 0,]
  
  #how many of each subcommunity are present in at least 50% of samples
  # sub.50 = final.df[final.df$PercSites >= 50,]
  # df1 = data.frame(table(factor(final.df[final.df$PercSites >= 50,]$Subcommunity)))
  # df1$PercentageOfSites = ">=50%"
  # PercentageOfASVs = c(df1$Freq[[1]]/nrow(abun), df1$Freq[[2]]/nrow(int))
  
  
  subcom = c("Abundant", "Intermediate", "Rare")
  num.25.keep = vector()
  num.50.keep = vector()
  num.75.keep = vector()
  num.99.keep = vector()
  
  for (i in 1:length(subcom)){
    
    com.class = subcom[[i]]
    
    cut.com = final.df[final.df$Subcommunity == com.class,]
    
    #percentage of ASV in over 25% of samples
    num.25.keep[[i]] = nrow(cut.com[cut.com$PercSites > 25,])/nrow(cut.com)
    
    #percentage of ASV in over 50% of samples
    num.50.keep[[i]] = nrow(cut.com[cut.com$PercSites > 50,])/nrow(cut.com)
    
    #percentage of ASV in over 75% of samples
    num.75.keep[[i]] = nrow(cut.com[cut.com$PercSites > 75,])/nrow(cut.com)
    
    #percentage of ASV in over 99% of samples
    num.99.keep[[i]] = nrow(cut.com[cut.com$PercSites > 99,])/nrow(cut.com)
    
  }
  
  #make into a data frame
  
  df.to.print = data.frame(Group = group, Pole = pole, Subcommunity = subcom, Perc25 = round(num.25.keep*100, 2), Perc50 = round(num.50.keep*100, 2), Perc75 = round(num.75.keep*100, 2), Perc99 = round(num.99.keep*100, 2))
  
  return(df.to.print)
  
  
  
}


out.df.1 = PercInSamples(pro.arc, pro.arc.abun, pro.arc.int, pro.arc.rare, "Arctic", "Prokaryote")
out.df.2 = PercInSamples(pro.ant, pro.ant.abun, pro.ant.int, pro.ant.rare, "Antarctic", "Prokaryote")
out.df.3 = PercInSamples(euk.arc, euk.arc.abun, euk.arc.int, euk.arc.rare, "Arctic", "Eukaryote")
out.df.4 = PercInSamples(euk.ant, euk.ant.abun, euk.ant.int, euk.ant.rare, "Antarctic", "Eukaryote")

final.out.df = rbind(out.df.1, out.df.2, out.df.3, out.df.4)

#save results

write.csv(final.out.df, "../results/percentage-in-sites.csv")
