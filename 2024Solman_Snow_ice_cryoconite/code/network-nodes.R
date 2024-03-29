#Script breakdown
#1. Calculate node-level properties
#2. Get taxonomy of nodes + report
#3. Find differences between subcommunities
#4. Find connections between subcommunities
#5. Identify keystone taxa

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

#positive and negative networks
g.sn2 = read.graph("../results/positive-and-negative-snow-network.gml", format="gml")
g.sp2 = read.graph("../results/positive-and-negative-springice-network.gml", format="gml")
g.sm2 = read.graph("../results/positive-and-negative-summerice-network.gml", format="gml")
g.cr2 = read.graph("../results/positive-and-negative-cryoconite-network.gml", format="gml")

#positive only
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

#1. Calculate node-level properties

# g = g.sn[[1]]

node_level_properties <- function(g, habitat){
  
  cor_g <- g
  
  #create a dataframe with 5 columns and a row for each node (ASV)
  df<-as.data.frame(matrix(NA,ncol=5,nrow=length(degree(cor_g))))
  
  #make ASV IDs row names of df
  rownames(df)<-names(degree(cor_g))
  
  #name the colums the node-level topological features
  colnames(df)<-c("degree","betweenness","closeness","eigenvector","habitat")
  
  #get betweenness
  btw<-betweenness(cor_g)
  #get closeness centrality
  cls<-closeness(cor_g)
  #get eigenvector centrality
  egv<-evcent(cor_g)
  
  #put the topological features into the dataframe
  df[,1]<-degree(cor_g)
  df[,2]<-btw
  df[,3]<-cls
  df[,4]<-egv$vector
  df[,5]<-habitat
  
  return(df)
  
}

sn.node = node_level_properties(g.sn[[1]], "Snow")
sp.node = node_level_properties(g.sp[[1]], "Spring Ice")
sm.node = node_level_properties(g.sm[[1]], "Summer Ice")
cr.node = node_level_properties(g.cr[[1]], "Cryoconite")

#combine and export node-level properties
x = rbind(sn.node, sp.node, sm.node, cr.node)
write.csv(x, "../results/node-properties.csv")

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

test = cr.net.tax[[1]]

#Snow phylum
sn.df.phy1 = sn.net.tax[[3]]
sn.df.phy2 <- sn.df.phy1[order(sn.df.phy1$perc, decreasing = TRUE), ]
sn.df.phy2$Habitat = "Snow"

#spring ice phylum
sp.df.phy1 = sp.net.tax[[3]]
sp.df.phy2 <- sp.df.phy1[order(sp.df.phy1$perc, decreasing = TRUE), ]
sp.df.phy2$Habitat = "Spring Ice"

#summer ice phylum
sm.df.phy1 = sm.net.tax[[3]]
sm.df.phy2 <- sm.df.phy1[order(sm.df.phy1$perc, decreasing = TRUE), ]
sm.df.phy2$Habitat = "Summer Ice"

#cryoconite phylum
cr.df.phy1 = cr.net.tax[[3]]
cr.df.phy2 <- cr.df.phy1[order(cr.df.phy1$perc, decreasing = TRUE), ]
cr.df.phy2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.phy2, sp.df.phy2, sm.df.phy2, cr.df.phy2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-phyla.csv")

#CLASS BREAKDOWN

#Snow class
sn.df.cla1 = sn.net.tax[[4]]
sn.df.cla2 <- sn.df.cla1[order(sn.df.cla1$perc, decreasing = TRUE), ]
sn.df.cla2$Habitat = "Snow"

#spring ice class
sp.df.cla1 = sp.net.tax[[4]]
sp.df.cla2 <- sp.df.cla1[order(sp.df.cla1$perc, decreasing = TRUE), ]
sp.df.cla2$Habitat = "Spring Ice"

#summer ice class
sm.df.cla1 = sm.net.tax[[4]]
sm.df.cla2 <- sm.df.cla1[order(sm.df.cla1$perc, decreasing = TRUE), ]
sm.df.cla2$Habitat = "Summer Ice"

#cryoconite class
cr.df.cla1 = cr.net.tax[[4]]
cr.df.cla2 <- cr.df.cla1[order(cr.df.cla1$perc, decreasing = TRUE), ]
cr.df.cla2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.cla2, sp.df.cla2, sm.df.cla2, cr.df.cla2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-class.csv")

#ORDER BREAKDOWN

#Snow order
sn.df.ord1 = sn.net.tax[[5]]
sn.df.ord2 <- sn.df.ord1[order(sn.df.ord1$perc, decreasing = TRUE), ]
sn.df.ord2$Habitat = "Snow"

#spring ice order
sp.df.ord1 = sp.net.tax[[5]]
sp.df.ord2 <- sp.df.ord1[order(sp.df.ord1$perc, decreasing = TRUE), ]
sp.df.ord2$Habitat = "Spring Ice"

#summer ice order
sm.df.ord1 = sm.net.tax[[5]]
sm.df.ord2 <- sm.df.ord1[order(sm.df.ord1$perc, decreasing = TRUE), ]
sm.df.ord2$Habitat = "Summer Ice"

#cryoconite order
cr.df.ord1 = cr.net.tax[[5]]
cr.df.ord2 <- cr.df.ord1[order(cr.df.ord1$perc, decreasing = TRUE), ]
cr.df.ord2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.ord2, sp.df.ord2, sm.df.ord2, cr.df.ord2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-order.csv")

#FAMILY BREAKDOWN

#Snow family
sn.df.fam1 = sn.net.tax[[6]]
sn.df.fam2 <- sn.df.fam1[order(sn.df.fam1$perc, decreasing = TRUE), ]
sn.df.fam2$Habitat = "Snow"

#spring ice family
sp.df.fam1 = sp.net.tax[[6]]
sp.df.fam2 <- sp.df.fam1[order(sp.df.fam1$perc, decreasing = TRUE), ]
sp.df.fam2$Habitat = "Spring Ice"

#summer ice family
sm.df.fam1 = sm.net.tax[[6]]
sm.df.fam2 <- sm.df.fam1[order(sm.df.fam1$perc, decreasing = TRUE), ]
sm.df.fam2$Habitat = "Summer Ice"

#cryoconite family
cr.df.fam1 = cr.net.tax[[6]]
cr.df.fam2 <- cr.df.fam1[order(cr.df.fam1$perc, decreasing = TRUE), ]
cr.df.fam2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.fam2, sp.df.fam2, sm.df.fam2, cr.df.fam2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-family.csv")

#GENUS BREAKDOWN

#Snow genus
sn.df.gen1 = sn.net.tax[[7]]
sn.df.gen2 <- sn.df.gen1[order(sn.df.gen1$perc, decreasing = TRUE), ]
sn.df.gen2$Habitat = "Snow"

#spring ice genus
sp.df.gen1 = sp.net.tax[[7]]
sp.df.gen2 <- sp.df.gen1[order(sp.df.gen1$perc, decreasing = TRUE), ]
sp.df.gen2$Habitat = "Spring Ice"

#summer ice genus
sm.df.gen1 = sm.net.tax[[7]]
sm.df.gen2 <- sm.df.gen1[order(sm.df.gen1$perc, decreasing = TRUE), ]
sm.df.gen2$Habitat = "Summer Ice"

#cryoconite genus
cr.df.gen1 = cr.net.tax[[7]]
cr.df.gen2 <- cr.df.gen1[order(cr.df.gen1$perc, decreasing = TRUE), ]
cr.df.gen2$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.df.gen2, sp.df.gen2, sm.df.gen2, cr.df.gen2)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-genus.csv")

#MICROFAUNA CONTRIBUTION
sn.mf = sn.net.tax[[3]][sn.net.tax[[3]]$Phylum == "Tardigrada" | sn.net.tax[[3]]$Phylum == "Rotifera" | sn.net.tax[[3]]$Phylum == "Nematozoa",]
sn.mf$Habitat = "Snow"
sp.mf = sp.net.tax[[3]][sp.net.tax[[3]]$Phylum == "Tardigrada" | sp.net.tax[[3]]$Phylum == "Rotifera" | sp.net.tax[[3]]$Phylum == "Nematozoa",]
sp.mf$Habitat = "Spring Ice"
sm.mf = sm.net.tax[[3]][sm.net.tax[[3]]$Phylum == "Tardigrada" | sm.net.tax[[3]]$Phylum == "Rotifera" | sm.net.tax[[3]]$Phylum == "Nematozoa",]
sm.mf$Habitat = "Summer Ice"
cr.mf = cr.net.tax[[3]][cr.net.tax[[3]]$Phylum == "Tardigrada" | cr.net.tax[[3]]$Phylum == "Rotifera" | cr.net.tax[[3]]$Phylum == "Nematozoa",]
cr.mf$Habitat = "Cryoconite"

#combine and export
x = rbind(sn.mf, sp.mf, sm.mf, cr.mf)
x$perc = round(x$perc,2)
write.csv(x, "../results/network-microfauna.csv")

#3. Find differences between subcommunities

# df1 = sn.node
# df2 = sp.node
# df3 = sm.node
# df4 = cr.node

signif_dif_nodes <- function(df1, df2, df3, df4){
  
  #bind the dfs into one
  df = rbind(df1, df2, df3, df4)
  
  # I reorder the groups order : I change the order of the factor 
  df$habitat <- factor(df$habitat , levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  #Degree
  my_comparisons <- list( c("Snow", "Spring Ice"), c("Snow", "Summer Ice"), c("Snow", "Cryoconite"),
                          c("Spring Ice", "Summer Ice"), c("Spring Ice", "Cryoconite"), c("Summer Ice", "Cryoconite"))
  
  #log10-transform the data
  df$log10degree = log10(df$degree)
  
  p1 <- ggboxplot(df, x = "habitat", y = "log10degree",
                  fill = "habitat")+ 
    stat_compare_means(comparisons = my_comparisons, size=5, p.adjust.method = "bonferroni")+
    #scale_fill_brewer(palette="Set1")
    scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
    ylab("log10 Degree")+
    theme(axis.title.y=element_text(size=20))
  p1 = p1 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p1)$y$range$range[2]*1.1, size=5)  # Add global p-value
  p1
  
  #Between
  
  #log10-transform the data
  df$log10betweenness = log10(df$betweenness)
  
  p2 <- ggboxplot(df, x = "habitat", y = "log10betweenness",
                  fill = "habitat")+
    #ylim(0, 3000)+
    stat_compare_means(comparisons = my_comparisons, size=5, p.adjust.method = "bonferroni")+
    #scale_fill_brewer(palette="Set1")
    scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
    ylab("log10 Betweenness")+
    theme(axis.title.y=element_text(size=20))
  p2 = p2 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p2)$y$range$range[2]*1.1, size=5)  # Add global p-value
  p2
  
  #Closeness
  # df.close = df[order(df$closeness, decreasing = TRUE), ]
  # df.close = df.close[27:nrow(df.close),]
  
  #log10-transform the data
  df$log10closeness = log10(df$closeness)
  
  p3 <- ggboxplot(df, x = "habitat", y = "log10closeness",
                  fill = "habitat")+ 
    #ylim(0, 0.6)+
    stat_compare_means(comparisons = my_comparisons, size=5, p.adjust.method = "bonferroni")+
    #scale_fill_brewer(palette="Set1")
    scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
    ylab("log10 Closeness")+
    theme(axis.title.y=element_text(size=20))
  p3 = p3 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p3)$y$range$range[2]*1.1, size=5)  # Add global p-value
  p3
  
  #Eigenvector
  
  #log10-transform the data
  df$log10eigenvector = log10(df$eigenvector)
  
  p4 <- ggboxplot(df, x = "habitat", y = "log10eigenvector",
                  fill= "habitat")+ 
    #ylim(0, 0.01)+
    stat_compare_means(comparisons = my_comparisons, size=5, p.adjust.method = "bonferroni")+ # Add pairwise comparisons p-value
    #scale_fill_brewer(palette="Set1")
    scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
    ylab("log10 Eigenvector Centrality")+
    theme(axis.title.y=element_text(size=20))
  p4 = p4 +
    stat_compare_means(label.y = layer_scales(p4)$y$range$range[2]*1.1, size=5)  # Add global p-value
  p4
  
  
  return(res.list=list(p1, p2, p3, p4))
  
}

sig.plots = signif_dif_nodes(sn.node, sp.node, sm.node, cr.node)

p1 = sig.plots[[1]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
p2 = sig.plots[[2]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
p3 = sig.plots[[3]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
p4 = sig.plots[[4]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

prow = (p1 | p2) / (p3 | p4)
prow

legend_b <- get_legend(sig.plots[[1]] +
                         guides(color = guide_legend(nrow = 1, override.aes = list(size = 3),title=NULL)) +
                         theme(legend.position = "bottom",
                               axis.text.y=element_text(size=25),
                               axis.title.y=element_text(size=25,face="bold"),
                               legend.text = element_text(size=25),
                               legend.title = element_blank())
)

final.p = plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

final.p

#export
pdf("../results/node-level-proporties-boxplot.pdf", height=10, width=13)
print(final.p)
dev.off()

#14. Identify keystone taxa

# ps = ps
# node_df = sn.node
# habitat = "Snow"

get_keys <- function(ps, node_df, habitat){
  
  #Get taxonomy of most highly connected/influential nodes
  tax.tab = data.frame(tax_table(ps))
  
  #make ASV ID rownames a column
  tax.tab <- tibble::rownames_to_column(tax.tab, "ASVID")
  node_df <- tibble::rownames_to_column(node_df, "ASVID")
  
  #use dplyr to bind the dataframes
  res <- inner_join(tax.tab, node_df, by = "ASVID")
  
  #only keep taxa in the top 40% of degree abundance
  #sort by top degree values and extract only the top 10%
  order.degree <- res[order(-res$degree),]
  max.degree <- head(order.degree,0.40*nrow(order.degree))
  
  #only keep taxa in the top 40% of closeness
  #sort by top closeness values and extract only the top 10%
  order.close <- res[order(-res$close),]
  max.close <- head(order.close,0.40*nrow(order.close))
  
  #only keep taxa in the bottom 40% of betweeness
  #sort by bottom betweeness values and extract only the bottom 10%
  order.between <- res[order(res$between),]
  min.between <- head(order.between,0.40*nrow(order.between))

  #Keystone taxa are those that overlap
  k.node = res[res$ASVID %in% max.degree$ASVID & res$ASVID %in% max.close$ASVID & res$ASVID %in% min.between$ASVID,]
  
  write.csv(k.node, paste0("../results/", habitat , "-keystone-taxa.csv"))
  
  return(k.node)
  
}

k.sn <- get_keys(ps, sn.node, "Snow")
k.sp <- get_keys(ps, sp.node, "SpringIce")
k.sm <- get_keys(ps, sm.node, "SummerIce")
k.cr <- get_keys(ps, cr.node, "Cryoconite")

#combine and export
x = rbind(k.sn, k.sp, k.sm, k.cr)
write.csv(x, "../results/network-keystone.csv")

#Which ASVs had the highest number of degrees?
#What was their taxonomy?
#Which ASVs had negative connections?

# g = g.sp[[2]]
# ps = ps
# node_df = sp.node

top_degree_and_neg_corr <- function(g, ps, node_df){
  
  df <- igraph::as_data_frame(g) #extract network data
  
  df2 = as.data.frame(sort(table(df$from), decreasing = TRUE))
  df3 = as.data.frame(sort(table(df$to), decreasing = TRUE))
  
  df4 = rbind(df2, df3)
  
  #combine rows with duplicate ASV IDs
  df5 = df4 %>% 
    group_by(Var1) %>%
    summarise_each(funs(sum))
  
  #get edited tax table
  net.tax <- network_tax(ps, node_df)
  net.tax.df <- net.tax[[1]]
  net.tax.df <- tibble::rownames_to_column(net.tax.df, "ASVID")
  
  #check this worked
  #Get taxonomy of nodes
  tax.tab = data.frame(tax_table(ps))
  
  #make ASV ID rownames a column
  tax.tab <- tibble::rownames_to_column(tax.tab, "ASVID")
  node_df <- tibble::rownames_to_column(node_df, "ASVID")
  
  #use dplyr to bind the dataframes
  res <- inner_join(tax.tab, node_df, by = "ASVID")
  
  #TEST
  t = sort(res$degree, decreasing=TRUE)
  t2 =sort(df5$Freq, decreasing=TRUE)
  t[1] == t2[1] #should equal TRUE
  
  #get the top 5 ASVs with the most correlations
  #order descending and get the first 5 rows
  top.5 = head(res[order(-res$degree),],5)
  
  #subset our degree weight dataframe by these top taxa
  df.cut <- df[df$weight <0,]
  names(df.cut) <- c("ASVID", "to", "weight")
  
  gen <- data.frame(cbind(net.tax.df$ASVID, net.tax.df$Genus))
  names(gen) <- c("ASVID", "Genus")
  df.gen <- merge(df.cut, gen, by="ASVID")
  names(df.gen) <- c("From", "ASVID", "weight", "FromGenus")
  df.gen2 <- merge(df.gen, gen, by="ASVID")
  names(df.gen2) <- c("From", "To", "weight", "FromGenus", "ToGenus")
  
  return(list=list(top.5,df.gen2))
}

sn.t.n <- top_degree_and_neg_corr(g.sn[[2]], ps, sn.node)
sn.t.n[[1]] #top 5 ASVs with high degree abundance
x1 = sn.t.n[[2]] #negative correlations
#get dataframe of number of negative correlations
x1.b = as.data.frame(sort(table(c(x1$FromGenus, x1$ToGenus)), decreasing = TRUE))
x1.b$Habitat = "Snow"

sp.t.n <- top_degree_and_neg_corr(g.sp[[2]], ps, sp.node)
sp.t.n[[1]] #top 5 ASVs with high degree abundance
x2 = sp.t.n[[2]] #negative correlations
#get dataframe of number of negative correlations
x2.b = as.data.frame(sort(table(c(x2$FromGenus, x2$ToGenus)), decreasing = TRUE))
x2.b$Habitat = "Spring Ice"

sm.t.n <- top_degree_and_neg_corr(g.sm[[2]], ps, sm.node)
sm.t.n[[1]] #top 5 ASVs with high degree abundance
x3 = sm.t.n[[2]] #negative correlations
#get dataframe of number of negative correlations
x3.b = as.data.frame(sort(table(c(x3$FromGenus, x3$ToGenus)), decreasing = TRUE))
x3.b$Habitat = "Summer Ice"

cr.t.n <- top_degree_and_neg_corr(g.cr[[2]], ps, cr.node)
cr.t.n[[1]] #top 5 ASVs with high degree abundance
x4 = cr.t.n[[2]] #negative correlations
#get dataframe of number of negative correlations
x4.b = as.data.frame(sort(table(c(x4$FromGenus, x4$ToGenus)), decreasing = TRUE))
x4.b$Habitat = "Cryoconite"

#combine and export
x = rbind(x1.b, x2.b, x3.b, x4.b)
write.csv(x, "../results/network-negative-associations.csv")

#Which pairs of ASVs commonly have negative interactions?
df.sn = data.frame(table(sn.t.n[[2]][,4:5]))
df.sn$Habitat = "Snow"
df.sp = data.frame(table(sp.t.n[[2]][,4:5]))
df.sp$Habitat = "Spring Ice"
df.sm = data.frame(table(sm.t.n[[2]][,4:5]))
df.sm$Habitat = "Summer Ice"
df.c = data.frame(table(cr.t.n[[2]][,4:5]))
df.c$Habitat = "Cryoconite"

#combine and export
x = rbind(df.sn, df.sp, df.sm, df.c)
x = x[x$Freq > 0,]

write.csv(x, "../results/network-negative-associations.csv")
