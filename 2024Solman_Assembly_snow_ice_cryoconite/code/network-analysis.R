#Script breakdown
#1. Clear workspace and load packages
#2. Load data
#3. Get abundance dataframes
#4. Merge 16S and 18S datasets
#5. Filter ASVs by percentage of samples
#6. Carry out correlation analysis


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

#get community data
#prokaryotes
ps.pro <- readRDS("../results/16S-ps-controls-removed.rds")
#eukaryotes without microfauna
ps.euk <- readRDS("../results/18S-ps-no-mm-controls-removed.rds")
#microfauna only
ps.mm <- readRDS("../results/18S-ps-mm-only-controls-removed.rds")

#remove the sample low read count samples as we did previously

#prokaryotes
ps.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro) #retain samples with >= num counts
ps.pro = filter_taxa(ps.pro, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts
#eukaryotes
ps.euk = prune_samples(sample_sums(ps.euk)>= 1900, ps.euk) #retain samples with >= num counts
ps.euk = filter_taxa(ps.euk, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts
#micrometazoans
ps.mm = prune_samples(sample_sums(ps.mm)>= 500, ps.mm) #retain samples with >= num counts
ps.mm = filter_taxa(ps.mm, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts

#3. Merge the datasets

merge_ps <- function(ps1, ps2, ps3){
  
  #Extract ASV tables
  tab1 = data.frame(t(otu_table(ps1)), check.names = FALSE)
  tab2 = data.frame(t(otu_table(ps2)), check.names = FALSE)
  tab3 = data.frame(t(otu_table(ps3)), check.names = FALSE)
  
  #sample name to column
  tab1 <- tibble::rownames_to_column(tab1, "SampleID")
  tab2 <- tibble::rownames_to_column(tab2, "SampleID")
  tab3 <- tibble::rownames_to_column(tab3, "SampleID")
  
  #use dplyr to join the ASV tables
  full.tab <- full_join(tab1, tab2, by="SampleID")
  full.tab[is.na(full.tab)] <- 0
  full.tab <- full_join(full.tab, tab3, by="SampleID")
  full.tab[is.na(full.tab)] <- 0
  
  #test if this worked
  ncol(tab1) + ncol(tab2) + ncol(tab3) - 3 == ncol(full.tab) - 1 #should say TRUE
  
  #get rid of our sample ID column
  full.tab <- tibble::column_to_rownames(full.tab, var="SampleID")
  
  #combine taxa tables
  full.tax = as.matrix(rbind(data.frame(tax_table(ps1)), data.frame(tax_table(ps2)), data.frame(tax_table(ps3))))
  
  #combine metadata
  meta1 = data.frame(sample_data(ps1))
  meta2 = data.frame(sample_data(ps2))
  meta3 = data.frame(sample_data(ps3))
  
  #get vector of all our samples
  all.samp = unique(c(meta1$SampleID, meta2$SampleID, meta3$SampleID))
  
  #bind the rows of our metadata
  full.meta = rbind(meta1[,-c(2)], meta2[,-c(2)], meta3[,-c(2)])
  
  #remove rownames
  rownames(full.meta) <- NULL
  
  #remove duplicate rows using dplyr
  full.meta2 <- full.meta %>% distinct()
  
  #check that this has worked
  all.samp == full.meta2$SampleID #should all equal TRUE
  
  #make rownames sample ID again
  rownames(full.meta2) <- full.meta2$SampleID
  
  #make new phyloseq object
  ASV = otu_table(t(full.tab), taxa_are_rows = TRUE)
  TAX = tax_table(full.tax)
  META = sample_data(full.meta2)
  new.phylo = phyloseq(ASV, TAX, META)
  
  return(new.phylo)
  
}

ps.f <- merge_ps(ps.pro, ps.euk, ps.mm)

#export combine phyloseq object
saveRDS(ps.f, "../results/combined-phylo-for-network-analysis.rds")

#4. Carry out correlation analysis
# ps = ps.f
# hab = data.frame(sample_data(subset_samples(ps.f, Habitat == "Snow")))$SampleID
# perc = 0.2
# habitat = "Snow"

get_my_network <- function(ps, hab, habitat, perc){
  
  set.seed(666) 
  
  #Separate habitat data
  h <- prune_samples(hab, ps)
  
  #filter by percentage abundance
  ps.ab.filt = h %>%
    tax_filter(min_prevalence=perc)
  
  #convert to relative abundance data
  ps.rel <- transform_sample_counts(ps.ab.filt, function(x) x/sum(x))
  
  #get our asv table
  asv.tab <- data.frame(t(otu_table(ps.rel)), check.names=FALSE)
  #check this has worked
  rowSums(asv.tab) #should all be 1
  
  #find the spearmans correlations between asvs
  cor_analysis <-rcorr(as.matrix(asv.tab),type="spearman")
  
  #assign the matrix of correlations and the p-values to r and p, respectively
  cor_r<-cor_analysis$r
  cor_p<-cor_analysis$P 
  
  #matrix diagonals - apply the value 1 to the diagonal values of the p matrix (instead of NA as they were before)
  diag(cor_p)<-1
  
  #adjust p-values using the false discovery rate method
  cor_p.adj <- p.adjust(cor_p, method="fdr")
  
  #how many of our adjusted p values are < 0.01
  length(cor_p.adj[cor_p.adj < 0.01]) 
  
  #make any value <= 0.01 = 1 and make any value of p > 0.01 = 0
  p.adj <- as.numeric(cor_p.adj <= 0.01)
  
  #let's look at the range of correlation coefficients
  hist(cor_r)
  
  #are there any with coefficients > | 0.8|?
  length(cor_r[cor_r < -0.8]) #there are  negative correlation coefficients
  length(cor_r[cor_r > 0.8]) #there  positive correlation coefficients
  
  #change the value of any correlation coefficients <  0.8  to 0
  cor_r_pos <- cor_r
  cor_r_pos[cor_r_pos < 0.8] <- 0
  cor_r_all <- cor_r
  cor_r_all[abs(cor_r_all) < 0.8] <- 0
  
  
  #multiple the correlation coefficients by adjusted p values so that insignificant coefficients become 0
  cor_all<-cor_r_all*p.adj
  cor_pos <- cor_r_pos*p.adj
  
  #are there any with coefficients > | 0.8| after removing non-significant correlations?
  neg <- length(cor_all[cor_all < -0.8]) #there are  negative correlation coefficients
  pos <- length(cor_all[cor_all > 0.8]) #there  positive correlation coefficients
  length(cor_all[cor_all > 0.8]) == length(cor_pos[cor_pos > 0.8]) #should be TRUE
  
  #create igraph graph from adjacency matrix
  cor_g_all <- graph.adjacency(cor_all, weighted=TRUE, mode="undirected")
  cor_g_pos <- graph.adjacency(cor_pos, weighted=TRUE, mode="undirected")
  
  #simplify the graph so it does not contain loops or multiple edges
  cor_g_all <- simplify(cor_g_all)
  cor_g_pos <- simplify(cor_g_pos)
  
  #delete vertices from the graph
  cor_g_all<-delete.vertices(cor_g_all,names(degree(cor_g_all)[degree(cor_g_all)==0]))
  cor_g_pos<-delete.vertices(cor_g_pos,names(degree(cor_g_pos)[degree(cor_g_pos)==0]))
  
  gorder(cor_g_all) #number of nodes
  gsize(cor_g_all) #number of edges
  gorder(cor_g_pos) #number of nodes
  gsize(cor_g_pos) #number of edges
  x = igraph::as_data_frame(cor_g_all)
  length(x$weight[x$weight < 0]) == gsize(cor_g_all) - gsize(cor_g_pos) #Should equal TRUE
  
  #export graph
  #positive and negative correlations
  write.graph(cor_g_all, paste0("../results/positive-and-negative-", habitat, "-network.gml"), format="gml")
  #positive correlations only
  write.graph(cor_g_pos, paste0("../results/positive-", habitat, "-network.gml"), format="gml")
  
  #transform negative values into positive values in full network so we can calculate our topographic properties
  cor_all_trans <- igraph::as_data_frame(cor_g_all) #extract network data
  cor_all_trans$weight = abs(cor_all_trans$weight) #convert to absolute value
  g.new = graph_from_data_frame(cor_all_trans, directed=FALSE) #make back into graph
  gorder(g.new) == gorder(cor_g_all) #should equal TRUE
  gsize(g.new) == gsize(cor_g_all) #should equal TRUE
  
  #export g.new so we can calculate topological properties in the next analysis scripts
  write.graph(g.new, paste0("../results/full-", habitat, "-network-for-calc-properties.gml"), format="gml")
  
  
  res.list = list(g.new, cor_g_all, cor_g_pos)
  
  return(res.list)
}

g.sn = get_my_network(ps = ps.f, 
                      hab = data.frame(sample_data(subset_samples(ps.f, Habitat == "Snow")))$SampleID, 
                      perc = 0.2,
                      habitat = "Snow")
g.sp = get_my_network(ps = ps.f, 
                      hab = data.frame(sample_data(subset_samples(ps.f, Habitat == "Spring Ice")))$SampleID, 
                      perc = 0.2,
                      habitat = "SpringIce")
g.sm = get_my_network(ps = ps.f, 
                      hab = data.frame(sample_data(subset_samples(ps.f, Habitat == "Summer Ice")))$SampleID, 
                      perc = 0.2,
                      habitat = "SummerIce")
g.cr = get_my_network(ps = ps.f, 
                      hab = data.frame(sample_data(subset_samples(ps.f, Habitat == "Cryoconite")))$SampleID, 
                      perc = 0.2,
                      habitat = "Cryoconite")
