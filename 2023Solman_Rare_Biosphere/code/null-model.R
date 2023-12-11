# Method 

# The relative contribution of stochastic and deterministic processes in shaping cryoconite communities was quantified using a phylogenetic null model (Stegen et al., 2013). β mean nearest taxon distance metric (βMNTD) was calculated using the comdistnt function of the picante package in R (Kembel et al., 2010). For each pair of communities, the null βMNTD was calculated by randomly shuffling the tips of the phylogeny and recalculating the βMNTD 999 times to get distribution of null values.
# 
# Observed βMNTD and the null βMNTD were used to calculate β nearest taxon index (βNTI). Lower values of βNTI ( < -2) indicates ASVs are phylogenetically clustered due to homogenous selection. Higher values of βNTI ( > +2) indicates ASVs are variable selection. Where |βNTI| < +2 indicates stochastic processes are predominantly driving community assembly. 
# 
# As stochastic processes include homogenizing dispersal, disprsal limitation and ecological drift, the Bray–Curtis-based Raup–Crick metric (RCbray) is calculated. Observed Bray-Curtis values are compared to a null distribution to give the RCbray metric. Dispersal limitation is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and RCbray > +0.95. Homogenizing dispersal is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and RCbray < -0.95. Ecological drift is quantified as the proportion of pairwise comparisons with |βNTI| < +2 and |RCbray| < +0.95.
# 
# Phylogenetic null models rely on significant phylogenetic signal within the data. Taxa that are more closely related must exhibit niche similarities. To test this, Mantel correlation plots between Cophenetic distances and Eucldean distances between niche preferences were generated with 999 random permutations for test of significance. Cophenetic distances between each ASV were calculated. The abundance-weighted mean of each environmental variable for each ASV was calculated using function wascores in ‘vegan’ package (Oksanen et al., 2022). Euclidean distances between niche preferences for each ASV were calculated. Mantel correlation plots could not be produced for Arctic dataset as environmental data was not available. 
# 
# Variation in community assembly processes along each environmental gradient was assessed using regression analysis of BNTI values and Euclidean distances of major environmental variables. Mantel tests with 9999 permutations were used to test statistical significance. This method was also used to assess the relationship between phylogenetic turnover and environmental factors after controlling for geographic distance. This was carried out using the mantel function in the ecodist package. Code adapted from 
# #https://github.com/seb369/landuse_comm_assembly.  

rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan) #wascores
library(ggplot2)
library(picante)
library(ecodist) #for distance() function
library(parallel)
library(svglite)
library(dplyr) #for %>%
library(scales) #for percentages

# Step 2: Read in data
pro <- readRDS("../results/16S-phylo-object-rarefied-var-trans.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied-var-trans.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) > 0, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) > 0, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) > 0, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) > 0, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.int <- prune_samples(sample_sums(euk.ant.int)>0, euk.ant.int)
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")

pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds") 
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds") 
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds") 
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds")

#Make each dataset mini so we can locally test the whole script
# pro = rarefy_even_depth(pro, sample.size = min(sample_sums(pro))/100)
# euk = rarefy_even_depth(euk, sample.size = min(sample_sums(euk))/100)
# 
# pro.arc = rarefy_even_depth(pro.arc, sample.size = min(sample_sums(pro.arc))/100)
# pro.ant = rarefy_even_depth(pro.ant, sample.size = min(sample_sums(pro.ant))/1000)
# pro.ant <- filter_taxa(pro.ant, function(x) sum(x) > 0, TRUE)
# euk.arc = rarefy_even_depth(euk.arc, sample.size = min(sample_sums(euk.arc))/100)
# euk.ant = rarefy_even_depth(euk.ant, sample.size = min(sample_sums(euk.ant))/200)
# 
# pro.arc.abun = rarefy_even_depth(pro.arc.abun, sample.size = min(sample_sums(pro.arc.abun))/100)
# pro.arc.int = rarefy_even_depth(pro.arc.int, sample.size = min(sample_sums(pro.arc.int))/10)
# pro.arc.rare = rarefy_even_depth(pro.arc.rare, sample.size = min(sample_sums(pro.arc.rare))/10)
# pro.ant.abun = rarefy_even_depth(pro.ant.abun, sample.size = min(sample_sums(pro.ant.abun))/100)
# pro.ant.int = rarefy_even_depth(pro.ant.int, sample.size = min(sample_sums(pro.ant.int))/100)
# pro.ant.rare = rarefy_even_depth(pro.ant.rare, sample.size = min(sample_sums(pro.ant.rare))/5)
# 
# euk.arc.abun = rarefy_even_depth(euk.arc.abun, sample.size = min(sample_sums(euk.arc.abun))/100)
# euk.arc.int = rarefy_even_depth(euk.arc.int, sample.size = min(sample_sums(euk.arc.int))/10)
# euk.arc.rare = rarefy_even_depth(euk.arc.rare, sample.size = min(sample_sums(euk.arc.rare))/10)
# euk.ant.abun = rarefy_even_depth(euk.ant.abun, sample.size = min(sample_sums(euk.ant.abun))/100)
# euk.ant.int = rarefy_even_depth(euk.ant.int, sample.size = min(sample_sums(euk.ant.int))*10)
# euk.ant.rare = rarefy_even_depth(euk.ant.rare, sample.size = min(sample_sums(euk.ant.rare))*10)


#phylo = pro.ant

mantel_correlogram_func <- function(phylo){
  
  #PROPORTIONAL TRANSFORMATION
  ps.rel <- transform_sample_counts(phylo, function(x) x/sum(x))
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(ps.rel))
  
  #list of variables we're interested in 
  keeps1 <- c("WaterDepth", "SedimentDepth", "TotalDepth", "pH", "DOC", "Cl", "SO4", "Mg","Ca", "Area")
  
  keeps2 <-c("DistanceToSea", "Elevation", "HCO3")
  
  samp_data1 <- samp_data[ , (names(samp_data) %in% keeps1)]
  samp_data2 <- samp_data[ , (names(samp_data) %in% keeps2)]
  
  #remove columns with more than 50% missing variables
  samp_data_trim1 <- samp_data1[ lapply( samp_data1, function(x) sum(is.na(x)) / length(x) ) < 0.5 ]
  samp_data_trim2 <- samp_data2[ lapply( samp_data2, function(x) sum(is.na(x)) / length(x) ) < 0.7 ]
  
  #Only keep complete cases
  meta1 <- samp_data_trim1[complete.cases(samp_data_trim1),]
  meta2 <- samp_data_trim2[complete.cases(samp_data_trim2),]
  
  #MANTEL CORRELOGRAM ONE
  ps.smol1 = prune_samples(rownames(meta1), ps.rel)
  # sample_df1 <- data.frame(sample_data(trim1))
  meta.trim1 <- meta1
  
  #remove ASVs with 0 relative abundance
  trim.filter1 <- filter_taxa(ps.smol1, function(x) sum(x) > 0, TRUE)
  
  #Get phylogenetic distances
  tree1 <- phy_tree(trim.filter1) #pull out phylogenetic tree
  tree.dist1 <- cophenetic(tree1)
  asvs1 <- tree1$tip.label
  phylo.dists1 <- tree.dist1[asvs1, asvs1]
  phylo.dists1[upper.tri(phylo.dists1, diag=TRUE)] = NA
  
  #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  #turn the dataframe into a matrix
  niches1 <- as.matrix(meta.trim1) 
  
  asv.table1 <- data.frame(t(otu_table(trim.filter1)), check.names = FALSE)
  asv.niche1 <- wascores(niches1, asv.table1) #this will find the weighted mean environmental parameter for each asv
  asv.niche.df1 <- data.frame(asv.niche1, check.names = FALSE)
  
  #generate euclidean distance matrix for each ASV using combined environmental parameters
  dist.out1 = as.matrix(dist(asv.niche.df1), labels=TRUE)
  
  #calculate mantel correlogram
  corlg1 <- mantel.correlog(dist.out1, phylo.dists1, r.type="spearman")
  
  #MANTEL CORRELOGRAM TWO
  ps.smol2 = prune_samples(rownames(meta2), ps.rel)
  # sample_df1 <- data.frame(sample_data(trim1))
  meta.trim2 <- meta2
  
  #remove ASVs with 0 relative abundance
  trim.filter2 <- filter_taxa(ps.smol2, function(x) sum(x) > 0, TRUE)
  
  #Get phylogenetic distances
  tree2 <- phy_tree(trim.filter2) #pull out phylogenetic tree
  tree.dist2 <- cophenetic(tree2)
  asvs2 <- tree2$tip.label
  phylo.dists2 <- tree.dist2[asvs2, asvs2]
  phylo.dists2[upper.tri(phylo.dists2, diag=TRUE)] = NA
  
  #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  #turn the dataframe into a matrix
  niches2 <- as.matrix(meta.trim2) 
  
  asv.table2 <- data.frame(t(otu_table(trim.filter2)), check.names = FALSE)
  asv.niche2 <- wascores(niches2, asv.table2) #this will find the weighted mean environmental parameter for each asv
  asv.niche.df2 <- data.frame(asv.niche2, check.names = FALSE)
  
  #generate euclidean distance matrix for each ASV using combined environmental parameters
  dist.out2 = as.matrix(dist(asv.niche.df2), labels=TRUE)
  
  #calculate mantel correlogram
  corlg2 <- mantel.correlog(dist.out2, phylo.dists2, r.type="spearman")
  
  # Prep data
  crlg1 <- data.frame(corlg1$mantel.res)
  crlg2 <- data.frame(corlg2$mantel.res)
  
  crlg1 <- crlg1 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  
  crlg2 <- crlg2 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  
  crlg1$D.cl = rownames(crlg1)
  crlg1$Analysis = "1"
  crlg2$D.cl = rownames(crlg2)
  crlg2$Analysis = "2"
  
  final.crlg = rbind(crlg1, crlg2)
  
  
  return(final.crlg)
}


pro.ant.mc <- mantel_correlogram_func(pro.ant)
pro.ant.abun.mc <- mantel_correlogram_func(pro.ant.abun)
pro.ant.int.mc <- mantel_correlogram_func(pro.ant.int)
pro.ant.rare.mc <- mantel_correlogram_func(pro.ant.rare)

euk.ant.mc <- mantel_correlogram_func(euk.ant)
euk.ant.abun.mc <- mantel_correlogram_func(euk.ant.abun)
euk.ant.int.mc <- mantel_correlogram_func(euk.ant.int)
euk.ant.rare.mc <- mantel_correlogram_func(euk.ant.rare)

#combine into single dataframe
pro.ant.mc$Group = "Prokaryote"
pro.ant.mc$Subcommunity = "Full"
pro.ant.abun.mc$Group = "Prokaryote"
pro.ant.abun.mc$Subcommunity = "Abundant"
pro.ant.int.mc$Group = "Prokaryote"
pro.ant.int.mc$Subcommunity = "Intermediate"
pro.ant.rare.mc$Group = "Prokaryote"
pro.ant.rare.mc$Subcommunity = "Rare"
euk.ant.mc$Group = "Eukaryote"
euk.ant.mc$Subcommunity = "Full"
euk.ant.abun.mc$Group = "Eukaryote"
euk.ant.abun.mc$Subcommunity = "Abundant"
euk.ant.int.mc$Group = "Eukaryote"
euk.ant.int.mc$Subcommunity = "Intermediate"
euk.ant.rare.mc$Group = "Eukaryote"
euk.ant.rare.mc$Subcommunity = "Rare"

plot.df = rbind(pro.ant.mc, pro.ant.abun.mc, pro.ant.int.mc, pro.ant.rare.mc, euk.ant.mc, euk.ant.abun.mc, euk.ant.int.mc, euk.ant.rare.mc)

#save the results dataframe
write.csv(plot.df, "../results/mantel-correlogram-results.csv")

#order our groups
df2plot = plot.df[plot.df$Subcommunity != "Full",]

df2plot$Group = factor(df2plot$Group, levels=c("Prokaryote", "Eukaryote"))
df2plot$Subcommunity = factor(df2plot$Subcommunity, levels = c("Rare", "Intermediate", "Abundant"))

p = ggplot(data=df2plot, aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis)) +
  geom_point(data=df2plot[df2plot$sig=="significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black", size=3, shape=16)+
  geom_point(data=df2plot[df2plot$sig=="non-significant",], aes(x=class.index, y=Mantel.cor, group=Analysis, color=Analysis), color = "black",size=3, shape=1)+
  geom_line(size=1)+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x = "Phylogenetic distance class", y="Mantel correlation")+
 # ylim(-0.15, 0.1)+
  theme_bw()+
  theme(axis.text = element_text(size=7),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=10),
        legend.position = "bottom",
        legend.text = element_text(size=10))+
  facet_grid(rows=vars(Group), cols=vars(Subcommunity), scales = "free_x")

print(p)

pdf("../results/mantel-correlogram.pdf")
print(p)
dev.off()


#Compute Beta-nearest Taxon Index (βNTI)

#Define function for calculating the βMNTD for each random null community

# Function for calculating the βMNTD for each random null community
bMNTD_null_func <- function(i, OTU.table, tree){
  tree$tip.label = sample(tree$tip.label)
  bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  A <- attr(bMNTD_s, "Size")
  B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
  if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
  bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_s),
                          rep=i)
  return(bMNTD_s.df)
}


#Define main function for calculating BNTI

# The main function for calculating βNTI
Phylo_turnover <- function(physeq, reps, nproc){
  # Extract OTU table
  OTU.table = t(otu_table(physeq))
  # Extract phylogenetic tree
  tree = phy_tree(physeq)
  # Get βMNTD between all communities
  bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  A <- attr(bMNTD_o, "Size")
  B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
  if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
  bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_o))
  
  # Get βMNTD for randomized null communities
  rep.list = seq(1, reps)
  bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)
  
  # Combine all data together and calculate βNTI for each sample pair
  bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
  bMNTD_s.means.df = bMNTD_s.df %>%
    group_by(Sample_1, Sample_2) %>%
    dplyr::summarize(mean_bMNTD = mean(bMNTD),
                     sd_bMNTD = sd(bMNTD))
  
  bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
    mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)
  return(bMNTD_o.df)
}

# C) Compute RC_Bray Model

# 1 Define function for calculating null values

# Function for calculating the distances in the null communities
RCbray_null_func <- function(i, freq.abd.df, alpha1, alpha2, N){
  # Get simulated communities and distance
  ## initally select OTUs weighted by their frequency. The number of OTUs selected should equal the richness of the samples.
  simcom1 = data.frame(table(sample(freq.abd.df$OTU, size=alpha1, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
  colnames(simcom1) = c("OTU","simcom1")
  simcom1$OTU = as.character(simcom1$OTU)
  simcom1 = inner_join(simcom1, freq.abd.df, by="OTU")
  simcom2 = data.frame(table(sample(freq.abd.df$OTU, size=alpha2, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
  colnames(simcom2) = c("OTU","simcom2")
  simcom2$OTU = as.character(simcom2$OTU)
  simcom2 = inner_join(simcom2, freq.abd.df, by="OTU")
  
  ## Now recruit OTUs based on their abundance in the metacommunity
  simcom1.abd = data.frame(table(sample(simcom1$OTU, size=N-alpha1, replace=T, prob=simcom1$p)), stringsAsFactors = F)
  colnames(simcom1.abd) = c("OTU","simcom1.abd")
  simcom1.abd$OTU = as.character(simcom1.abd$OTU)
  simcom1 = full_join(simcom1, simcom1.abd, by="OTU") %>%
    mutate(simcom1.abd = ifelse(is.na(simcom1.abd), 1, simcom1.abd)) %>%
    select(OTU, simcom1.abd)
  
  simcom2.abd = data.frame(table(sample(simcom2$OTU, size=N-alpha2, replace=T, prob=simcom2$p)), stringsAsFactors = F)
  colnames(simcom2.abd) = c("OTU","simcom2.abd")
  simcom2.abd$OTU = as.character(simcom2.abd$OTU)
  simcom2 = full_join(simcom2, simcom2.abd, by="OTU") %>%
    mutate(simcom2.abd = ifelse(is.na(simcom2.abd), 1, simcom2.abd)) %>%
    select(OTU, simcom2.abd)
  
  
  simcom = full_join(simcom1, simcom2, by="OTU")
  simcom[is.na(simcom)] = 0
  rownames(simcom) = simcom$OTU
  simcom$OTU = NULL
  
  #relative abundance communities
  simcom.rel = simcom
  simcom.rel$simcom1.abd = simcom$simcom1.abd/sum(simcom$simcom1.abd)
  simcom.rel$simcom2.abd = simcom$simcom2.abd/sum(simcom$simcom2.abd)

  null.dist.rel = vegdist(t(simcom.rel), method="bray")[1]
  return(null.dist.rel)
}


#2 Define main function for calculating RC Bray

# Main function for calculating RCbray
Calc_RCbray <- function(physeq, reps, nproc){
  # Get OTU table from phyloseq object
  otu.table = otu_table(physeq)
  
  # Get alpha diversity for each sample
  otu.PA.table = otu.table
  otu.PA.table[otu.PA.table > 0] = 1
  alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)
  
  #transform data for Bray-Curtis calculation
  physeq.rel = transform_sample_counts(physeq, function(x) x/sum(x))
  
  #get community data table
  rel.tab = data.frame(t(otu_table(physeq.rel)))
  
  # Get beta diversity matrix
  beta.table = as.matrix(vegdist(rel.tab), method="bray", diag=TRUE, upper=TRUE)
  
  ## Get metacommunity
  # Calculate the number of individuals in the meta community (Average read depth)
  N <- mean(apply(t(otu.table), 1, sum))
  
  # Calculate the average relative abundance of each taxa across communities
  p.m <- apply(t(otu.table), 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  
  # Calculate the occurrence frequency of each taxa across communities
  otu.table.bi <- 1*(t(otu.table)>0)
  freq <- apply(otu.table.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # Combine
  freq.abd.df = data.frame(p=p, freq=freq) %>%
    tibble::rownames_to_column(var="OTU") %>%
    filter(p != 0, freq != 0) %>%
    arrange(p)
  
  # For each pair of samples run the RCbray analysis
  comps = combn(alpha.df$Sample_ID, m=2, simplify = F)
  RCb.df = data.frame(Site1 = character(), Site2 = character(), RCb = numeric(), stringsAsFactors = F)
  for (j in seq(1, length(comps))){
    sam = comps[[j]]
    alpha1 = alpha.df[alpha.df$Sample_ID == sam[1],]$OTU.n
    alpha2 = alpha.df[alpha.df$Sample_ID == sam[2],]$OTU.n
    # Permute "reps" many times
    rep.list = seq(1, reps)
    null.list = mclapply(rep.list, RCbray_null_func, freq.abd.df=freq.abd.df, alpha1=alpha1, alpha2=alpha2, N=N, mc.cores=nproc)
    
    RCb = (length(null.list[null.list > beta.table[sam[1], sam[2]]]) + (0.5*length(null.list[null.list == beta.table[sam[1], sam[2]]])))/reps
    RCb = (RCb - 0.5)*2
    
    RCb.df = rbind(RCb.df, data.frame(Site1=sam[1], Site2=sam[2], RCb=RCb, stringsAsFactors = F))
  }
  
  RCb.df
  return(RCb.df)
}

##RUN THE MODEL

get_bnti_res <- function(phylo, subcom, group, pole){
  
  
  df = Phylo_turnover(transform_sample_counts(phylo, function(x) x/sum(x)), 1000, 10) #should be 1000 10
  df$Group = group
  df$Subcommunity = subcom
  df$Pole = pole
  return(df)
}

pro.ant.bNTI.df = get_bnti_res(pro.ant, "Full", "Prokaryote", "Antarctic")
pro.ant.abun.bNTI.df = get_bnti_res(pro.ant.abun, "Abundant", "Prokaryote", "Antarctic")
pro.ant.int.bNTI.df = get_bnti_res(pro.ant.int, "Intermediate", "Prokaryote", "Antarctic")
pro.ant.rare.bNTI.df = get_bnti_res(pro.ant.rare, "Rare", "Prokaryote", "Antarctic")
pro.arc.bNTI.df = get_bnti_res(pro.arc, "Full", "Prokaryote", "Arctic")
pro.arc.abun.bNTI.df = get_bnti_res(pro.arc.abun, "Abundant", "Prokaryote", "Arctic")
pro.arc.int.bNTI.df = get_bnti_res(pro.arc.int, "Intermediate", "Prokaryote", "Arctic")
pro.arc.rare.bNTI.df = get_bnti_res(pro.arc.rare, "Rare", "Prokaryote", "Arctic")
euk.ant.bNTI.df = get_bnti_res(euk.ant, "Full", "Eukaryote", "Antarctic")
euk.ant.abun.bNTI.df = get_bnti_res(euk.ant.abun, "Abundant", "Eukaryote", "Antarctic")
euk.ant.int.bNTI.df = get_bnti_res(euk.ant.int, "Intermediate", "Eukaryote", "Antarctic")
euk.ant.rare.bNTI.df = get_bnti_res(euk.ant.rare, "Rare", "Eukaryote", "Antarctic")
euk.arc.bNTI.df = get_bnti_res(euk.arc, "Full", "Eukaryote", "Arctic")
euk.arc.abun.bNTI.df = get_bnti_res(euk.arc.abun, "Abundant", "Eukaryote", "Arctic")
euk.arc.int.bNTI.df = get_bnti_res(euk.arc.int, "Intermediate", "Eukaryote", "Arctic")
euk.arc.rare.bNTI.df = get_bnti_res(euk.arc.rare, "Rare", "Eukaryote", "Arctic")

full.bNTI.df = rbind(pro.ant.bNTI.df, pro.ant.abun.bNTI.df, pro.ant.int.bNTI.df, pro.ant.rare.bNTI.df, 
                     pro.arc.bNTI.df, pro.arc.abun.bNTI.df, pro.arc.int.bNTI.df, pro.arc.rare.bNTI.df,
                     euk.ant.bNTI.df, euk.ant.abun.bNTI.df, euk.ant.int.bNTI.df, euk.ant.rare.bNTI.df, 
                     euk.arc.bNTI.df,euk.arc.abun.bNTI.df, euk.arc.int.bNTI.df, euk.arc.rare.bNTI.df)

write.csv(full.bNTI.df, "../results/bNTI-results-table.csv")

get_rcbray_res <- function(phylo, subcom, group, pole){
  
  
  df = Calc_RCbray(phylo, 999, 20) #should be 999 20
  df$Group = group
  df$Subcommunity = subcom
  df$Pole = pole
  return(df)
  
}

pro.ant.rcbray.df = get_rcbray_res(pro.ant, "Full", "Prokaryote", "Antarctic")
pro.ant.abun.rcbray.df = get_rcbray_res(pro.ant.abun, "Abundant", "Prokaryote", "Antarctic")
pro.ant.int.rcbray.df = get_rcbray_res(pro.ant.int, "Intermediate", "Prokaryote", "Antarctic")
pro.ant.rare.rcbray.df = get_rcbray_res(pro.ant.rare, "Rare", "Prokaryote", "Antarctic")
pro.arc.rcbray.df = get_rcbray_res(pro.arc, "Full", "Prokaryote", "Arctic")
pro.arc.abun.rcbray.df = get_rcbray_res(pro.arc.abun, "Abundant", "Prokaryote", "Arctic")
pro.arc.int.rcbray.df = get_rcbray_res(pro.arc.int, "Intermediate", "Prokaryote", "Arctic")
pro.arc.rare.rcbray.df = get_rcbray_res(pro.arc.rare, "Rare", "Prokaryote", "Arctic")
euk.ant.rcbray.df = get_rcbray_res(euk.ant, "Full", "Eukaryote", "Antarctic")
euk.ant.abun.rcbray.df = get_rcbray_res(euk.ant.abun, "Abundant", "Eukaryote", "Antarctic")
euk.ant.int.rcbray.df = get_rcbray_res(euk.ant.int, "Intermediate", "Eukaryote", "Antarctic")
euk.ant.rare.rcbray.df = get_rcbray_res(euk.ant.rare, "Rare", "Eukaryote", "Antarctic")
euk.arc.rcbray.df = get_rcbray_res(euk.arc, "Full", "Eukaryote", "Arctic")
euk.arc.abun.rcbray.df = get_rcbray_res(euk.arc.abun, "Abundant", "Eukaryote", "Arctic")
euk.arc.int.rcbray.df = get_rcbray_res(euk.arc.int, "Intermediate", "Eukaryote", "Arctic")
euk.arc.rare.rcbray.df = get_rcbray_res(euk.arc.rare, "Rare", "Eukaryote", "Arctic")

full.rcbray.df = rbind(pro.ant.rcbray.df, pro.ant.abun.rcbray.df, pro.ant.int.rcbray.df, pro.ant.rare.rcbray.df, 
                     pro.arc.rcbray.df, pro.arc.abun.rcbray.df, pro.arc.int.rcbray.df, pro.arc.rare.rcbray.df,
                     euk.ant.rcbray.df, euk.ant.abun.rcbray.df, euk.ant.int.rcbray.df, euk.ant.rare.rcbray.df, 
                     euk.arc.rcbray.df,euk.arc.abun.rcbray.df, euk.arc.int.rcbray.df, euk.arc.rare.rcbray.df)

write.csv(full.rcbray.df, "../results/RCbray-results-table.csv")


#GET RESULTS READY FOR PLOTTING

#create dataframes
#make sure everything is formatted correctly for merging
RC_bray <- full.rcbray.df
colnames(RC_bray) <- c("Sample_2", "Sample_1", "RCb", "Group", "Subcommunity", "Pole")
RC_bray$Sample_1 <- as.integer(RC_bray$Sample_1)
RC_bray$Sample_2 <- as.integer(RC_bray$Sample_2)
full.bNTI.df$Sample_1 <- as.integer(full.bNTI.df$Sample_1)
full.bNTI.df$Sample_2 <- as.integer(full.bNTI.df$Sample_2)

turnover.df = inner_join(full.bNTI.df, RC_bray)

turnover.df = turnover.df %>%
  mutate(process = ifelse(abs(bNTI) < 2, #if absolute value of BNTI is less than 2
                          ifelse(abs(RCb) < 0.95, "Drift", #and absolute value of RCb is less than 0.95 = Drift
                                 ifelse(RCb >= 0.95, "Dispersal Limited", #if RCb is > 0.95 = Dispersal limited
                                        ifelse(RCb <= -0.95, "Homogenizing Dispersal", "ERROR"))), #< -0.95 = homo disp
                          ifelse(bNTI >= 2, "Variable Selection", #more than 2 = variable selection
                                 ifelse(bNTI <= -2, "Homogeneous Selection", "ERROR")))) #< -2 = homogenous selection

turnover.df$process = factor(turnover.df$process, levels = c("Drift", 
                                                                     "Dispersal Limited", "Homogenizing Dispersal", 
                                                                     "Variable Selection", "Homogeneous Selection"))

#for each abundance class
#find the number of site pairs
#and the percentage of those pairs that shows that process

abundance.list <- unique(turnover.df$Subcommunity)

#Arctic Prokaryotes
arc_pro_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Abundance=character())
new.df <- subset(turnover.df, Pole == "Arctic" & Group == "Prokaryote")

for (i in 1:length(abundance.list)){
  x <- new.df[new.df$Subcommunity == abundance.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Subcommunity<- abundance.list[i]
  arc_pro_final_data <- rbind(arc_pro_final_data, x)
}
arc_pro_final_data$Pole = "Arctic"
arc_pro_final_data$Group = "Prokaryote"

#Arctic Eukaryotes
arc_euk_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Abundance=character())
new.df <- subset(turnover.df, Pole == "Arctic" & Group == "Eukaryote")

for (i in 1:length(abundance.list)){
  x <- new.df[new.df$Subcommunity == abundance.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Subcommunity<- abundance.list[i]
  arc_euk_final_data <- rbind(arc_euk_final_data, x)
}
arc_euk_final_data$Pole = "Arctic"
arc_euk_final_data$Group = "Eukaryote"

#Antarctic Prokaryotes
ant_pro_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Abundance=character())
new.df <- subset(turnover.df, Pole == "Antarctic" & Group == "Prokaryote")

for (i in 1:length(abundance.list)){
  x <- new.df[new.df$Subcommunity == abundance.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Subcommunity<- abundance.list[i]
  ant_pro_final_data <- rbind(ant_pro_final_data, x)
}
ant_pro_final_data$Pole = "Antarctic"
ant_pro_final_data$Group = "Prokaryote"

#Antarctic eukaryotes
ant_euk_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Abundance=character())
new.df <- subset(turnover.df, Pole == "Antarctic" & Group == "Eukaryote")

for (i in 1:length(abundance.list)){
  x <- new.df[new.df$Subcommunity == abundance.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Subcommunity<- abundance.list[i]
  ant_euk_final_data <- rbind(ant_euk_final_data, x)
}
ant_euk_final_data$Pole = "Antarctic"
ant_euk_final_data$Group = "Eukaryote"

#bind results together
final_data = rbind(arc_pro_final_data, arc_euk_final_data, ant_pro_final_data, ant_euk_final_data)

#add percentages column
final_data.edit <- final_data %>% dplyr::group_by(Subcommunity, Pole, Group) %>%
  dplyr::mutate(frac = n_sites / sum(n_sites))
final_data.edit$frac = round(final_data.edit$frac,2)

#PLOT RESULTS
#final_data$round_perc <- round(final_data$perc, 2)

#remove rows with less than 1%
#x.cut <- x[which(x$frac > 0.01),]

#remove the full communitiy from plotting
final_data.edit.plot = subset(final_data.edit, Subcommunity != "Full")

#make group a factor
final_data.edit.plot$Group = factor(final_data.edit.plot$Group, levels=c("Prokaryote", "Eukaryote"))

  p1 = ggplot(final_data.edit.plot, aes(x = Subcommunity,y = n_sites, 
                              fill = factor(process, levels=c("Homogeneous Selection", "Variable Selection", "Homogenizing Dispersal", "Dispersal Limited", "Drift")))) + 
  geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels=scales::percent)+
          geom_text(
    aes(label = percent(frac)), position = position_fill(0.5), size=5) +
  facet_grid(~ Group + Pole)+
    theme_bw()+
    ylab("Relative Contribution %")+
    theme(text=element_text(size=12), 
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.title.x=element_blank(), 
          axis.text=element_text(size=12), 
          legend.text = element_text(size=15))+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6", "#6ca8ab", "#517e80"))
  
   print(p1)
   
  pdf("../results/null-model.pdf", width=14, height=8)
  print(p1)
  dev.off()