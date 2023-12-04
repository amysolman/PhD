#A) MANTEL CORRELOGRAM FUNCTION 

# phylo = ps
# 
# #make ps small to test function
# phylo2 <- subset_samples(phylo, Glacier=="Commonwealth")
# phylo2 = filter_taxa(phylo2, function(x) sum(x) >= 1, TRUE)
# phylo3= filter_taxa(phylo2, function(x) sum(x) < 30, TRUE)
# 
 #phylo = ps1

mantel_correlogram_func <- function(phylo){
  
  #PROPORTIONAL TRANSFORMATION
  ps.rel <- transform_sample_counts(phylo, function(x) x/sum(x))
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(ps.rel))
  
  # keeps = c("Latitude.z", "DistanceToSea", "Altitude", "WaterDepth", "TotalDepth", "SedimentDepth", "Area", "Ice.lid","Temp","EC","pH","pCO2", "DOC","DIC","C.N","Cl.age","DO","N","C","H","HCO3","NH4","NO2","NO3","TN","TIN","TON","DON","TDN","TDP","DOP","PO4","SiO2","Cl","SO4","Na","K","Mg","Ca","F","DRP")
  
  #list of variables we're interested in 
  keeps1 = c("Ca", "Mg", "K", "Na", "SO4", "Cl", "NO3", "EC", "pH", "Area", "TotalDepth")

  keeps2 = c("DistanceToSea", "Altitude", "HCO3")
  
  #samp_data.trim <- samp_data[ , (names(samp_data) %in% keeps)]
  
  samp_data1 <- samp_data[ , (names(samp_data) %in% keeps1)]
  samp_data2 <- samp_data[ , (names(samp_data) %in% keeps2)]
  
  #remove columns with more than 30% missing variables
  #samp_data.trim.trim <- samp_data.trim[ lapply(samp_data.trim, function(x) sum(is.na(x)) / length(x) ) < 0.1 ]
  samp_data_trim1 <- samp_data1[ lapply( samp_data1, function(x) sum(is.na(x)) / length(x) ) < 0.5 ]
  samp_data_trim2 <- samp_data2[ lapply( samp_data2, function(x) sum(is.na(x)) / length(x) ) < 0.7 ]
  
  #Only keep complete cases
  # meta <- samp_data.trim.trim[complete.cases(samp_data.trim.trim),]
  meta1 <- samp_data_trim1[complete.cases(samp_data_trim1),]
  meta2 <- samp_data_trim2[complete.cases(samp_data_trim2),]
  
  # #MANTEL CORRELOGRAM
  # ps.smol = prune_samples(rownames(meta1), ps.rel)
  # #sample_df1 <- data.frame(sample_data(trim1))
  # meta.trim <- meta1
  # 
  # #remove ASVs with 0 relative abundance
  # trim.filter <- filter_taxa(ps.smol, function(x) sum(x) > 0, TRUE)
  # 
  # #Get phylogenetic distances
  # tree <- phy_tree(trim.filter) #pull out phylogenetic tree
  # tree.dist <- cophenetic(tree)
  # asvs <- tree$tip.label
  # phylo.dists <- tree.dist[asvs, asvs]
  # phylo.dists[upper.tri(phylo.dists, diag=TRUE)] = NA
  # 
  # #Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
  # #turn the dataframe into a matrix
  # niches <- as.matrix(meta.trim) 
  # 
  # asv.table <- data.frame(t(otu_table(trim.filter)), check.names = FALSE)
  # asv.niche <- wascores(niches, asv.table) #this will find the weighted mean environmental parameter for each asv
  # asv.niche.df <- data.frame(asv.niche, check.names = FALSE)
  # 
  # #generate euclidean distance matrix for each ASV using combined environmental parameters
  # dist.out = as.matrix(dist(asv.niche.df), labels=TRUE)
  # 
  # #calculate mantel correlogram
  # corlg <- mantel.correlog(dist.out, phylo.dists, r.type="spearman")
  
  # #MANTEL CORRELOGRAM ONE
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
  # crlg <- data.frame(corlg$mantel.res)
  crlg1 <- data.frame(corlg1$mantel.res)
  crlg2 <- data.frame(corlg2$mantel.res)
  
  # crlg <- crlg %>%
  #   mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
  #   filter(!(is.na(Pr.corrected.)))
  
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


#B) Compute Beta-nearest Taxon Index (βNTI)

#B.1) Define function for calculating the βMNTD for each random null community

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


#B.2) Define main function for calculating BNTI

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

#C) Compute RC_Bray Model

#C.1) Define function for calculating null values

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


#C.2) Define main function for calculating RC Bray

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
