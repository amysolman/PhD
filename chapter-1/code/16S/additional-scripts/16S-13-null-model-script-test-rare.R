
rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan) #wascores
library(dplyr) #for %>%
library(ggplot2)
library(picante)
#detach("package:phyloseq", unload=TRUE)
library(ecodist) #for distance() function
#library(tibble) #rownames_to_column
library(parallel)

#TOTAL DATASET
#ps <- readRDS("../../results/16S/phylo-objects/16S-phyloseq-object-rarefied-decontam.rds")

# #Total rare
#ps.abun <- readRDS("../../results/16S/phylo-objects/16S-total-abundant.rds")
# #Total rare
# ps.int <- readRDS("../../results/16S/phylo-objects/16S-total-intermediate.rds")
# #Total Rare
ps.rare <- readRDS("../../results/16S/phylo-objects/16S-total-rare.rds")


set.seed(6666)
#I know my data is even sampling depth so I'm just going to get relative abundances
#get relative abundances
#rarefy to small sample depth so I can run this code on my computer
new.ps <- rarefy_even_depth(ps.rare, sample.size = 2)

rel_abund <- transform_sample_counts(new.ps, function(x) x/sum(x))


meta_df <- function(phylo){
  
  #Get metadata of samples
  samp_data <- data.frame(sample_data(phylo))
  
  #list of variables we're interested in 
  keeps <- c("Distance_To_Sea", "Elevation", "Water_Depth", "Sediment_Depth", "Total_Depth","Conductivity", "pH", "DOC_mgL.1",
             "Cl_merge", "SO4_merge", "Na_merge", "K_merge", "Mg_merge","Ca_merge", "HCO3_merge", "Radius", "EW", "NS")
  samp_data <- samp_data[ , (names(samp_data) %in% keeps)]
  
  #get cryoconite hole areas
  area1 <- pi*samp_data$Radius^2
  area2 <- pi*(samp_data$NS/2)*(samp_data$EW/2)
  samp_data$Area <- coalesce(area1,area2)
  
  drops <- c("EW", "NS", "Radius")
  samp_data <- samp_data[ , !(names(samp_data) %in% drops)]
  
  #make sure data are numeric
  samp_data[1:16] <- data.frame(lapply(samp_data[1:16],as.numeric))
  
  #change sample names 
  names(samp_data) <- c("Distance_To_Sea", "Elevation", "Water_Depth", "Sediment_Depth", "Total_Depth","Conductivity", "pH", "DOC_mgL.1", "Cl", "SO4", "Na", "K", "Mg","Ca", "HCO3", "Area")
  
  #remove columns with more than 50% missing variables
  samp_data_trim <- samp_data[ lapply( samp_data, function(x) sum(is.na(x)) / length(x) ) < 0.5 ]
  
  #Only keep complete cases
  samp_data_trim_complete <- samp_data_trim[complete.cases(samp_data_trim),]
  
  return(samp_data_trim_complete)
}

#get meta data
meta <- meta_df(rel_abund)

#make sure metadata and phyloseq object have the sample samples
samples_to_keep <- as.numeric(rownames(meta))
rel.abun.trim = subset_samples(rel_abund, SampleID %in% samples_to_keep)
sample_df <- data.frame(sample_data(rel.abun.trim))
meta.trim <- subset(meta, rownames(meta) %in% sample_df$SampleID)

#remove ASVs with 0 relative abundance
rel.abun.trim.filter <- filter_taxa(rel.abun.trim, function(x) sum(x) > 0, TRUE)

#Get phylogenetic distances
tree <- phy_tree(rel.abun.trim.filter) #pull out phylogenetic tree
tree.dist <- cophenetic(tree)
asvs <- tree$tip.label
phylo.dists <- tree.dist[asvs, asvs]
phylo.dists[upper.tri(phylo.dists, diag=TRUE)] = NA

#give meta sample names as row names
#rownames(meta) <- rownames(metadata) 

#turn the dataframe into a matrix
niches <- as.matrix(meta.trim) 

asv.table <- data.frame(t(otu_table(rel.abun.trim.filter)), check.names = FALSE)
asv.niche <- wascores(niches, asv.table) #this will find the weighted mean environmental parameter for each asv
asv.niche.df <- data.frame(asv.niche, check.names = FALSE)


#Conductivity
Cond <- asv.niche.df$Conductivity #get Dist2Sea niche preferences
names(Cond) <- rownames(asv.niche.df) #match them with the asv rowname
Cond.dist <- as.matrix(dist(Cond), labels=TRUE) #convert into a matrix
x <- Cond.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Cond.corlg <- mantel.correlog(x, phylo.dists)

# #Elevation
# Elevation <- asv.niche.df$Elevation #get niche preferences
# names(Elevation) <- rownames(asv.niche.df) #match them with the asv rowname
# Elevation.dist <- as.matrix(dist(Elevation), labels=TRUE) #convert into a matrix
# x <- Elevation.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# Elevation.corlg <- mantel.correlog(x, phylo.dists)

# #Total_Depth
# Total_Depth <- asv.niche.df$Total_Depth #get niche preferences
# names(Total_Depth) <- rownames(asv.niche.df) #match them with the asv rowname
# Total_Depth.dist <- as.matrix(dist(Total_Depth), labels=TRUE) #convert into a matrix
# x <- Total_Depth.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# Total_Depth.corlg <- mantel.correlog(x, phylo.dists)

#pH
pH <- asv.niche.df$pH #get niche preferences
names(pH) <- rownames(asv.niche.df) #match them with the asv rowname
pH.dist <- as.matrix(dist(pH), labels=TRUE) #convert into a matrix
x <- pH.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
pH.corlg <- mantel.correlog(x, phylo.dists)

#Cl
Cl <- asv.niche.df$Cl #get niche preferences
names(Cl) <- rownames(asv.niche.df) #match them with the asv rowname
Cl.dist <- as.matrix(dist(Cl), labels=TRUE) #convert into a matrix
x <- Cl.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Cl.corlg <- mantel.correlog(x, phylo.dists)

#SO4
SO4 <- asv.niche.df$SO4 #get niche preferences
names(SO4) <- rownames(asv.niche.df) #match them with the asv rowname
SO4.dist <- as.matrix(dist(SO4), labels=TRUE) #convert into a matrix
x <- SO4.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
SO4.corlg <- mantel.correlog(x, phylo.dists)

#Na
Na <- asv.niche.df$Na #get niche preferences
names(Na) <- rownames(asv.niche.df) #match them with the asv rowname
Na.dist <- as.matrix(dist(Na), labels=TRUE) #convert into a matrix
x <- Na.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Na.corlg <- mantel.correlog(x, phylo.dists)

#K
K <- asv.niche.df$K #get niche preferences
names(K) <- rownames(asv.niche.df) #match them with the asv rowname
K.dist <- as.matrix(dist(K), labels=TRUE) #convert into a matrix
x <- K.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
K.corlg <- mantel.correlog(x, phylo.dists)

#Mg
Mg <- asv.niche.df$Mg #get niche preferences
names(Mg) <- rownames(asv.niche.df) #match them with the asv rowname
Mg.dist <- as.matrix(dist(Mg), labels=TRUE) #convert into a matrix
x <- Mg.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Mg.corlg <- mantel.correlog(x, phylo.dists)

#Ca
Ca <- asv.niche.df$Ca #get niche preferences
names(Ca) <- rownames(asv.niche.df) #match them with the asv rowname
Ca.dist <- as.matrix(dist(Ca), labels=TRUE) #convert into a matrix
x <- Ca.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Ca.corlg <- mantel.correlog(x, phylo.dists)

# #DOC
# DOC <- asv.niche.df$DOC #get niche preferences
# names(DOC) <- rownames(asv.niche.df) #match them with the asv rowname
# DOC.dist <- as.matrix(dist(DOC), labels=TRUE) #convert into a matrix
# x <- DOC.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# DOC.corlg <- mantel.correlog(x, phylo.dists)
# 
# #TDN
# TDN <- asv.niche.df$TDN #get niche preferences
# names(TDN) <- rownames(asv.niche.df) #match them with the asv rowname
# TDN.dist <- as.matrix(dist(TDN), labels=TRUE) #convert into a matrix
# x <- TDN.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# TDN.corlg <- mantel.correlog(x, phylo.dists)
# 
# #TDP
# TDP <- asv.niche.df$TDP #get niche preferences
# names(TDP) <- rownames(asv.niche.df) #match them with the asv rowname
# TDP.dist <- as.matrix(dist(TDP), labels=TRUE) #convert into a matrix
# x <- TDP.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# TDP.corlg <- mantel.correlog(x, phylo.dists)
# 
# #C.N
# C.N <- asv.niche.df$C.N #get niche preferences
# names(C.N) <- rownames(asv.niche.df) #match them with the asv rowname
# C.N.dist <- as.matrix(dist(C.N), labels=TRUE) #convert into a matrix
# x <- C.N.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# C.N.corlg <- mantel.correlog(x, phylo.dists)
# 
# #NO3.N_mueq
# NO3.N_mueq <- asv.niche.df$NO3.N_mueq #get niche preferences
# names(NO3.N_mueq) <- rownames(asv.niche.df) #match them with the asv rowname
# NO3.N_mueq.dist <- as.matrix(dist(NO3.N_mueq), labels=TRUE) #convert into a matrix
# x <- NO3.N_mueq.dist[asvs, asvs]
# x[upper.tri(x, diag = TRUE)] = NA
# NO3.N_mueq.corlg <- mantel.correlog(x, phylo.dists)


# Prep data
# Dist2Sea.crlg <- data.frame(Dist2Sea.corlg$mantel.res) %>%
#   mutate(property = "Dist2Sea")
# Elevation.crlg <- data.frame(Elevation.corlg$mantel.res) %>%
#   mutate(property = "Elevation")
# Total_Depth.crlg <- data.frame(Total_Depth.corlg$mantel.res) %>%
#   mutate(property = "Total_Depth")
Cond.crlg <- data.frame(Cond.corlg$mantel.res) %>%
  mutate(property = "Cond")
pH.crlg <- data.frame(pH.corlg$mantel.res) %>%
  mutate(property = "pH")
Cl.crlg <- data.frame(Cl.corlg$mantel.res) %>%
  mutate(property = "Cl")
SO4.crlg <- data.frame(SO4.corlg$mantel.res) %>%
  mutate(property = "SO4")
Na.crlg <- data.frame(Na.corlg$mantel.res) %>%
  mutate(property = "Na")
K.crlg <- data.frame(K.corlg$mantel.res) %>%
  mutate(property = "K")
Mg.crlg <- data.frame(Mg.corlg$mantel.res) %>%
  mutate(property = "Mg")
Ca.crlg <- data.frame(Ca.corlg$mantel.res) %>%
  mutate(property = "Ca")

crlg <- rbind(Cond.crlg, pH.crlg, Cl.crlg, SO4.crlg, Na.crlg, K.crlg, Mg.crlg, Ca.crlg) %>%
  mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
  filter(!(is.na(Pr.corrected.)))

#Save crlg table
write.table(crlg, "rare-crlg.txt")

p <- ggplot(data=crlg, aes(x=class.index, y=Mantel.cor)) +
  geom_point(data=crlg[crlg$sig=="significant",], color = "black", size=2, shape=16) +
  geom_point(data=crlg[crlg$sig=="non-significant",], color = "black",size=2, shape=1) +
  geom_line(data=crlg, aes(color=property)) +
  geom_hline(yintercept = 0, linetype=2) +
  labs(x = "Phylogenetic distance class", y="Mantel correlation", color="Cryoconite property")+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=14),
        strip.text.x = element_text(size = 14))+
  ggtitle("Rare Community Mantel Correlogram")

#save pdf
pdf("rare-mantel-correlogram.pdf")
print(p)
dev.off()


# Rarefy to an even depth
set.seed(72)  # setting seed for reproducibility
# bulk.physeq.rare = rarefy_even_depth(bulk.physeq)
# Normalize read counts (this gives relative abundance)
new.ps.norm = transform_sample_counts(new.ps, function(x) x/sum(x))

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

# Get permutation beta null deviations
full.bNTI.df = Phylo_turnover(new.ps.norm, 10, 10)
#Save bNTI table
write.table(full.bNTI.df, "rare-full_bNTI.txt")



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
  
  null.dist = vegdist(t(simcom), method="bray")[1]
  return(null.dist)
}


# Main function for calculating RCbray
Calc_RCbray <- function(physeq, reps, nproc){
  # Get OTU table from phyloseq object
  otu.table = otu_table(physeq)
  
  # Get alpha diversity for each sample
  otu.PA.table = otu.table
  otu.PA.table[otu.PA.table > 0] = 1
  alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)
  
  # Get beta diversity matrix
  beta.table = as.matrix(vegdist(t(otu.PA.table), method="bray", diag=TRUE, upper=TRUE))
  
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

#get RC Bray results
RCb.df = Calc_RCbray(new.ps.norm, 2, 2)
#Save RCb.df table
write.table(RCb.df, "rare-RCb.df.txt")

#Plot results
x.1 = data.frame(sample_data(new.ps)) %>%
  select(SampleID, Glacier) %>%
  rename(Sample_1 = SampleID, Glacier_1 = Glacier) #get the corresponding glacier for each sample in sample_1

x.2 = data.frame(sample_data(new.ps)) %>%
  select(SampleID, Glacier) %>%
  rename(Sample_2 = SampleID, Glacier_2 = Glacier) #get the corresponding glacier for each sample in sample_2

#create dataframes
#make sure everything is formatted correctly for merging
RC_bray <- RCb.df
colnames(RC_bray) <- c("Sample_2", "Sample_1", "RCb")
RC_bray$Sample_1 <- as.integer(RC_bray$Sample_1)
RC_bray$Sample_2 <- as.integer(RC_bray$Sample_2)
full.bNTI.df$Sample_1 <- as.integer(full.bNTI.df$Sample_1)
full.bNTI.df$Sample_2 <- as.integer(full.bNTI.df$Sample_2)


bNTI.df = inner_join(full.bNTI.df, x.1) %>% #join the BNTI results with the right glacier names for samples1 and sample 2
  inner_join(x.2)
# Merge both datasets
turnover.df = inner_join(bNTI.df, RC_bray) #match the BNTI and RCbray results by sample name and merge
head(turnover.df)


#We need to assign each BNTI result to a particular glacier
# Make bNTI figure
#here we make sure we only draw comparisons between samples from the same site
#that's why Diamond Glacier and Upper Wright are removed from this analysis
within.bNTI.df = bNTI.df %>% #so get your BNTI dataframe
  filter(Glacier_1 == Glacier_2) %>% ##where the glaciers are the sample (i.e. the BNTI results are from samples within the same glacier)
  mutate(Glacier = Glacier_1) #give it an overall glacier name of that glacier. So this omits comparisons with samples from different sites

#Save within.bNTI.df
write.table(within.bNTI.df, "rare-within.bNTI.df.txt")

#plot the results!
eco.bNTI.plot = ggplot(within.bNTI.df, aes(x=Glacier, y=bNTI)) +
  geom_boxplot(outlier.shape=1) +
  geom_hline(yintercept = 2, linetype=2, size=0.5) +
  geom_hline(yintercept = -2, linetype=2, size=0.5) +
  labs(x="Glacier", y="BNTI") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=14))+
  ylim(-2, 10)+
  ggtitle("Rare Community BNTI Results")


#save pdf
pdf("rare-eco.bNTI.plot.pdf")
print(eco.bNTI.plot)
dev.off()


#Make the community analysis plot
#do the same as with did with JUST the BNTI dataframe with the BNTI and RCb dataframe
eco.turnover.df = turnover.df %>%
  filter(Glacier_1 == Glacier_2) %>%
  mutate(Glacier = Glacier_1)

eco.turnover.df = eco.turnover.df %>%
  mutate(process = ifelse(abs(bNTI) < 2, #if absolute value of BNTI is less than 2
                          ifelse(abs(RCb) < 0.95, "Drift", #and absolute value of RCb is less than 0.95 = Drift
                                 ifelse(RCb >= 0.95, "Dispersal Limited", #if RCb is > 0.95 = Dispersal limited
                                        ifelse(RCb <= -0.95, "Homogenizing Dispersal", "ERROR"))), #< -0.95 = homo disp
                          ifelse(bNTI >= 2, "Variable Selection", #more than 2 = variable selection
                                 ifelse(bNTI <= -2, "Homogeneous Selection", "ERROR")))) #< -2 = homogenous selection

eco.turnover.df$process = factor(eco.turnover.df$process, levels = c("Drift", 
                                                                     "Dispersal Limited", "Homogenizing Dispersal", 
                                                                     "Variable Selection", "Homogeneous Selection"))

#Save eco.turnover.df
write.table(eco.turnover.df, "rare-eco.turnover.df.txt")

#for each glacier
#find the number of number of site pairs
#and the percentage of those pairs that shows that process

glacier.list <- unique(eco.turnover.df$Glacier)

final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Glacier=character())

for (i in 1:length(glacier.list)){
  x <- eco.turnover.df[eco.turnover.df$Glacier == glacier.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Glacier <- glacier.list[i]
  final_data <- rbind(final_data, x)
}

#Save final_data
write.table(final_data, "rare-final_data.txt")

eco.turnover.plot = ggplot(final_data, aes(x=Glacier, y=perc, fill=process)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("white", "grey75", "grey50", "grey30", "black")) + 
  labs(x="Glacier", y="Percent of site pairs", fill="Process") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=14),
        legend.key.size = unit(10, "mm"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))+
  theme(aspect.ratio=4/10)+
  ggtitle("Rare Community RC Bray Results")


#save pdf
pdf("rare-eco.turnover.plot.pdf")
print(eco.turnover.plot)
dev.off()

# Merge the plots
eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot, 
                              rel_widths=c(0.6, 1), labels=c("A", "B"))
#save pdf
pdf("rare-eco.plot.pdf")
print(eco.plot)
dev.off() 


