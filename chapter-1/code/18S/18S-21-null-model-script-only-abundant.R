# ---
#   title: "Null Model Analysis"
# author: "Amy Solman"
# date: "12/10/2021"
# output: html_document
# ---
#   This script will carry out null model analysis on my communities. Code taken from https://github.com/seb369/landuse_comm_assembly 

# How do you carry out null model analysis?
  
#   A) Test for phylogenetic signal
# 1 Rarefy data (if this hasn't been done already) and get relative abundances
# 2 Edit metadata and fill missing data with median values
# 3 Calculate phylogenetic distances between ASVs
# 4 Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences
# 5 Construct mantel correlograms of phylogenetic distances and niche preferences
# 6 Plot to identify significant correlations
# 
# B) Compute Beta-nearest Taxon Index (βNTI)
# Rarefy data and get relative abundances 
# Define function for calculating the βMNTD for each random null community
# Define main function for calculating βNTI
# Carry out the βNTI calculation
#  
# C) Compute RC_Bray Model
# 1 Define function for calculating null values
# 2 Define main function for calculating RC Bray
# 3 Calculate RC_Bray
# 
# D) Plot results
# 1 Bind BNTI and RC_bray results into single dataframe
# 2 Plot BNTI
# 3 Prep data to plot RC Bray
# 4 Plot RC Bray
# 5 Merge plots


# Install Packages

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


#Read in data 

#ps.glaciers.merged.RA <- readRDS("../data-output/18S-my-phyloseq-object-rarefied.rds") #load rarefied phyloseq objectcol

#All glaciers merged abundances
ps.glaciers.merged.RA <- readRDS("../data-output/18S-all-glaciers-merged-rarefied-abundant.rds") #merged glaciers abundant phyloseq object
# ps.glaciers.merged.RR <- readRDS("../data-output/18S-all-glaciers-merged-rarefied-rare.rds") #merged glaciers rare phyloseq object
# ps.glaciers.merged.I <- readRDS("../data-output/18S-all-glaciers-merged-rarefied-intermediate.rds") #merged glaciers intermediate phyloseq object


# A) Test for phylogenetic signal
 
# 1 Rarefy data (if this hasn't been done already) and get relative abundances

set.seed(6666)
#I know my data is even sampling depth so I'm just going to get relative abundances
#get relative abundances
#rarefy to small sample depth so I can run this code on my computer
ps.glaciers.merged.RA <- rarefy_even_depth(ps.glaciers.merged.RA, sample.size = 2)

rel_abund <- transform_sample_counts(ps.glaciers.merged.RA, function(x) x/sum(x))


# 2 Edit metadata and fill missing data with median values

metadata <- data.frame(sample_data(ps.glaciers.merged.RA))

#keep only variables we want to include in our mode
meta_trim <- metadata[,17:92]
meta_trim <- meta_trim[,c(1:3, 5, 10, 12, 15, 18:33, 36, 39, 42, 45, 48, 51, 54:62)]
#meta_trim <- meta_trim[,1:3]
#calcualte area
meta_trim$area <- pi*(metadata$NS/2)*(metadata$EW/2)
#remove diamter info
drops <- c("NS","EW")
meta_trim <- meta_trim[ , !(names(meta_trim) %in% drops)]
#remove more unneccessary columns
meta_trim <- meta_trim[,-c(64:74)]
#make dataframe numeric
meta_trim[] <- lapply(meta_trim, as.numeric)

#Replace NAs with median values
f=function(x){
  x<-as.numeric(as.character(x)) #first convert each column into numeric if it is from factor
  x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}

meta=as.data.frame(apply(meta_trim,2,f))


# 3 Calculate phylogenetic distances between ASVs

#Get phylogenetic distances
tree <- phy_tree(rel_abund) #pull out phylogenetic tree
tree.dist <- cophenetic(tree)
asvs <- tree$tip.label
phylo.dists <- tree.dist[asvs, asvs]
phylo.dists[upper.tri(phylo.dists, diag=TRUE)] = NA

# 4 Calculate the weighted mean environmental parameter for each asv as a proxy for niche preferences

#give meta sample names as row names
rownames(meta) <- rownames(metadata) 

#turn the dataframe into a matrix
niches <- as.matrix(meta) 

asv.table <- data.frame(t(otu_table(rel_abund)), check.names = FALSE)
asv.niche <- wascores(niches, asv.table) #this will find the weighted mean environmental parameter for each asv
asv.niche.df <- data.frame(asv.niche, check.names = FALSE)


# 5 Construct mantel correlograms of phylogenetic distances and niche preferences

#Distance to Sea
Dist2Sea <- asv.niche.df$Distance_To_Sea #get Dist2Sea niche preferences
names(Dist2Sea) <- rownames(asv.niche.df) #match them with the asv rowname
Dist2Sea.dist <- as.matrix(dist(Dist2Sea), labels=TRUE) #convert into a matrix
x <- Dist2Sea.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Dist2Sea.corlg <- mantel.correlog(x, phylo.dists)

#Elevation
Elevation <- asv.niche.df$Elevation #get niche preferences
names(Elevation) <- rownames(asv.niche.df) #match them with the asv rowname
Elevation.dist <- as.matrix(dist(Elevation), labels=TRUE) #convert into a matrix
x <- Elevation.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Elevation.corlg <- mantel.correlog(x, phylo.dists)

#Total_Depth
Total_Depth <- asv.niche.df$Total_Depth #get niche preferences
names(Total_Depth) <- rownames(asv.niche.df) #match them with the asv rowname
Total_Depth.dist <- as.matrix(dist(Total_Depth), labels=TRUE) #convert into a matrix
x <- Total_Depth.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
Total_Depth.corlg <- mantel.correlog(x, phylo.dists)

#pH
pH <- asv.niche.df$pH #get niche preferences
names(pH) <- rownames(asv.niche.df) #match them with the asv rowname
pH.dist <- as.matrix(dist(pH), labels=TRUE) #convert into a matrix
x <- pH.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
pH.corlg <- mantel.correlog(x, phylo.dists)

#DOC
DOC <- asv.niche.df$DOC #get niche preferences
names(DOC) <- rownames(asv.niche.df) #match them with the asv rowname
DOC.dist <- as.matrix(dist(DOC), labels=TRUE) #convert into a matrix
x <- DOC.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
DOC.corlg <- mantel.correlog(x, phylo.dists)

#TDN
TDN <- asv.niche.df$TDN #get niche preferences
names(TDN) <- rownames(asv.niche.df) #match them with the asv rowname
TDN.dist <- as.matrix(dist(TDN), labels=TRUE) #convert into a matrix
x <- TDN.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
TDN.corlg <- mantel.correlog(x, phylo.dists)

#TDP
TDP <- asv.niche.df$TDP #get niche preferences
names(TDP) <- rownames(asv.niche.df) #match them with the asv rowname
TDP.dist <- as.matrix(dist(TDP), labels=TRUE) #convert into a matrix
x <- TDP.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
TDP.corlg <- mantel.correlog(x, phylo.dists)

#C.N
C.N <- asv.niche.df$C.N #get niche preferences
names(C.N) <- rownames(asv.niche.df) #match them with the asv rowname
C.N.dist <- as.matrix(dist(C.N), labels=TRUE) #convert into a matrix
x <- C.N.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
C.N.corlg <- mantel.correlog(x, phylo.dists)

#NO3.N_mueq
NO3.N_mueq <- asv.niche.df$NO3.N_mueq #get niche preferences
names(NO3.N_mueq) <- rownames(asv.niche.df) #match them with the asv rowname
NO3.N_mueq.dist <- as.matrix(dist(NO3.N_mueq), labels=TRUE) #convert into a matrix
x <- NO3.N_mueq.dist[asvs, asvs]
x[upper.tri(x, diag = TRUE)] = NA
NO3.N_mueq.corlg <- mantel.correlog(x, phylo.dists)

# 6 Plot to identify significant correlations

# Prep data
Dist2Sea.crlg <- data.frame(Dist2Sea.corlg$mantel.res) %>%
  mutate(property = "Dist2Sea")
Elevation.crlg <- data.frame(Elevation.corlg$mantel.res) %>%
  mutate(property = "Elevation")
Total_Depth.crlg <- data.frame(Total_Depth.corlg$mantel.res) %>%
  mutate(property = "Total_Depth")
pH.crlg <- data.frame(pH.corlg$mantel.res) %>%
  mutate(property = "pH")
DOC.crlg <- data.frame(DOC.corlg$mantel.res) %>%
  mutate(property = "DOC")
TDN.crlg <- data.frame(TDN.corlg$mantel.res) %>%
  mutate(property = "TDN")
TDP.crlg <- data.frame(TDP.corlg$mantel.res) %>%
  mutate(property = "TDP")
C.N.crlg <- data.frame(C.N.corlg$mantel.res) %>%
  mutate(property = "C.N")
NO3.N_mueq.crlg <- data.frame(NO3.N_mueq.corlg$mantel.res) %>%
  mutate(property = "NO3.N_mueq")

crlg <- rbind(Dist2Sea.crlg, Elevation.crlg, Total_Depth.crlg, pH.crlg, DOC.crlg, TDN.crlg, TDP.crlg, C.N.crlg, NO3.N_mueq.crlg) %>%
  mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
  filter(!(is.na(Pr.corrected.)))

#Save crlg table
write.table(crlg, "results/abundant-crlg.txt")

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
        strip.text.x = element_text(size = 14))

#save pdf
pdf("results/abundant-mantel-correlogram.pdf")
print(p)
dev.off()


#B) Compute Beta-nearest Taxon Index (βNTI)

#Rarefy data and get relative abundances 

# Rarefy to an even depth
#set.seed(72)  # setting seed for reproducibility
# bulk.physeq.rare = rarefy_even_depth(bulk.physeq)
# Normalize read counts (this gives relative abundance)
#ps.glaciers.merged.RA.norm = transform_sample_counts(ps.glaciers.merged.RA, function(x) x/sum(x))


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

# Carry out the BNTI calculation
# Get permutation beta null deviations
full.bNTI.df = Phylo_turnover(rel_abund, 10, 10)
# full.bNTI.df = Phylo_turnover(rel_abund, 1000, 10)
#Save bNTI table
write.table(full.bNTI.df, "results/abundant-full-bNTI.txt")

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
  
  null.dist = vegdist(t(simcom), method="bray")[1]
  return(null.dist)
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

# 3 Calculate RC_Bray

#get RC Bray results
RCb.df = Calc_RCbray(rel_abund, 5, 10)
# RCb.df = Calc_RCbray(rel_abund, 999, 20)
#Save RCb.df table
write.table(RCb.df, "results/abundant-RCb-df.txt")

#D) Plot Results

# 1 Bind BNTI and RC_bray results into single dataframe

x.1 = data.frame(sample_data(rel_abund)) %>%
  select(SampleID, Glacier) %>%
  rename(Sample_1 = SampleID, Glacier_1 = Glacier) #get the corresponding glacier for each sample in sample_1

x.2 = data.frame(sample_data(rel_abund)) %>%
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

# 2 Plot BNTI

#We need to assign each BNTI result to a particular glacier
# Make bNTI figure
#here we make sure we only draw comparisons between samples from the same site
#that's why Diamond Glacier and Upper Wright are removed from this analysis
within.bNTI.df = bNTI.df %>% #so get your BNTI dataframe
  filter(Glacier_1 == Glacier_2) %>% ##where the glaciers are the sample (i.e. the BNTI results are from samples within the same glacier)
  mutate(Glacier = Glacier_1) #give it an overall glacier name of that glacier. So this omits comparisons with samples from different sites

#Save within.bNTI.df
write.table(within.bNTI.df, "results/abundant-within-bNTI-df.txt")

#plot the results!
eco.bNTI.plot = ggplot(within.bNTI.df, aes(x=Glacier, y=bNTI)) +
  geom_boxplot(outlier.shape=1) +
  geom_hline(yintercept = 2, linetype=2, size=0.5) +
  geom_hline(yintercept = -2, linetype=2, size=0.5) +
  labs(x="Glacier", y=paste0("\u03B2NTI")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=14))+
  ylim(-2, 10)+
  ggtitle(paste0("Abundant Taxa \u03B2NTI Results"))


#save pdf
# pdf("results/abundant-eco-bNTI-plot.pdf")
# print(eco.bNTI.plot)
# dev.off()

ggsave(file="results/abundant-eco-bNTI-plot.svg", plot=eco.bNTI.plot)

# 3 Prep data to plot RC Bray

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
write.table(eco.turnover.df, "results/abundant-eco-turnover-df.txt")

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
write.table(final_data, "results/abundant-final-data.txt")

# 4 Plot RC Bray

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
  #theme(aspect.ratio=4/10)+
  ggtitle("Abundant Taxa RC Bray Results")


#save pdf
# pdf("results/abundant-eco-turnover-plot.pdf")
# print(eco.turnover.plot)
# dev.off()

ggsave(file="results/abundant-eco-turnover-plot.svg", plot=eco.turnover.plot)

#5 Merge plots

# Merge the plots
# eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot, 
#                               rel_widths=c(0.6, 1), labels=c("A", "B"))

eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot, 
                              labels=c("A", "B"), ncol=1)
#save pdf
# pdf("results/abundant-eco-plot.pdf")
# print(eco.plot)
# dev.off()

ggsave(file="results/abundant-eco-plot.svg", plot=eco.plot)
