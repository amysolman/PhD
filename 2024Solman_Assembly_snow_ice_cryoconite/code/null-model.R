# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(scales) #for percentages
library(picante) #for faith's pd and comdistnt functions
library(parallel)
library(svglite)

#Code taken from Barnett et al., (2020) was used for the model fitting (https://github.com/seb369/landuse_comm_assembly).  

ps.pro <- readRDS("../results/16S-phylo-object.rds") 
ps.euk <- readRDS("../results/18S-phylo-object-micro-remove.rds") 
# ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

#reduce the size of the ps objects for testing
# ps.pro = rarefy_even_depth(ps.pro, sample.size = min(sample_sums(ps.pro))/500, rngseed = 666)
# ps.euk = rarefy_even_depth(ps.euk, sample.size = min(sample_sums(ps.euk))/500, rngseed = 666)
# ps.mm = rarefy_even_depth(ps.mm, sample.size = min(sample_sums(ps.mm))/50, rngseed = 666)

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
# ps = ps.mm
# hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID
# group = "Micrometazoans"
# habitat = "Summer Ice"

get_bnti_res <- function(ps, hab, group, habitat){
  phylo = prune_samples(hab, ps)
  df = Phylo_turnover(transform_sample_counts(phylo, function(x) x/sum(x)), 1000, 10) #should be 1000 10
  df$Group = group
  df$Habitat = habitat
  return(df)
}

pro.sn.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
pro.sp.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice")
pro.sm.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice")
pro.cr.bNTI.df = get_bnti_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite")
euk.sn.bNTI.df = get_bnti_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Eukaryote", "Snow")
euk.sp.bNTI.df = get_bnti_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Eukaryote", "Spring Ice")
euk.sm.bNTI.df = get_bnti_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Eukaryote", "Summer Ice")
euk.cr.bNTI.df = get_bnti_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Eukaryote", "Cryoconite")
# mm.sn.bNTI.df = get_bnti_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Micrometazoan", "Snow")
# mm.sp.bNTI.df = get_bnti_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Micrometazoan", "Spring Ice")
# mm.sm.bNTI.df = get_bnti_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Micrometazoan", "Summer Ice")
# mm.cr.bNTI.df = get_bnti_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Micrometazoan", "Cryoconite")


full.bNTI.df = rbind(pro.sn.bNTI.df, pro.sp.bNTI.df, pro.sm.bNTI.df, pro.cr.bNTI.df,
                     euk.sn.bNTI.df, euk.sp.bNTI.df, euk.sm.bNTI.df, euk.cr.bNTI.df)

write.csv(full.bNTI.df, "../results/bNTI-results-table.csv")

get_rcbray_res <- function(ps, hab, group, habitat){
  
  phylo = prune_samples(hab, ps)
  df = Calc_RCbray(phylo, 999, 20) #should be 999 20
  df$Group = group
  df$Habitat = habitat
  return(df)
  
}

pro.sn.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, "Prokaryote", "Snow")
pro.sp.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, "Prokaryote", "Spring Ice")
pro.sm.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, "Prokaryote", "Summer Ice")
pro.cr.rcbray.df = get_rcbray_res(ps.pro, data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, "Prokaryote", "Cryoconite")
euk.sn.rcbray.df = get_rcbray_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, "Eukaryote", "Snow")
euk.sp.rcbray.df = get_rcbray_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, "Eukaryote", "Spring Ice")
euk.sm.rcbray.df = get_rcbray_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, "Eukaryote", "Summer Ice")
euk.cr.rcbray.df = get_rcbray_res(ps.euk, data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, "Eukaryote", "Cryoconite")
# mm.sn.rcbray.df = get_rcbray_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, "Micrometazoan", "Snow")
# mm.sp.rcbray.df = get_rcbray_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, "Micrometazoan", "Spring Ice")
# mm.sm.rcbray.df = get_rcbray_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, "Micrometazoan", "Summer Ice")
# mm.cr.rcbray.df = get_rcbray_res(ps.mm, data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, "Micrometazoan", "Cryoconite")

full.rcbray.df = rbind(pro.sn.rcbray.df, pro.sp.rcbray.df, pro.sm.rcbray.df, pro.cr.rcbray.df,
                       euk.sn.rcbray.df, euk.sp.rcbray.df, euk.sm.rcbray.df, euk.cr.rcbray.df)
# mm.sn.rcbray.df, mm.sp.rcbray.df, mm.sm.rcbray.df, mm.cr.rcbray.df)

write.csv(full.rcbray.df, "../results/RCbray-results-table.csv")


#GET RESULTS READY FOR PLOTTING

#create dataframes
#make sure everything is formatted correctly for merging
RC_bray <- full.rcbray.df
colnames(RC_bray) <- c("Sample_2", "Sample_1", "RCb", "Group", "Habitat")

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

#for each habitat
#find the number of number of site pairs
#and the percentage of those pairs that shows that process

hab.list <- unique(turnover.df$Habitat)

#Prokaryotes
pro_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Habitat=character())
new.df <- subset(turnover.df, Group == "Prokaryote")

for (i in 1:length(hab.list)){
  x <- new.df[new.df$Habitat == hab.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Habitat<- hab.list[i]
  pro_final_data <- rbind(pro_final_data, x)
}
pro_final_data$Group = "Prokaryote"

#Eukaryotes
euk_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Habitat=character())
new.df <- subset(turnover.df, Group == "Eukaryote")

for (i in 1:length(hab.list)){
  x <- new.df[new.df$Habitat == hab.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Habitat<- hab.list[i]
  euk_final_data <- rbind(euk_final_data, x)
}
euk_final_data$Group = "Eukaryote"

#Micrometazoans
# mm_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Habitat=character())
# new.df <- subset(turnover.df, Group == "Micrometazoan")
# 
# for (i in 1:length(hab.list)){
#   x <- new.df[new.df$Habitat == hab.list[i], ]
#   x <- x %>%
#     group_by(process) %>%
#     dplyr::summarize(n_sites = n(),
#                      perc=(n()/nrow(x))*100) %>%
#     as.data.frame
#   x$Habitat<- hab.list[i]
#   mm_final_data <- rbind(mm_final_data, x)
# }
# mm_final_data$Group = "Micrometazoan"


#bind results together
final_data = rbind(pro_final_data, euk_final_data) #, mm_final_data)

#add percentages column
final_data.edit <- final_data %>% dplyr::group_by(Habitat, Group) %>%
  dplyr::mutate(frac = n_sites / sum(n_sites))
final_data.edit$frac = round(final_data.edit$frac,2)

#save dataframe
write.csv(final_data.edit, "../results/null-model-results.csv")

#PLOT RESULTS

# final_data.edit = final_data.edit[final_data.edit$Group != "Micrometazoan",]

#make sure everything is in the right order
final_data.edit$Group = factor(final_data.edit$Group, 
                               levels = c("Prokaryote", "Eukaryote"))
final_data.edit$Habitat = factor(final_data.edit$Habitat,
                                 levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))


p1 = ggplot(final_data.edit, aes(x = Habitat,y = n_sites, 
                                 fill = factor(process, levels=c("Homogeneous Selection", "Variable Selection", "Homogenizing Dispersal", "Dispersal Limited", "Drift")))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels=scales::percent)+
  geom_text(
    aes(label = percent(frac)), position = position_fill(0.5), size=4) +
  facet_grid(~ Group)+
  theme_bw()+
  ylab("Relative Contribution %")+
  theme(text=element_text(size=10), legend.position = "bottom", legend.title = element_blank(), axis.title.x=element_blank(), axis.text=element_text(size=10), legend.text = element_text(size=10), strip.text.x = element_text(size = 10))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  scale_fill_manual(values=c("darkorange", "darkorange3", "cyan", "cyan3", "cyan4"))

print(p1)

pdf("../results/null-model.pdf", height = 10, width = 12)
print(p1)
dev.off()
