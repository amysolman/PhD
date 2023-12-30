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


#Read in the results of our model
full.bNTI.df = read.csv("../results/bNTI-results-table.csv")
full.rcbray.df = read.csv("../results/RCbray-results-table.csv")

#remove first column
full.bNTI.df = full.bNTI.df[,!names(full.bNTI.df) %in% c("X")]
full.rcbray.df = full.rcbray.df[,!names(full.rcbray.df) %in% c("X")]

#GET RESULTS READY FOR PLOTTING

#create dataframes
#make sure everything is formatted correctly for merging
RC_bray <- full.rcbray.df
colnames(RC_bray) <- c("Sample_2", "Sample_1", "RCb", "Group", "Subcommunity", "Pole")
RC_bray$Sample_1 <- as.integer(RC_bray$Sample_1)
RC_bray$Sample_2 <- as.integer(RC_bray$Sample_2)
full.bNTI.df$Sample_1 <- as.integer(full.bNTI.df$Sample_1)
full.bNTI.df$Sample_2 <- as.integer(full.bNTI.df$Sample_2)

#make sure group columns are the same data type
full.bNTI.df$Group = as.character(full.bNTI.df$Group)
RC_bray$Group = as.character(RC_bray$Group)

#join dataframes
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

#save results
write.csv(final_data.edit, "../results/null-modell-plotting-results.csv")

#PLOT RESULTS
#final_data$round_perc <- round(final_data$perc, 2)

#remove rows with less than 1%
#x.cut <- x[which(x$frac > 0.01),]

#remove the full communitiy from plotting
final_data.edit.plot = subset(final_data.edit, Subcommunity != "Full")

#make group a factor
final_data.edit.plot$Group = factor(final_data.edit.plot$Group, levels=c("Prokaryote", "Eukaryote"))
#make subcommunity a factor
final_data.edit.plot$Subcommunity = factor(final_data.edit.plot$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))

#omit rows where the fraction is 0
final_data.edit.plot = final_data.edit.plot[final_data.edit.plot$frac > 0,]

p1 = ggplot(final_data.edit.plot, aes(x = Subcommunity,y = n_sites, 
                                      fill = factor(process, levels=c("Homogeneous Selection", "Variable Selection", "Homogenizing Dispersal", "Dispersal Limited", "Drift")))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels=scales::percent)+
  geom_text(
    aes(label = percent(frac)), position = position_fill(0.5), size=5) +
  facet_grid(~ Group + Pole)+
  theme_bw()+
  ylab("Relative Contribution %")+
  theme(text=element_text(size=15), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.title.x=element_blank(), 
        axis.text=element_text(size=12), 
        legend.text = element_text(size=20))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6", "#6ca8ab", "#517e80"))

print(p1)

pdf("../results/null-model.pdf", width=14, height=8)
print(p1)
dev.off()