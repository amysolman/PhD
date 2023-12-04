rm(list=ls())
library(ggplot2)
library(phyloseq)
library(patchwork)
library(ggpubr)
library(glue)
library(tibble) # for add_column
library(reshape2)
library(cowplot)
library(betapart)
library(dplyr)

# 3. Import data and separate into Arctic and Antarctic samples
# 2. Import data
pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds") 

###################################################################################################
###################################################################################################

#Beta diversity partitioning analysis 
#calculate total beta diversity and beta diversity due to richness and turnover and plot

# phylo = pro_res[[7]] #arctic abundant community
# subcom = "Rare"
# data = "Prokaryote"
# pole = "Arctic"

beta_part <- function(phylo, subcom, data, pole){
  
  set.seed(666)
  
  #normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)
  
  #for each subcommunity get the count table with rows as samples and species as columns
  counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)
  
  #max(rowSums(counts)) - min(rowSums(counts))
  
  #beta partition analysis
  out = bray.part(counts)
  
  #convert the matrices into a dataframe to plot
  out.turn = melt(as.matrix(out$bray.bal))
  names(out.turn) = c("Sample1", "Sample2", "Beta")
  out.turn$Part = "Turnover"
  out.rich = melt(as.matrix(out$bray.gra))
  names(out.rich) = c("Sample1", "Sample2", "Beta")
  out.rich$Part = "Richness"
  out.whole = melt(as.matrix(out$bray))
  names(out.whole) = c("Sample1", "Sample2", "Beta")
  out.whole$Part = "Whole"
  
  #bind data frame
  three.df = rbind(out.whole, out.turn, out.rich)
  
  #remove samples being compared to themseves
  three.df.trim = three.df[three.df['Sample1'] != three.df['Sample2'],]
  
  #add abundance class
  three.df.trim$Subcommunity = subcom
  
  #and dataset
  three.df.trim$Data = data
  three.df.trim$Pole = pole
  
  #get proportions
  two.df = rbind(out.turn, out.rich)
  #remove samples being compared to themseves
  two.df.trim = two.df[two.df['Sample1'] != two.df['Sample2'],]
  #get proportions
  two.df.trim$Proportion = two.df.trim$Beta/sum(two.df.trim$Beta)
  sum(two.df.trim$Proportion)
  two.df.trim$Subcommunity = subcom
  two.df.trim$Data = data
  two.df.trim$Pole = pole
  
  res.list = list(three.df.trim, two.df.trim)
  
  return(res.list)
  
}

ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.abun.beta.part.arc = beta_part(ps, "Abundant", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.int.beta.part.arc = beta_part(ps, "Intermediate", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.rare.beta.part.arc = beta_part(ps, "Rare", "Prokaryote", "Arctic")

ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.abun.beta.part.ant = beta_part(ps, "Abundant", "Prokaryote", "Antarctic")

ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.int.beta.part.ant = beta_part(ps, "Intermediate", "Prokaryote", "Antarctic")

ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
pro.rare.beta.part.ant = beta_part(ps, "Rare", "Prokaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.abun.beta.part.arc = beta_part(ps, "Abundant", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.int.beta.part.arc = beta_part(ps, "Intermediate", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.rare.beta.part.arc = beta_part(ps, "Rare", "Eukaryote", "Arctic")

ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.abun.beta.part.ant = beta_part(ps, "Abundant", "Eukaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.int.beta.part.ant = beta_part(ps, "Intermediate", "Eukaryote", "Antarctic")

ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
euk.rare.beta.part.ant = beta_part(ps, "Rare", "Eukaryote", "Antarctic")

#make into one dataframe
final.beta = rbind(pro.abun.beta.part.arc[[1]], pro.int.beta.part.arc[[1]], pro.rare.beta.part.arc[[1]],
                   pro.abun.beta.part.ant[[1]], pro.int.beta.part.ant[[1]], pro.rare.beta.part.ant[[1]],
                   euk.abun.beta.part.arc[[1]], euk.int.beta.part.arc[[1]], euk.rare.beta.part.arc[[1]],
                   euk.abun.beta.part.ant[[1]], euk.int.beta.part.ant[[1]], euk.rare.beta.part.ant[[1]])

final.beta.barplot = rbind(pro.abun.beta.part.arc[[2]], pro.int.beta.part.arc[[2]], pro.rare.beta.part.arc[[2]],
                           pro.abun.beta.part.ant[[2]], pro.int.beta.part.ant[[2]], pro.rare.beta.part.ant[[2]],
                           euk.abun.beta.part.arc[[2]], euk.int.beta.part.arc[[2]], euk.rare.beta.part.arc[[2]],
                           euk.abun.beta.part.ant[[2]], euk.int.beta.part.ant[[2]], euk.rare.beta.part.ant[[2]])


# ggplot() + 
#    geom_histogram(aes(Petal.Width))+ 
#    facet_grid(Species~.)

#using KRuskal-Wallist Test are there significant differences in the proportion of dissimilarity account for Turnover and Richness between the subcommunities?
#Test Turnover
turn.df.pro.arc = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(turn.df.pro.arc$Beta, turn.df.pro.arc$Subcommunity,
                     p.adjust.method = "BH")
turn.df.pro.ant = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(turn.df.pro.ant$Beta, turn.df.pro.ant$Subcommunity,
                     p.adjust.method = "BH")
turn.df.euk.arc = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(turn.df.euk.arc$Beta, turn.df.euk.arc$Subcommunity,
                     p.adjust.method = "BH")
turn.df.euk.ant = final.beta[final.beta$Part == "Turnover" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(turn.df.euk.ant$Beta, turn.df.euk.ant$Subcommunity,
                     p.adjust.method = "BH")

rich.df.pro.arc = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(rich.df.pro.arc$Beta, rich.df.pro.arc$Subcommunity,
                     p.adjust.method = "BH")
rich.df.pro.ant = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Prokaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(rich.df.pro.ant$Beta, rich.df.pro.ant$Subcommunity,
                     p.adjust.method = "BH")
rich.df.euk.arc = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Arctic",]
pairwise.wilcox.test(rich.df.euk.arc$Beta, rich.df.euk.arc$Subcommunity,
                     p.adjust.method = "BH")
rich.df.euk.ant = final.beta[final.beta$Part == "Richness" & final.beta$Data == "Eukaryote" & final.beta$Pole == "Antarctic",]
pairwise.wilcox.test(rich.df.euk.ant$Beta, rich.df.euk.ant$Subcommunity,
                     p.adjust.method = "BH")

final.beta %>%
  filter(Part == "Whole") %>%
  mutate(across(Data, factor, levels=c("Prokaryote","Eukaryote"))) %>%
  mutate(across(Pole, factor, levels=c("Arctic","Antarctic"))) %>%
  ggplot(aes(x=Subcommunity, y=Beta, fill=Part)) +
  geom_boxplot(position=position_dodge(1))+
  facet_wrap(~ Data + Pole, scales = "free_x")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Bray-Curtis Dissimilarity")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        legend.text = element_text(size=10), 
        legend.title = element_blank())

p = final.beta.barplot %>%
  mutate(across(Data, factor, levels=c("Prokaryote","Eukaryote"))) %>%
  mutate(across(Pole, factor, levels=c("Arctic","Antarctic"))) %>%
  mutate(across(Subcommunity, factor, levels=c("Rare", "Intermediate", "Abundant"))) %>%
  group_by(Part, Subcommunity, Pole, Data) %>%
  dplyr::summarise(across(c(Proportion), sum)) %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.y=element_text(size=10,face="bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())

pdf("../results/beta-diversity-partitioning.pdf")
print(p)
dev.off()


#get percentages from final.beta.barplot dataframe

#merge proportion by part, subcommunity, data and pole

df = final.beta.barplot %>% group_by(Part, Subcommunity, Data, Pole) %>% summarise_each(funs(sum))

#round(df[df$Part == "Richness" & df$Pole == "Arctic" & df$Subcommunity == "Abundant" & df$Data == "Prokaryote",]$Proportion*100, 2)

x = data.frame(cbind(df$Part, df$Subcommunity, df$Data, df$Pole, round(df$Proportion*100, 2)))
names(x) = c("Partition", "Subcommunity", "Dataset", "Pole", "Proportion")

rownames(x) <- NULL
knitr::kable(x, caption = 'Table 7. Beta-diversity partitioning analysis results.')

###Beta diversity partitioning analysis

df.keep = data.frame()

#Abundant Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Arctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())






#Abundant Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Antarctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.x=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())






#Abundant Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Abundant" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Eukaryote"
pole = "Arctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Intermediate" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Eukaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Eukaryotes
ps = subset_samples(euk, Subcommunity == "Rare" & Pole == "Arctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Eukaryote"
pole = "Arctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)

x = df.keep[df.keep$Data == "Eukaryote" & df.keep$Part == "Richness",]

#barchart

p = df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())
pdf("../results/test-barplot.pdf", height = 30)
print(p)
dev.off()


#Abundant Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Abundant" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Abundant"
data = "Prokaryote"
pole = "Antarctic"

set.seed(666)

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)






#Intermediate Antarctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Intermediate" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Intermediate"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)





#Rare Arctic Prokaryotes
ps = subset_samples(pro, Subcommunity == "Rare" & Pole == "Antarctic")
ps = filter_taxa(ps, function(x) sum(x) >= 1, TRUE)
ps = prune_samples(sample_sums(ps)>=1, ps)
phylo = ps
subcom ="Rare"
data = "Prokaryote"
pole = "Antarctic"

#normalisation of the data does not appear to work with this analyses. Have looked at other implementations of bray.part and the subcommunities were not normalised prior to this step (https://github.com/shuojiao/Rare-and-abundant-bacteria-/blob/master/Beta-diversity%20and%20CRT.R)

#get the count table with rows as samples and species as columns
counts <- data.frame(t(otu_table(phylo)), check.names=FALSE)

#beta partition analysis
out = bray.part(counts)

#convert the matrices into a dataframe to plot
out.turn = melt(as.matrix(out$bray.bal))
names(out.turn) = c("Sample1", "Sample2", "Beta")
out.turn$Part = "Turnover"
out.rich = melt(as.matrix(out$bray.gra))
names(out.rich) = c("Sample1", "Sample2", "Beta")
out.rich$Part = "Richness"
out.whole = melt(as.matrix(out$bray))
names(out.whole) = c("Sample1", "Sample2", "Beta")
out.whole$Part = "Whole"

#this give us the WHOLE beta diversity (Bray-Curtis) between each pair of samples,
#the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)

#bind data frame
df = rbind(out.whole, out.turn, out.rich) 
#remove samples being compared to themseves
df.trim = df[df['Sample1'] != df['Sample2'],]
#add abundance class
df.trim$Subcommunity = subcom
#and dataset
df.trim$Data = data
df.trim$Pole = pole


#this give us the beta diversity attributed to changes in richness (e.g. the same species just more or less)
#the beta diversity attributed to turnover (e.g. different species)
#and the proportion of the dataset attributed to those things.

#get proportions
df2 = rbind(out.turn, out.rich)
#remove samples being compared to themseves
df2.trim = df2[df2['Sample1'] != df2['Sample2'],]
#get proportions
df2.trim$Proportion = df2.trim$Beta/sum(df2.trim$Beta)
sum(df2.trim$Proportion) == 1 #should be TRUE
df2.trim$Subcommunity = subcom
df2.trim$Data = data
df2.trim$Pole = pole

df.keep = rbind(df.keep, df2.trim)


#barchart

df.keep %>%
  ggplot(aes(x=Subcommunity, y=Proportion, fill=Part)) +
  facet_wrap(~ Data + Pole, scales = "free_x")+
  geom_bar(stat="identity")+
  theme_bw()+
  guides(fill=guide_legend(title="Partition"))+
  ylab("Proportion %")+
  scale_colour_brewer(palette="Set1")+
  theme(legend.position = "bottom",
        axis.text=element_text(size=10),
        axis.title.x=element_blank(),
        legend.text = element_text(size=10), 
        legend.title = element_blank())