#1. clear workspace and load packages
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
require("fantaxtic")

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-norm-rare.rds") 
#remove control samples
pro = subset_samples(pro, Habitat != "Control")
#remove taxa with zero countrs
pro = filter_taxa(pro, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
pro = prune_samples(sample_sums(pro)>=1, pro)

euk <- readRDS("../results/18S-phylo-object-norm-rare.rds") 
#remove control samples
euk = subset_samples(euk, Habitat != "Control")
#remove taxa with zero countrs
euk = filter_taxa(euk, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
euk = prune_samples(sample_sums(euk)>=1, euk)


#16S

#1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.
#2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
#3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
#4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Step One: #1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.

#subset to No-PMA samples only
noPMA = subset_samples(pro, Treatment == "No-PMA")
#remove taxa with zero countrs
noPMA = filter_taxa(noPMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
noPMA = prune_samples(sample_sums(noPMA)>=1, noPMA)

#get our data
data <- noPMA %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out = data %>%
  group_by(SampleID, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out2 = out %>%
  group_by(SampleID) %>%
  dplyr::summarise(NoPMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out2$ASV = out$OTU

#sanity check
test = out2[out2$SampleID == "S21.25",]
sum(test$NoPMARelAbun) #should equal 1

#Step Two: #2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.

#subset to PMA samples only
PMA = subset_samples(pro, Treatment == "PMA")
#remove taxa with zero countrs
PMA = filter_taxa(PMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
PMA = prune_samples(sample_sums(PMA)>=1, PMA)

#get our data
data <- PMA %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out3 = data %>%
  group_by(SampleID, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out4 = out3 %>%
  group_by(SampleID) %>%
  dplyr::summarise(PMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out4$ASV = out3$OTU

#sanity check
test = out4[out4$SampleID == "S21.25",]
sum(test$PMARelAbun) #should equal 1


#Step Three: #3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.

#combine dataframes
df = full_join(out2, out4, by=c("SampleID", "ASV"))
#set NA values to zero
df[is.na(df)] <- 0
#get viability ratios
df$viability = df$PMARelAbun / df$NoPMARelAbun

#if an ASV has zero relative abundance in untreated samples omit them from the viability ratios

df$viability[df$NoPMARelAbun == 0] <- NA

#Step Four: #4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Which are the 15 most abundant ASVs in the non-treated samples?

#firstly add ASV ID as a column in the tax table
tax <- data.frame(tax_table(noPMA))
tax$ASV = rownames(tax)
tax_table(noPMA) <- as.matrix(tax)

#find the 15 most abundant ASVs
top.phy <- top_taxa(noPMA, tax_level = "ASV", 
                    n_taxa = 15)
top15 = top.phy[[2]]

#only keep abundant ASVs
df.chop = df[df$ASV %in% top15$taxid,]

#give ASVs simple names
df.chop.id <- transform(df.chop, id=match(ASV, unique(ASV)))

#change viability ratios greater than 1 to 1
df.chop.id$viability[df.chop.id$viability > 1] <- 1

#change viability ratios NaN to NA
df.chop.id$viability[is.nan(df.chop.id$viability)] <- NA

#log transform untreated relative abundances
df.chop.id$logNoPMARelAbun = log(df.chop.id$NoPMARelAbun)

#log transform viability ratios
# df.chop.id$logViability = log(df.chop.id$viability)

#add Habitat data
df.chop.id$Habitat <- data$Habitat[match(df.chop.id$SampleID, data$SampleID)]

#set habitat as factor
df.chop.id$Habitat <- factor(df.chop.id$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

# plot as heat map
p1 = ggplot(df.chop.id, aes(as.factor(id), SampleID, fill= logNoPMARelAbun)) + 
  geom_tile() +
  facet_grid(rows = vars(Habitat), 
             scales = "free", 
             space = "free", 
             labeller = label_wrap_gen(width=10),
             switch = "y")+
  scale_fill_gradient(low="red", high="blue", na.value="white",name=expression("Relative Abundance Log"[10]))+
  scale_x_discrete(position = "top") +
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("ASV")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p1


p2 = ggplot(df.chop.id, aes(as.factor(id), SampleID, fill= viability)) + 
  geom_tile() +
  facet_grid(rows = vars(Habitat), 
             scales = "free", 
             space = "free", 
             labeller = label_wrap_gen(width=10))+
  scale_fill_gradient(low="grey", high="blue", na.value="white",
                      name="Viability Ratio iDNA/tDNA")+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right")+
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  xlab("ASV")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p2


p3 = p1 + p2 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
p3

pdf("../results/16S-dna-ratio.pdf")
print(p3)
dev.off()

#18S

#1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.
#2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
#3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
#4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Step One: #1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.

#subset to No-PMA samples only
noPMA = subset_samples(euk, Treatment == "No-PMA")
#remove taxa with zero countrs
noPMA = filter_taxa(noPMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
noPMA = prune_samples(sample_sums(noPMA)>=1, noPMA)

#get our data
data <- noPMA %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out = data %>%
  group_by(SampleID, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out2 = out %>%
  group_by(SampleID) %>%
  dplyr::summarise(NoPMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out2$ASV = out$OTU

#sanity check
test = out2[out2$SampleID == "S21.25",]
sum(test$NoPMARelAbun) #should equal 1

#Step Two: #2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.

#subset to PMA samples only
PMA = subset_samples(euk, Treatment == "PMA")
#remove taxa with zero countrs
PMA = filter_taxa(PMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
PMA = prune_samples(sample_sums(PMA)>=1, PMA)

#get our data
data <- PMA %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out3 = data %>%
  group_by(SampleID, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out4 = out3 %>%
  group_by(SampleID) %>%
  dplyr::summarise(PMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out4$ASV = out3$OTU

#sanity check
test = out4[out4$SampleID == "S21.25",]
sum(test$PMARelAbun) #should equal 1


#Step Three: #3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.

#combine dataframes
df = full_join(out2, out4, by=c("SampleID", "ASV"))
#set NA values to zero
df[is.na(df)] <- 0
#get viability ratios
df$viability = df$PMARelAbun / df$NoPMARelAbun

#if an ASV has zero relative abundance in untreated samples omit them from the viability ratios

df$viability[df$NoPMARelAbun == 0] <- NA

#Step Four: #4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Which are the 15 most abundant ASVs in the non-treated samples?

#firstly add ASV ID as a column in the tax table
tax <- data.frame(tax_table(noPMA))
tax$ASV = rownames(tax)
tax_table(noPMA) <- as.matrix(tax)

#find the 15 most abundant ASVs
top.phy <- top_taxa(noPMA, tax_level = "ASV", 
                    n_taxa = 15)
top15 = top.phy[[2]]

#only keep abundant ASVs
df.chop = df[df$ASV %in% top15$taxid,]

#give ASVs simple names
df.chop.id <- transform(df.chop, id=match(ASV, unique(ASV)))

#change viability ratios greater than 1 to 1
df.chop.id$viability[df.chop.id$viability > 1] <- 1

#change viability ratios NaN to NA
df.chop.id$viability[is.nan(df.chop.id$viability)] <- NA

#log transform untreated relative abundances
df.chop.id$logNoPMARelAbun = log(df.chop.id$NoPMARelAbun)

#add Habitat data
df.chop.id$Habitat <- data$Habitat[match(df.chop.id$SampleID, data$SampleID)]

#set habitat as factor
df.chop.id$Habitat <- factor(df.chop.id$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

# plot as heat map
p4 = ggplot(df.chop.id, aes(as.factor(id), SampleID, fill= logNoPMARelAbun)) + 
  geom_tile() +
  facet_grid(rows = vars(Habitat), 
             scales = "free", 
             space = "free", 
             labeller = label_wrap_gen(width=10),
             switch = "y")+
  scale_fill_gradient(low="red", high="blue", na.value="white",name=expression("Relative Abundance Log"[10]))+
  scale_x_discrete(position = "top") +
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  xlab("ASV")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p4


p5 = ggplot(df.chop.id, aes(as.factor(id), SampleID, fill= viability)) + 
  geom_tile() +
  facet_grid(rows = vars(Habitat), 
             scales = "free", 
             space = "free", 
             labeller = label_wrap_gen(width=10))+
  scale_fill_gradient(low="grey", high="blue", na.value="white",
                      name="Viability Ratio iDNA/tDNA")+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right")+
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  xlab("ASV")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p5


p6 = p4 + p5 + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
p6

pdf("../results/18S-dna-ratio.pdf")
print(p6)
dev.off()



#Combine plots

p1a = p1 + theme(legend.position = "none")
p2a = p2 + theme(legend.position = "none")
p7 = (p1a + p2a) / (p4 + p5) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))

pdf("../results/16S-18S-dna-ratio.pdf", height=10, width=8)
print(p7)
dev.off()
