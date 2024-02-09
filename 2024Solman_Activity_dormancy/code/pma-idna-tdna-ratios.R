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
library(gtools)


# 2. Import data
pro <- readRDS("../results/pma-16S-phylo-object-norm-rare.rds") 
euk <- readRDS("../results/pma-18S-phylo-object-norm-rare.rds") 


#16S

#1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.
#2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
#3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
#4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Step One: #1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.

#subset to tDNA samples only
noPMA = subset_samples(pro, Treatment == "tDNA")
#remove taxa with zero countrs
noPMA = filter_taxa(noPMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
noPMA = prune_samples(sample_sums(noPMA)>=1, noPMA)

#modify taxa names
noPMA.m = noPMA
x = data.frame(tax_table(noPMA.m)) %>% 
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

tax_table(noPMA.m) = as.matrix(y)

#get our data
data <- noPMA.m %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out = data %>%
  group_by(SampleID, Genus, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out2 = out %>%
  group_by(SampleID) %>%
  dplyr::summarise(NoPMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out2$ASV = out$OTU
out2$Genus = out$Genus

#sanity check
test = out2[out2$SampleID == "S21.25",]
sum(test$NoPMARelAbun) #should equal 1

#Step Two: #2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.

#subset to PMA samples only
PMA = subset_samples(pro, Treatment == "iDNA")
#remove taxa with zero countrs
PMA = filter_taxa(PMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
PMA = prune_samples(sample_sums(PMA)>=1, PMA)

#modify taxa names
PMA.m = PMA
x = data.frame(tax_table(PMA.m)) %>% 
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

tax_table(PMA.m) = as.matrix(y)

#get our data
data <- PMA.m %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out3 = data %>%
  group_by(SampleID, Genus, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out4 = out3 %>%
  group_by(SampleID) %>%
  dplyr::summarise(PMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out4$ASV = out3$OTU
out4$Genus = out3$Genus

#sanity check
test = out4[out4$SampleID == "S21.25",]
sum(test$PMARelAbun) #should equal 1


#Step Three: #3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.

#combine dataframes
df = full_join(out2, out4, by=c("SampleID", "ASV"))
#set NA values to zero
df$NoPMARelAbun[is.na(df$NoPMARelAbun)] <- 0
df$PMARelAbun[is.na(df$PMARelAbun)] <- 0
#get viability ratios
df$viability = df$PMARelAbun / df$NoPMARelAbun

#if an ASV has zero relative abundance in untreated samples (total community) omit them from the viability ratios
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
top15$ASV = top15$taxid

#remove any ASVs not in the top 15
top15 = top15[top15$tax_rank < 16,]

#only keep abundant ASVs
df.chop = df[df$ASV %in% top15$taxid,]

#give ASVs simple names
# df.chop.id <- transform(df.chop, id=match(ASV, unique(ASV)))
# df.chop.id = df.chop
df.chop.id = full_join(df.chop, top15[,c("tax_rank", "ASV")], by="ASV")

#change viability ratios greater than 1 to 1
df.chop.id$viability[df.chop.id$viability > 1] <- 1

#change viability ratios NaN to NA
df.chop.id$viability[is.nan(df.chop.id$viability)] <- NA

#log transform untreated relative abundances
df.chop.id$logNoPMARelAbun = log10(df.chop.id$NoPMARelAbun)

#log transform viability ratios
# df.chop.id$logViability = log(df.chop.id$viability)

#add Habitat data
df.chop.id$Habitat <- data$Habitat[match(df.chop.id$SampleID, data$SampleID)]

#set habitat as factor
df.chop.id$Habitat <- factor(df.chop.id$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#give them a new name
df.chop.id$name = paste0(df.chop.id$tax_rank, ": ", df.chop.id$Genus.x)

#make sure the names are in the right order
df.chop.id$name <- factor(df.chop.id$name, levels=c(mixedsort(unique(df.chop.id$name))))

# plot as heat map
p1 = ggplot(df.chop.id, aes(as.factor(name), SampleID, fill= logNoPMARelAbun)) + 
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
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-0.05),
        axis.title.x = element_blank(),
        plot.margin = margin(1,3,1,1, "cm"))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p1


p2 = ggplot(df.chop.id, aes(as.factor(name), SampleID, fill= viability)) + 
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
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-0.05),
        axis.title.x = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p2


p3 = (p1 + ggtitle("Prokaryote") + theme(plot.title = element_text(size=35))) + p2
p3

pdf("../results/pma-16S-dna-ratio.pdf", width=20, height=10)
print(p3)
dev.off()

#18S

#1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.
#2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.
#3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.
#4. For the 15 most abundant ASVs: Plot relative abundances from untreated samples and viability ratios.

#Step One: #1. Calculate the relative abundances of each ASV in untreated samples from each sample in each habitat type.

#subset to tDNA samples only
noPMA = subset_samples(euk, Treatment == "tDNA")
#remove taxa with zero countrs
noPMA = filter_taxa(noPMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
noPMA = prune_samples(sample_sums(noPMA)>=1, noPMA)

#modify taxa names
noPMA.m = noPMA
x = data.frame(tax_table(noPMA.m)) %>% 
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

tax_table(noPMA.m) = as.matrix(y)

#get our data
data <- noPMA.m %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out = data %>%
  group_by(SampleID, Genus, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out2 = out %>%
  group_by(SampleID) %>%
  dplyr::summarise(NoPMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out2$ASV = out$OTU
out2$Genus = out$Genus

#sanity check
test = out2[out2$SampleID == "S21.25",]
sum(test$NoPMARelAbun) #should equal 1

#Step Two: #2. Calculate the relative abundances of each ASV in treated samples from each sample in each habitat type.

#subset to PMA samples only
PMA = subset_samples(euk, Treatment == "iDNA")
#remove taxa with zero countrs
PMA = filter_taxa(PMA, function(x) sum(x) > 0, TRUE)
#remove samples with zero counts
PMA = prune_samples(sample_sums(PMA)>=1, PMA)

#modify taxa names
PMA.m = PMA
x = data.frame(tax_table(PMA.m)) %>% 
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

tax_table(PMA.m) = as.matrix(y)

#get our data
data <- PMA.m %>%
  psmelt() %>%
  as_tibble()

#group data by SampleID and ASV and 
out3 = data %>%
  group_by(SampleID, Genus, OTU) %>% #for each ASV in each sample
  dplyr::summarise(Abundance = sum(Abundance)) #get the total abundance

#from this get the relative abundance of each 
out4 = out3 %>%
  group_by(SampleID) %>%
  dplyr::summarise(PMARelAbun = Abundance/sum(Abundance))

#give ASV IDs
out4$ASV = out3$OTU
out4$Genus = out3$Genus

#sanity check
test = out4[out4$SampleID == "S21.25",]
sum(test$PMARelAbun) #should equal 1


#Step Three: #3. Divide relative abundance from treated samples by relative abundance from untreated samples to get viability ratio.

#combine dataframes
df = full_join(out2, out4, by=c("SampleID", "ASV"))
#set NA values to zero
df$NoPMARelAbun[is.na(df$NoPMARelAbun)] <- 0
df$PMARelAbun[is.na(df$PMARelAbun)] <- 0
#get viability ratios
df$viability = df$PMARelAbun / df$NoPMARelAbun

#if an ASV has zero relative abundance in untreated samples (total community) omit them from the viability ratios
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
top15$ASV = top15$taxid

#remove any ASVs not in the top 15
top15 = top15[top15$tax_rank < 16,]

#only keep abundant ASVs
df.chop = df[df$ASV %in% top15$taxid,]

#give ASVs simple names
# df.chop.id <- transform(df.chop, id=match(ASV, unique(ASV)))
# df.chop.id = df.chop
df.chop.id = full_join(df.chop, top15[,c("tax_rank", "ASV")], by="ASV")

#change viability ratios greater than 1 to 1
df.chop.id$viability[df.chop.id$viability > 1] <- 1

#change viability ratios NaN to NA
df.chop.id$viability[is.nan(df.chop.id$viability)] <- NA

#log transform untreated relative abundances
df.chop.id$logNoPMARelAbun = log10(df.chop.id$NoPMARelAbun)

#log transform viability ratios
# df.chop.id$logViability = log(df.chop.id$viability)

#add Habitat data
df.chop.id$Habitat <- data$Habitat[match(df.chop.id$SampleID, data$SampleID)]

#set habitat as factor
df.chop.id$Habitat <- factor(df.chop.id$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#give them a new name
df.chop.id$name = paste0(df.chop.id$tax_rank, ": ", df.chop.id$Genus.x)

#make sure the names are in the right order
df.chop.id$name <- factor(df.chop.id$name, levels=c(mixedsort(unique(df.chop.id$name))))

# plot as heat map
p4 = ggplot(df.chop.id, aes(as.factor(name), SampleID, fill= logNoPMARelAbun)) + 
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
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-0.05),
        axis.title.x = element_blank(),
        plot.margin = margin(1,3,1,1, "cm"))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1))
p4


p5 = ggplot(df.chop.id, aes(as.factor(name), SampleID, fill= viability)) + 
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
        strip.text.y = element_blank(),
        axis.text.x = element_text(angle=45, hjust=-0.05),
        axis.title.x = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.1, override.aes = list(size = 20)))
p5


p6 = (p4 + ggtitle("B: Microbial Eukaryote") + theme(plot.title = element_text(size=40),
                                                  strip.text = element_text(size=20),
                                                  legend.text = element_text(size=25),
                                                  legend.title = element_text(size=25),
                                                  legend.key.size = unit(3,"line"),
                                                  axis.text = element_text(size=18),
                                                  axis.title = element_text(size=18))) + p5 + theme(plot.title = element_text(size=40),
                                                                                          strip.text = element_text(size=20),
                                                                                          legend.text = element_text(size=25),
                                                                                          legend.title = element_text(size=25),
                                                                                          legend.key.size = unit(3,"line"),
                                                                                          axis.text = element_text(size=18),
                                                                                          axis.title = element_text(size=18))
p6

pdf("../results/pma-168S-dna-ratio.pdf", width=20, height=10)
print(p6)
dev.off()



#Combine plots

# p1a = p1 + theme(legend.position = "none")
# p2a = p2 + theme(legend.position = "none")
# p7 = (p1a + p2a) / (p4 + p5) + plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 20))

# p7 = ((p1 + ggtitle("Prokaryote") + theme(plot.title = element_text(size=35), legend.position = "none")) 
#       + (p2 + theme(legend.position = "none"))) / p6
# 
# pdf("../results/pma-16S-18S-dna-ratio.pdf", height=20, width=15)
# print(p7)
# dev.off()

p7 = ((p1 + ggtitle("A: Prokaryote") + theme(plot.title = element_text(size=40),
                                        strip.text = element_text(size=20),
                                        legend.text = element_text(size=16),
                                        axis.text = element_text(size=18),
                                        axis.title = element_text(size=18),
                                        legend.position = "none")) + (p2 + theme(plot.title = element_text(size=40),
                                                                                   strip.text = element_text(size=20),
                                                                                   legend.text = element_text(size=16),
                                                                                   axis.text = element_text(size=18),
                                                                                   axis.title = element_text(size=18),
                                                                                   legend.position = "none"))) / p6 

pdf("../results/pma-16S-18S-dna-ratio.pdf", height=25, width=18)
print(p7)
dev.off()
