# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

source("00-solman-functions.R")
library(phyloseq)
# library(RColorBrewer)
# library(tidyr) #trans data between wide and long format
library(ggplot2)
library(vegan) #for vegdist
# source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
# library(dplyr)
library(cowplot) #plot_grid
library(reshape2) #for function melt
# library(funrar) #for make relative
# library(stringr) #mutate function, to split columns
# library(gridExtra) #for exporting as pdf
# library(scales)
# library(microbiome) #for summarize_phyloseq
# library(picante) #for faith's pd
library(ggpmisc) #for adding lm stats to the plot
library(broom) #for adding stats to the plot
library(fossil) #for earth.dist function
library(patchwork)
library(ggpubr)
library(dplyr)



#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds") 

# ps = ps.pro
# sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID
# col = "#d8a4ff"

ddr <- function(ps, sam, col){
  
  #Separate habitat data
  h <- prune_samples(sam, ps)
  
  #remove empty samples or ASVs
  h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
  h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)
  
  #Calculate Bray-Curtis dissimilarities between sites
  bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")
  
  #Calculate geographic distance between sites
  ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
  geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers
  
  #transform these data into columns
  comm.dist = melt(as.matrix(bc))
  names(comm.dist) = c("Sample1", "Sample2", "Beta")
  comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]
  
  geo.dist = melt(as.matrix(geo.dist))
  names(geo.dist) = c("Sample1", "Sample2", "km")
  geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]
  
  #make into dataframe
  data = comm.dist.trim 
  data$km = geo.dist.trim$km
  
  plot(data$km, data$Beta) #we can see a somewhat linear relationship
  cor(data$km, data$Beta, method="spearman") #spearman's correlation
  
  #log transform our data before fitting the linear model
  data$log10km <- log10(data$km)
  data$log10Beta <- log10(data$Beta)
  
  #fit the model  
  model <- lm(log10Beta ~ log10km, data = data )
  
  # #model results
  summary(model)
  
  plot(data$log10km, data$log10Beta)
  
  #plot
  
  p  =  ggplot(data, aes(x=log10km, y=log10Beta)) +
    geom_point(data=data , aes(x=log10km, y=log10Beta, fill=col),
               alpha=.4, size=3, shape=21)+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE, show.legend = FALSE)+
    theme_bw()+
    xlab("log10(km)")+
    ylab("log10(Bray-Curtis Dissimilarity)")+
    # scale_colour_brewer(palette="Set1")+
    #     theme(legend.position = "bottom",
    #       axis.text=element_text(size=10),
    #       axis.title=element_text(size=10,face="bold"),
    #       legend.text = element_text(size=10), 
    #         legend.title = element_blank())+
    scale_fill_manual(values=c(col))+
    guides(color = guide_legend(override.aes = list(size = 4, alpha=.5) ) )+
    stat_poly_eq(aes(label =  paste(stat(eq.label), "*\" with \"*", 
                                    #stat(rr.label), "*\", \"*", 
                                    stat(adj.rr.label), "*\", and \"*", 
                                    #stat(f.value.label), "*\", and \"*",
                                    stat(p.value.label), "*\".\"",
                                    sep = "")),
                 formula = y ~ x, parse = TRUE, size = 3, label.y = "bottom")
  
  print(p)
  
  res = list(model, p)
  
  return(res)
  
}


pro.snow <- ddr(ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, col="#fa9f99")
pro.sp <- ddr(ps=ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, col="#a4c64d")
pro.sum <- ddr(ps=ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, col="#4dd2d6")
pro.cry <- ddr(ps=ps.pro, sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, col="#d8a4ff")

summary(pro.snow[[1]])
summary(pro.sp[[1]])
summary(pro.sum[[1]])
summary(pro.cry[[1]])

print(pro.snow[[2]])
print(pro.sp[[2]])
print(pro.sum[[2]])
print(pro.cry[[2]])

# pro.p <- plot_grid(pro.snow[[2]],
#                    pro.sp[[2]],
#                    pro.sum[[2]],
#                    pro.cry[[2]],
#                    labels="AUTO")
# pro.p

pro.p = egg::ggarrange(pro.snow[[2]] + ggtitle("Snow") + ylim(-1.8, 0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none"),
                       pro.sp[[2]] + ggtitle("Spring Ice") + ylim(-1.8, 0) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none"),
                       pro.sum[[2]] + ggtitle("Summer Ice") + ylim(-1.8, 0) + theme(legend.position = "none"),
                       pro.cry[[2]] + ggtitle("Cryoconite") + ylim(-1.8, 0)+ theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"))

pdf("../results/prokaryote-distance-decay.pdf", width=10, height=6)
print(pro.p)
dev.off()

euk.snow <- ddr(ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, col="#fa9f99")
euk.sp <- ddr(ps=ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, col="#a4c64d")
euk.sum <- ddr(ps=ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, col="#4dd2d6")
euk.cry <- ddr(ps=ps.euk, sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, col="#d8a4ff")

summary(euk.snow[[1]])
summary(euk.sp[[1]])
summary(euk.sum[[1]])
summary(euk.cry[[1]])

print(euk.snow[[2]])
print(euk.sp[[2]])
print(euk.sum[[2]])
print(euk.cry[[2]])

# euk.p <- plot_grid(euk.snow[[2]],
#                    euk.sp[[2]],
#                    euk.sum[[2]],
#                    euk.cry[[2]],
#                    labels="AUTO")
# euk.p

euk.p = egg::ggarrange(euk.snow[[2]] + ggtitle("Snow") + ylim(-1.8, 0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none"),
                       euk.sp[[2]] + ggtitle("Spring Ice") + ylim(-1.8, 0) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none"),
                       euk.sum[[2]] + ggtitle("Summer Ice") + ylim(-1.8, 0) + theme(legend.position = "none"),
                       euk.cry[[2]] + ggtitle("Cryoconite") + ylim(-1.8, 0)+ theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"))

pdf("../results/eukaryote-distance-decay.pdf", width=10, height=6)
print(euk.p)
dev.off()

mm.snow <- ddr(ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, col="#fa9f99")
mm.sp <- ddr(ps=ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, col="#a4c64d")
mm.sum <- ddr(ps=ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, col="#4dd2d6")
mm.cry <- ddr(ps=ps.mm, sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, col="#d8a4ff")


summary(mm.snow[[1]])
summary(mm.sp[[1]])
summary(mm.sum[[1]])
summary(mm.cry[[1]])

print(mm.snow[[2]])
print(mm.sp[[2]])
print(mm.sum[[2]])
print(mm.cry[[2]])

# mm.p <- plot_grid(mm.snow[[2]] + ggtitle("Snow") + ylim(-1.8, 0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
#                    mm.sp[[2]] + ggtitle("Spring Ice") + ylim(-1.8, 0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
#                    mm.sum[[2]] + ggtitle("Summer Ice") + ylim(-1.8, 0),
#                    mm.cry[[2]] + ggtitle("Cryoconite") + ylim(-1.8, 0))
# mm.p

mm.p = egg::ggarrange(mm.snow[[2]] + ggtitle("Snow") + ylim(-1.8, 0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none"),
                      mm.sp[[2]] + ggtitle("Spring Ice") + ylim(-1.8, 0) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), legend.position = "none"),
                      mm.sum[[2]] + ggtitle("Summer Ice") + ylim(-1.8, 0) +theme(legend.position = "none"),
                      mm.cry[[2]] + ggtitle("Cryoconite") + ylim(-1.8, 0)+ theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none"))

pdf("../results/micrometazoan-distance-decay.pdf", width=10, height=6)
print(mm.p)
dev.off()

#Combine plots


#as_ggplot(pro.p) | as_ggplot(euk.p) | as_ggplot(mm.p)

thm <- theme(plot.title = element_text(face = 2, size = 14))
first_plot      <- wrap_elements(as_ggplot(pro.p) + plot_annotation(title = "Prokaryote", theme = thm))
second_plot   <- wrap_elements(as_ggplot(euk.p) + plot_annotation(title = "Microbial eukaryote", theme = thm))
third_plot   <- wrap_elements(as_ggplot(mm.p) + plot_annotation(title = "Microfauna", theme = thm))

all.p = first_plot | second_plot | third_plot
all.p

pdf("../results/all-distance-decay.pdf", width=25, height=10)
print(all.p)
dev.off()

#Looking at different spatial scales
#Prokaryotes

#####################PRO SNOW#####################
ps = ps.pro
sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#Calculate environmental dissimilarity between sites
chem <- data.frame(sample_data(h.f))
chem <- chem[,c(9,10,13:ncol(chem))]
#z-score transform the environmental data
chem.z = data.frame(scale(chem))
#Get get distance matrix for environmental factors
chem.dist<-vegdist(chem.z, method="euclidean", na.rm=TRUE)

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

chem.dist = melt(as.matrix(chem.dist))
names(chem.dist) = c("Sample1", "Sample2", "ChemDist")
chem.dist.trim = chem.dist[chem.dist['Sample1'] != chem.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000
data$ChemDist = chem.dist.trim$ChemDist

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation
cor(data$m, data$ChemDist, method="spearman") #spearman's correlation
cor(data$ChemDist, data$Beta, method="spearman") #spearman's correlation


####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
snow.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
snow.test.data$Habitat = "Snow"



#####################PRO SPRING ICE#####################
ps = ps.pro
sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#Calculate environmental dissimilarity between sites
chem <- data.frame(sample_data(h.f))
chem <- chem[,c(9,10,13:ncol(chem))]
#z-score transform the environmental data
chem.z = data.frame(scale(chem))
#Get get distance matrix for environmental factors
chem.dist<-vegdist(chem.z, method="euclidean", na.rm=TRUE)

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

chem.dist = melt(as.matrix(chem.dist))
names(chem.dist) = c("Sample1", "Sample2", "ChemDist")
chem.dist.trim = chem.dist[chem.dist['Sample1'] != chem.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000
data$ChemDist = chem.dist.trim$ChemDist

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation
cor(data$m, data$ChemDist, method="spearman") #spearman's correlation
cor(data$ChemDist, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sp.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sp.test.data$Habitat = "Spring Ice"


#####################PRO SUMMER ICE#####################
ps = ps.pro
sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#Calculate environmental dissimilarity between sites
chem <- data.frame(sample_data(h.f))
chem <- chem[,c(9,10,13:ncol(chem))]
#z-score transform the environmental data
chem.z = data.frame(scale(chem))
#Get get distance matrix for environmental factors
chem.dist<-vegdist(chem.z, method="euclidean", na.rm=TRUE)

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

chem.dist = melt(as.matrix(chem.dist))
names(chem.dist) = c("Sample1", "Sample2", "ChemDist")
chem.dist.trim = chem.dist[chem.dist['Sample1'] != chem.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000
data$ChemDist = chem.dist.trim$ChemDist

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation
cor(data$m, data$ChemDist, method="spearman") #spearman's correlation
cor(data$ChemDist, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sum.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sum.test.data$Habitat = "Summer Ice"

#####################PRO CRYO#####################

ps = ps.pro
sam = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#Calculate environmental dissimilarity between sites
chem <- data.frame(sample_data(h.f))
chem <- chem[,c(9,10,13:ncol(chem))]
#z-score transform the environmental data
chem.z = data.frame(scale(chem))
#Get get distance matrix for environmental factors
chem.dist<-vegdist(chem.z, method="euclidean", na.rm=TRUE)

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

chem.dist = melt(as.matrix(chem.dist))
names(chem.dist) = c("Sample1", "Sample2", "ChemDist")
chem.dist.trim = chem.dist[chem.dist['Sample1'] != chem.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000
data$ChemDist = chem.dist.trim$ChemDist

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation
cor(data$m, data$ChemDist, method="spearman") #spearman's correlation
cor(data$ChemDist, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
cryo.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
cryo.test.data$Habitat = "Cryoconite"

#put all our dataframes together
pro.all.habs = rbind(snow.test.data, sp.test.data, sum.test.data, cryo.test.data)

#STEP ONE ORDER THE LEVELS TO BE TESTED
pro.all.habs$Group <- ordered(pro.all.habs$Group,
                              levels = c(">55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                         "1000-1250m", "1250-1500m", "1500-2200m"))

pro.all.habs$Habitat <- ordered(pro.all.habs$Habitat,
                                levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#####PLOT THE DATA######
pro.p<-ggplot(pro.all.habs, aes(x=Group, y=Beta, fill=Habitat)) +
  geom_boxplot(position=position_dodge(1))+
  theme_bw()+
  theme(legend.position="bottom", legend.title = element_blank())+
  scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  xlab("Distance between sample pairs (m)")+
  ylab("Bray-Curtis Dissimilarity")
pro.p


# STEP TWO COMPUTE SUMMARY STATISTICS
pro.sum.stats = group_by(pro.all.habs, Group, Habitat) %>%
  summarise(
    count = n(),
    Beta_mean = mean(Beta, na.rm = TRUE),
    Beta_sd = sd(Beta, na.rm = TRUE),
    Chem_mean = mean(ChemDist, na.rm = TRUE),
    Chem_sd = sd(ChemDist, na.rm = TRUE))


####A DIFFERENT PLOT
pro.snow.df = pro.sum.stats[pro.sum.stats$Habitat == "Snow",]

scale=20
p1 = ggplot(pro.snow.df, aes(x=Group, y=Beta_mean)) +
  geom_line(aes(group = 1), color="#fa9f99") +
  geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#fa9f99") +
  geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
  geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#fa9f99")+
  geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Snow")
p1


pro.sp.df = pro.sum.stats[pro.sum.stats$Habitat == "Spring Ice",]
p2 = ggplot(pro.sp.df, aes(x=Group, y=Beta_mean)) +
  geom_line(aes(group = 1), color="#a4c64d") +
  geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#a4c64d") +
  geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
  geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#a4c64d")+
  geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Spring Ice")
p2


pro.sum.df = pro.sum.stats[pro.sum.stats$Habitat == "Summer Ice",]
p3 = ggplot(pro.sum.df, aes(x=Group, y=Beta_mean)) +
  geom_line(aes(group = 1), color="#4dd2d6") +
  geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#4dd2d6") +
  geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
  geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#4dd2d6")+
  geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Summer Ice")
p3

pro.cryo.df = pro.sum.stats[pro.sum.stats$Habitat == "Cryoconite",]
p4 = ggplot(pro.cryo.df, aes(x=Group, y=Beta_mean)) +
  geom_line(aes(group = 1), color="#d8a4ff") +
  geom_line(aes(group=1, y = Chem_mean/scale), color="#000000") +
  scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Environmental Euclidean Distance")) +
  expand_limits(y = 0)+
  theme_bw()+
  geom_errorbar( aes(ymin = Beta_mean-Beta_sd, ymax = Beta_mean+Beta_sd, color=Beta_mean),width = 0.2, color="#d8a4ff") +
  geom_errorbar( aes(ymin = (Chem_mean/scale)-(Chem_sd/scale), ymax = (Chem_mean/scale)+(Chem_sd/scale)),width = 0.2, color="#000000") +
  geom_point(aes(x=Group, y=Beta_mean), size = 2, color="#d8a4ff")+
  geom_point(aes(x=Group, y=Chem_mean/scale), size = 2, color="black")+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylab("Community Bray-Curtis Dissimilarity")+
  ggtitle("Cryoconite")
p4

pro.line.p = p1 + p2 + p3 + p4
pro.line.p

pdf("../results/prokaryote-line-distance.pdf")
print(pro.line.p)
dev.off()


#Eukaryotes

#####################euk SNOW#####################
ps = ps.euk
sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
snow.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
snow.test.data$Habitat = "Snow"



#####################euk SPRING ICE#####################
ps = ps.euk
sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sp.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sp.test.data$Habitat = "Spring Ice"


#####################euk SUMMER ICE#####################
ps = ps.euk
sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sum.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sum.test.data$Habitat = "Summer Ice"

#####################euk CRYO#####################

ps = ps.euk
sam = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
cryo.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
cryo.test.data$Habitat = "Cryoconite"

#put all our dataframes together
euk.all.habs = rbind(snow.test.data, sp.test.data, sum.test.data, cryo.test.data)

#STEP ONE ORDER THE LEVELS TO BE TESTED
euk.all.habs$Group <- ordered(euk.all.habs$Group,
                              levels = c(">55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                         "1000-1250m", "1250-1500m", "1500-2200m"))

euk.all.habs$Habitat <- ordered(euk.all.habs$Habitat,
                                levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#####PLOT THE DATA######
euk.p<-ggplot(euk.all.habs, aes(x=Group, y=Beta, fill=Habitat)) +
  geom_boxplot(position=position_dodge(1))+
  theme_bw()+
  theme(legend.position="bottom", legend.title = element_blank())+
  scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  xlab("Distance between sample pairs (m)")+
  ylab("Bray-Curtis Dissimilarity")
euk.p




# STEP TWO COMPUTE SUMMARY STATISTICS
euk.sum.stats = group_by(euk.all.habs, Group, Habitat) %>%
  summarise(
    count = n(),
    mean = mean(Beta, na.rm = TRUE),
    sd = sd(Beta, na.rm = TRUE)
  )

####A DIFFERENT PLOT
euk.snow.df = euk.sum.stats[euk.sum.stats$Habitat == "Snow",]
p1 = ggplot(euk.snow.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color=c("#fa9f99"))) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color=c("#fa9f99")),width = 0.2) +
  geom_point(size = 2, color=c("#fa9f99"))+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylim(0,1)+
  scale_color_manual(values = c("#fa9f99"))+
  ylab("Bray-Curtis Dissimilarity")+
  ggtitle("Snow")

euk.sp.df = euk.sum.stats[euk.sum.stats$Habitat == "Spring Ice",]
p2 = ggplot(euk.sp.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#a4c64d")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#a4c64d"), width = 0.2) +
  geom_point(size = 2, color="#a4c64d")+
  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank())+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  scale_color_manual(values = c("#a4c64d"))+
  ggtitle("Spring Ice")


euk.sum.df = euk.sum.stats[euk.sum.stats$Habitat == "Summer Ice",]
p3 = ggplot(euk.sum.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#4dd2d6")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#4dd2d6"), width = 0.2) +
  geom_point(size = 2, color="#4dd2d6")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  scale_color_manual(values = c("#4dd2d6"))+
  xlab("Distance between sample pairs (m)")+
  ggtitle("Summer Ice")

euk.cryo.df = euk.sum.stats[euk.sum.stats$Habitat == "Cryoconite",]
p4 = ggplot(euk.cryo.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#d8a4ff")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#d8a4ff"), width = 0.2) +
  geom_point(size = 2, color="#d8a4ff")+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  xlab("Distance between sample pairs (m)")+
  scale_color_manual(values = c("#d8a4ff"))+
  ggtitle("Cryoconite")

#plot_grid(p1, p2, p3, p4)

euk.line.p = p1 + p2 + p3 + p4

pdf("../results/eukkaryote-line-distance.pdf")
print(euk.line.p)
dev.off()

#Microfauna

#####################mm SNOW#####################
ps = ps.mm
sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
snow.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
snow.test.data$Habitat = "Snow"



#####################mm SPRING ICE#####################
ps = ps.mm
sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sp.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sp.test.data$Habitat = "Spring Ice"


#####################mm SUMMER ICE#####################
ps = ps.mm
sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
sum.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
sum.test.data$Habitat = "Summer Ice"

#####################mm CRYO#####################

ps = ps.mm
sam = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID

#Separate habitat data
h <- prune_samples(sam, ps)

#remove empty samples or ASVs
h.f = phyloseq::filter_taxa(h, function(x) sum(x) > 0, TRUE)
h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)

#Calculate Bray-Curtis dissimilarities between sites
bc <- vegdist(as.matrix(t(otu_table(h.f))), method="bray")

#Calculate geographic distance between sites
ll <- data.frame(sample_data(h.f ))[c("Longitude", "Latitude")]
geo.dist <- earth.dist(ll, dist=TRUE) #distance in kilometers

#transform these data into columns
comm.dist = melt(as.matrix(bc))
names(comm.dist) = c("Sample1", "Sample2", "Beta")
comm.dist.trim  = comm.dist[comm.dist['Sample1'] != comm.dist['Sample2'],]

geo.dist = melt(as.matrix(geo.dist))
names(geo.dist) = c("Sample1", "Sample2", "km")
geo.dist.trim = geo.dist[geo.dist['Sample1'] != geo.dist['Sample2'],]

#make into dataframe
data = comm.dist.trim 
data$km = geo.dist.trim$km
data$m = data$km*1000

plot(data$m, data$Beta) #we can see a somewhat linear relationship
abline(v=50, col="blue")
abline(v=200, col="blue")
abline(v=500, col="blue")
abline(v=1000, col="blue")

cor(data$m, data$Beta, method="spearman") #spearman's correlation

####USING A ONE WAY ANOVA TO TEST DIFFERENCES IN THE MEAN BRAY-CURTIS VALUES OF SAMPLES 0-50M APART AND 50-200M APART

#is there a significant difference in the bray-curtis values between 0-50m and 50-200m?
df_55 = data[data$m <= 55,]
df_55$Group = ">55m"
df_200 = data[data$m > 55 & data$m <= 200,]
df_200$Group = "55-200m"
df_300 = data[data$m > 200 & data$m <= 300,]
df_300$Group = "200-300m"
df_400 = data[data$m > 300 & data$m <= 400,]
df_400$Group = "300-400m"
df_500 = data[data$m > 400 & data$m <= 500,]
df_500$Group = "400-500m"
df_750 = data[data$m > 500 & data$m <= 750,]
df_750$Group = "500-750m"
df_1000 = data[data$m > 750 & data$m <= 1000,]
df_1000$Group = "750-1000m"
df_1250 = data[data$m > 1000 & data$m <= 1250,]
df_1250$Group = "1000-1250m"
df_1500 = data[data$m > 1250 & data$m <= 1500,]
df_1500$Group = "1250-1500m"
df_2200 = data[data$m > 1500 & data$m <= 2200,]
df_2200$Group = "1500-2200m"

#DATA FOR TESTING AND PLOTTING
cryo.test.data = rbind(df_55, df_200, df_300, df_400, df_500, df_750, df_1000, df_1250, df_1500, df_2200)
cryo.test.data$Habitat = "Cryoconite"

#put all our dataframes together
mm.all.habs = rbind(snow.test.data, sp.test.data, sum.test.data, cryo.test.data)

#STEP ONE ORDER THE LEVELS TO BE TESTED
mm.all.habs$Group <- ordered(mm.all.habs$Group,
                             levels = c(">55m", "55-200m", "200-300m", "300-400m", "400-500m", "500-750m", "750-1000m",
                                        "1000-1250m", "1250-1500m", "1500-2200m"))

mm.all.habs$Habitat <- ordered(mm.all.habs$Habitat,
                               levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

#####PLOT THE DATA######
mm.p<-ggplot(mm.all.habs, aes(x=Group, y=Beta, fill=Habitat)) +
  geom_boxplot(position=position_dodge(1))+
  theme_bw()+
  theme(legend.position="bottom", legend.title = element_blank())+
  scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  xlab("Distance between sample pairs (m)")+
  ylab("Bray-Curtis Dissimilarity")
mm.p




# STEP TWO COMPUTE SUMMARY STATISTICS
mm.sum.stats = group_by(mm.all.habs, Group, Habitat) %>%
  summarise(
    count = n(),
    mean = mean(Beta, na.rm = TRUE),
    sd = sd(Beta, na.rm = TRUE)
  )

####A DIFFERENT PLOT
mm.snow.df = mm.sum.stats[mm.sum.stats$Habitat == "Snow",]
p1 = ggplot(mm.snow.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color=c("#fa9f99"))) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color=c("#fa9f99")),width = 0.2) +
  geom_point(size = 2, color=c("#fa9f99"))+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x=element_blank())+
  ylim(0,1)+
  scale_color_manual(values = c("#fa9f99"))+
  ylab("Bray-Curtis Dissimilarity")+
  ggtitle("Snow")

mm.sp.df = mm.sum.stats[mm.sum.stats$Habitat == "Spring Ice",]
p2 = ggplot(mm.sp.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#a4c64d")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#a4c64d"), width = 0.2) +
  geom_point(size = 2, color="#a4c64d")+
  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), axis.ticks=element_blank())+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  scale_color_manual(values = c("#a4c64d"))+
  ggtitle("Spring Ice")


mm.sum.df = mm.sum.stats[mm.sum.stats$Habitat == "Summer Ice",]
p3 = ggplot(mm.sum.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#4dd2d6")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#4dd2d6"), width = 0.2) +
  geom_point(size = 2, color="#4dd2d6")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  scale_color_manual(values = c("#4dd2d6"))+
  xlab("Distance between sample pairs (m)")+
  ggtitle("Summer Ice")

mm.cryo.df = mm.sum.stats[mm.sum.stats$Habitat == "Cryoconite",]
p4 = ggplot(mm.cryo.df, aes(x=Group, y=mean)) +
  geom_line(aes(group = 1, color="#d8a4ff")) +
  theme_bw()+
  geom_errorbar( aes(ymin = mean-sd, ymax = mean+sd, color="#d8a4ff"), width = 0.2) +
  geom_point(size = 2, color="#d8a4ff")+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylim(0,1)+
  ylab("Bray-Curtis Dissimilarity")+
  xlab("Distance between sample pairs (m)")+
  scale_color_manual(values = c("#d8a4ff"))+
  ggtitle("Cryoconite")

#plot_grid(p1, p2, p3, p4)

mm.line.p = p1 + p2 + p3 + p4

pdf("../results-proportional/mmkaryote-line-distance.pdf")
print(mm.line.p)
dev.off()

#Combine Plots + save results tables


all.p = plot_grid(pro.p + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()), 
                  euk.p+ theme(legend.position="none", axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()), mm.p,
                  ncol=1, rel_heights = c(0.7,0.7,1.2))

pdf("../results/boxplots-distances.pdf", width=35, height=20)
print(all.p)
dev.off()


write.csv(pro.sum.stats, "../results/prokaryote-distance-stats.csv")
write.csv(euk.sum.stats, "../results/eukaryote-distance-stats.csv")
write.csv(mm.sum.stats, "../results/microfauna-distance-stats.csv")



pro.line.p | euk.line.p | mm.line.p


thm <- theme(plot.title = element_text(face = 2, size = 14))
first_plot      <- wrap_elements(pro.line.p + plot_annotation(title = "Prokaryote", theme = thm))
second_plot   <- wrap_elements(euk.line.p + plot_annotation(title = "Microbial eukaryote", theme = thm))
third_plot   <- wrap_elements(mm.line.p + plot_annotation(title = "Microfauna", theme = thm))

all.line.p = first_plot | second_plot | third_plot
all.line.p

pdf("../results/all-dissimilarity-line-plots.pdf", width=25, height=10)
print(all.line.p)
dev.off()