# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

source("00-solman-functions.R")
library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist, cca and mantel functions
library(dplyr) #for %>% and full_join functions
library(ggcorrplot)
library(psych)
library(patchwork)

#prokaryotes
ps.pro <- readRDS("../results/16S-ps-norm.rds")
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-no-mm-ps-norm.rds")
#micrometazoans only
ps.mm <- readRDS("../results/18S-mm-only-ps-norm.rds")

#Mantal test function

# ps = ps.pro
# hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID
# habitat = "Cryoconite"

mantel_single_var <- function(ps, hab, habitat){
  
  #subset data by habitat
  ps.hab <- prune_samples(hab, ps)
  
  #get community data
  comm <- data.frame(t(otu_table(ps.hab)), check.names=FALSE)
  
  #get metadata
  meta = data.frame(sample_data(ps.hab))
  meta = meta[,c(9:ncol(meta))]
  #remove columns with all NA values
  meta.trim <- meta[ , colSums(is.na(meta)) < nrow(meta)] 
  #remove columns with less than 10 values > 0
  # c <- meta
  # c[c > 0 | c < 0] <- 1
  # c[is.na(c)] <- 0
  # meta.trim <- meta[,colSums(c) >=10]
  
  #save output
  variname = vector()
  rho = vector()
  pval = vector()
  
  for (i in 1:ncol(meta.trim)){
    
    # Get Bray-Curtis dissimilarity matrix for community data
    commdist<-vegdist(comm, method="bray")
    
    vardist<-vegdist(meta.trim[,i], method="euclidean", na.rm = TRUE)
    
    set.seed(666)
    
    comm_env <- vegan::mantel(commdist, vardist, method="spearman", permutations=999, na.rm=TRUE)
    comm_env
    
    variname = c(variname, names(meta.trim)[[i]])
    rho = c(rho, round(comm_env$statistic, 2))
    pval = c(pval, round(comm_env$signif, 4))
    
  }
  
  res.df = data.frame(Variable = variname, Rho = rho, Pval = pval)
  names(res.df) = c("Variable", paste0(habitat, " MantelR"), paste0(habitat, " P val"))
  
  return(res.df)
  
  
}


#Run the function

pro.snow <- mantel_single_var(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID, habitat="Snow")

pro.sp = mantel_single_var(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Spring Ice")))$SampleID, habitat="Spring Ice")

pro.sum = mantel_single_var(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Summer Ice")))$SampleID, habitat="Summer Ice")

pro.cry = mantel_single_var(ps = ps.pro, hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite")

euk.snow = mantel_single_var(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Snow")))$SampleID, habitat="Snow")

euk.sp = mantel_single_var(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Spring Ice")))$SampleID, habitat="Spring Ice")

euk.sum = mantel_single_var(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Summer Ice")))$SampleID, habitat="Summer Ice")

euk.cry = mantel_single_var(ps = ps.euk, hab = data.frame(sample_data(subset_samples(ps.euk, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite")

mm.snow = mantel_single_var(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID, habitat="Snow")

mm.sp = mantel_single_var(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID, habitat="Spring Ice")

mm.sum = mantel_single_var(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID, habitat="Summer Ice")

mm.cry = mantel_single_var(ps = ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID, habitat="Cryoconite")


#bind results into dataframes
pro.df = full_join(pro.snow, pro.sp, by="Variable")
pro.df = full_join(pro.df, pro.sum, by="Variable")
pro.df = full_join(pro.df, pro.cry, by="Variable")

euk.df = full_join(euk.snow, euk.sp, by="Variable")
euk.df = full_join(euk.df, euk.sum, by="Variable")
euk.df = full_join(euk.df, euk.cry, by="Variable")

mm.df = full_join(mm.snow, mm.sp, by="Variable")
mm.df = full_join(mm.df, mm.sum, by="Variable")
mm.df = full_join(mm.df, mm.cry, by="Variable")

#put all dataframes together and reorder for export
full.df = full_join(pro.df, euk.df, by="Variable")
full.df = full_join(full.df, mm.df, by="Variable")
vars.ord = c("pH", "Conductivity_muS", "Temp", "Area", "Depth_cm",
             "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
             "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl",
             "Li", "Al", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu",
             "Zn", "Rb", "Sr", "Mo", "Ag", "Cd", "Ba", "Lu", 
             "Pb", "Zr", "Nb", "Sn", "Cs", "Ta", "Re", "U",
             "Y", "Hf", "La", "Ce", "Pr", "Nd", "Sm", "Eu",
             "Gd", "Tb", "Dy", "Ho", "Er", "Yb")

#set the correct order
full.df.ord = full.df[match(vars.ord, full.df$Variable),]
write.csv(full.df.ord, "../results/mantel-single-var-results.csv")

df = full_join(pro.snow, pro.sp, by="Variable")
df = full_join(df, pro.sum, by="Variable")
df = full_join(df, pro.cry, by="Variable")
rownames(df) = df$Variable

#subset to major physical factors, nutrients and ions only
keep = c("pH", "Conductivity_muS", "Area", "Depth_cm",
         "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")
df.sub = df[rownames(df) %in% keep,]

#set the correct order
df.sub.ord = df.sub[match(keep, rownames(df.sub)),]

#rename rows
rownames(df.sub.ord) = c("pH", "Conductivity", "Area", "Depth",
                         "TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]", "SO[4]", 
                         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")

#remove some rows because they don't have any significant correlations
df.sub.ord = df.sub.ord[! rownames(df.sub.ord) %in% c("TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]"),]

r.df = df.sub.ord[,c(8,6,4,2)]
names(r.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")
p.df = df.sub.ord[,c(9,7,5,3)]
names(p.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")

pro.p = ggcorrplot::ggcorrplot(r.df, method = "square", lab=TRUE, sig.level = 0.05, p.mat = p.df, ggtheme = ggplot2::theme_bw, 
                               show.legend = TRUE,
                               insig="blank", lab_size = 4,
                               colors=c("#88d2d6", "#FFFFFF", "#f8c4be"))+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(face="bold"))+
  ggtitle("Prokaryote") +
  scale_x_discrete(labels = scales::parse_format())
pro.p

pdf("../results/mantel-correlation-tests-prokaryotes.pdf", height=10, width=10)
print(pro.p)
dev.off()



#Eukaryotes

df = full_join(euk.snow, euk.sp, by="Variable")
df = full_join(df, euk.sum, by="Variable")
df = full_join(df, euk.cry, by="Variable")
rownames(df) = df$Variable

#subset to major nutrients and ions only
keep = c("pH", "Conductivity_muS", "Area", "Depth_cm",
         "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")
df.sub = df[rownames(df) %in% keep,]
#set the correct order
df.sub.ord = df.sub[match(keep, rownames(df.sub)),]
#rename rows
rownames(df.sub.ord) = c("pH", "Conductivity", "Area", "Depth",
                         "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
                         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")

#remove some rows because they don't have any significant correlations
df.sub.ord = df.sub.ord[! rownames(df.sub.ord) %in% c("TC", "DOC", "TN", "NO2", "NO3", "PO4"),]

r.df = df.sub.ord[,c(8,6,4,2)]
names(r.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")
p.df = df.sub.ord[,c(9,7,5,3)]
names(p.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")

euk.p = ggcorrplot::ggcorrplot(r.df, method = "square", lab=TRUE, sig.level = 0.05, p.mat = p.df, ggtheme = ggplot2::theme_bw, show.legend = TRUE,
                               insig="blank", lab_size = 4,
                               colors=c("#88d2d6", "#FFFFFF", "#f8c4be"))+
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x=element_blank(), plot.title = element_text(face="bold"))+
  ggtitle("Microbial Eukaryote")
euk.p

pdf("../results/mantel-correlation-tests-eukaryotes.pdf", height=10, width=10)
print(euk.p)
dev.off()


#Microfauna

df = full_join(mm.snow, mm.sp, by="Variable")
df = full_join(df, mm.sum, by="Variable")
df = full_join(df, mm.cry, by="Variable")
rownames(df) = df$Variable

#subset to major nutrients and ions only
keep = c("pH", "Conductivity_muS", "Area", "Depth_cm",
         "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")
df.sub = df[rownames(df) %in% keep,]
#set the correct order
df.sub.ord = df.sub[match(keep, rownames(df.sub)),]
#rename rows
rownames(df.sub.ord) = c("pH", "Conductivity", "Area", "Depth",
                         "TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]", "SO[4]", 
                         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")

#remove some rows because they don't have any significant correlations
df.sub.ord = df.sub.ord[! rownames(df.sub.ord) %in% c("TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]"),]

r.df = df.sub.ord[,c(8,6,4,2)]
names(r.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")
p.df = df.sub.ord[,c(9,7,5,3)]
names(p.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")

mm.p = ggcorrplot::ggcorrplot(r.df, method = "square", lab=TRUE, sig.level = 0.05, p.mat = p.df, ggtheme = ggplot2::theme_bw, show.legend = TRUE,
                              insig="blank", lab_size = 4,
                              colors=c("#88d2d6", "#FFFFFF", "#f8c4be"))+
  theme(legend.position = "none", plot.title = element_text(face="bold"))+
  ggtitle("Microfauna")+
  scale_x_discrete(labels = scales::parse_format())
mm.p

pdf("../results/mantel-correlation-tests-microfauna.pdf", height=10, width=10)
print(mm.p)
dev.off()


full.p = pro.p / euk.p / mm.p
full.p

pdf("../results/mantel-correlation-tests-all.pdf", height=8, width=8)
print(full.p)
dev.off()

#GET PLOT WITHOUT MICROFAUNA
df = full_join(euk.snow, euk.sp, by="Variable")
df = full_join(df, euk.sum, by="Variable")
df = full_join(df, euk.cry, by="Variable")
rownames(df) = df$Variable

#subset to major nutrients and ions only
keep = c("pH", "Conductivity_muS", "Area", "Depth_cm",
         "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", 
         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")
df.sub = df[rownames(df) %in% keep,]
#set the correct order
df.sub.ord = df.sub[match(keep, rownames(df.sub)),]
#rename rows
rownames(df.sub.ord) = c("pH", "Conductivity", "Area", "Depth",
                         "TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]", "SO[4]", 
                         "Cl", "Fe", "Na", "Mg", "K", "Ca", "Br", "Fl")

#remove some rows because they don't have any significant correlations
df.sub.ord = df.sub.ord[! rownames(df.sub.ord) %in% c("TC", "DOC", "TN", "NO[2]", "NO[3]", "PO[4]", "Fl"),]

r.df = df.sub.ord[,c(8,6,4,2)]
names(r.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")
p.df = df.sub.ord[,c(9,7,5,3)]
names(p.df) = c("Cryoconite", "Summer Ice", "Spring Ice", "Snow")

full.p2 = pro.p / (ggcorrplot::ggcorrplot(r.df, method = "square", lab=TRUE, sig.level = 0.05, p.mat = p.df, ggtheme = ggplot2::theme_bw, show.legend = TRUE,
                                          insig="blank", lab_size = 4,
                                          colors=c("#88d2d6", "#FFFFFF", "#f8c4be"))+
                     theme(legend.position = "none", plot.title = element_text(face="bold"))+
                     ggtitle("Microbial Eukaryote")+
                     scale_x_discrete(labels = scales::parse_format()))
 
full.p2

pdf("../results/mantel-correlation-tests-pro-euk.pdf", height=6, width=8)
print(full.p2)
dev.off()