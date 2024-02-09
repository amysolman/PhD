#Differential Abundance Analysis with Wilcoxon rank sum, DESeq2 and Aldex2

rm(list=ls())
library(dplyr)
library(ggpubr)
library(car) #for levin's test of homogeneity of variances
library(broom)
library(tidyr)
library(stringr) #to split the column
source("00-solman-functions.R")
library(cowplot)
library(knitr)
library(phyloseq)
library(vegan)
# BiocManager::install("DESeq2")
library(DESeq2)
# require(plyr) #for mean and standard error
library(dplyr)
# BiocManager::install("ALDEx2")
library(ALDEx2)
library(tibble)
library(patchwork)
library(plyr)

#first we want to use our count data so re-import it from before it was proportionally transfromed

#eukkaryotes
ps.pro <- readRDS("../results/16S-ps-controls-removed.rds")
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-ps-no-mm-controls-removed.rds")

#remove the sample low read count samples as we did previously

#prokaryotes
ps.p.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro) #retain samples with >= num counts
ps.f.pro = filter_taxa(ps.p.pro, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts
#eukaryotes
ps.p.euk = prune_samples(sample_sums(ps.euk)>= 1900, ps.euk) #retain samples with >= num counts
ps.f.euk = filter_taxa(ps.p.euk, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts



#Wilcoxon Rank Sum Test

# phylo = ps.f.pro
# group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,
#           data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID)
# lev = "Family"
# hab1 = "Spring Ice"
# hab2 = "Cryoconite"

wilcoxon.daa <- function(phylo, group, lev){
  
  #Separate habitat data
  ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(ps, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  
  #clr transformation
  ps.clr <- microbiome::transform(ps.trim, "clr")
  
  #extract the data
  aa_data <- data.frame(t(otu_table(ps.clr)), check.names = FALSE)
  #counts <- data.frame(t(otu_table(ps.trim)))
  
  #add info to the data
  aa_data <- cbind(aa_data, 
                   Habitat = data.frame(sample_data(ps.clr))$Habitat)
  
  agg <- names(aa_data[, !names(aa_data) %in% "Habitat"])
  
  wilcoxon_p <- c() # Initialize empty vector for p-values
  foldchange <- c() # Initiate empty vector for effect size (difference between log means) - see https://microbiome.github.io/tutorials/all.html
  
  for (i in agg) {
    
    # test to see if the distribution of clr transformed abundances is different between the two habitats
    result <- wilcox.test(aa_data[, i] ~ Habitat,
                          data = aa_data) 
    
    # get the z (slope) coefficient of a linear model of clr transformed abundances between the two habitats (strength of the relationship)
    e <- coef(lm(aa_data[,i] ~ Habitat, data = aa_data))
    
    # Stores p-value to the vector with this column name
    wilcoxon_p[[i]]  <- result$p.value
    foldchange[[i]] <- e[[2]]
    
  }
  
  wilcoxon_p <- data.frame(taxa =  names(wilcoxon_p),
                           p_raw = unlist(wilcoxon_p))
  
  foldchange <- data.frame(taxa =  names(foldchange),
                           effect = unlist(foldchange))
  test.res = full_join(wilcoxon_p, foldchange, by = "taxa")
  
  #adjust p-values
  test.res$p_adjusted <- p.adjust(test.res$p_raw, method = "fdr")
  
  #only keep rows with significant p-values
  #and absolute effect sizes (linear rgeression slopes) > 1
  test.res_sig <- test.res[which(test.res$p_adjusted < 0.05 & abs(test.res$effect) > 1),]
  
  #add taxonomic info
  tax.info = data.frame(tax_table(ps.clr))
  tax.info$taxa = rownames(tax.info)
  test.res_sig_tax = left_join(test.res_sig, tax.info, by="taxa")
  
  #remove those that say uncultured or of unknown placement
  test.res_sig_tax_rm = test.res_sig_tax[test.res_sig_tax$Family != "uncultured",]
  test.res_sig_tax_rm = test.res_sig_tax[test.res_sig_tax$Family != "Incertae_Sedis",]
  
  return(test.res_sig_tax_rm)
  
}

# DESeq2
# http://127.0.0.1:19618/library/DESeq2/doc/DESeq2.html

# phylo = ps.f.pro
# # group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,
# #           data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
# lev = "Family"
# # hab1 = "Snow"
# # hab2 = "Spring Ice"

DESeq2 <- function(phylo, lev){
  
  #Separate habitat data
  #ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(phylo, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  
  #convert phylo object to deseq2 object
  #no need to normalise the data
  dds = phyloseq_to_deseq2(ps.trim, ~ Habitat)
  dds #look at features of our deseq2 object 
  
  #side note - look at the distribution of our count data
  #we can see that it is NOT normally distributed 
  #therefore we cannot use methods for normally distributed data to predict the likelihood
  #of events or features like the mean/sd etc
  #poisson distributions are more appropriate to fit for microbiome/RNA-Seq data
  
  # hist(rowSums(data.frame(otu_table(ps.trim))))
  
  #because our count data is sparce DESeq has issues calculating the geometric mean 
  #so we need to calculate the size factors using a zero-tolerant method
  #calculate geometric means prior to estimate size factors 
  #(https://github.com/joey711/phyloseq/issues/387)
  
  # gm_mean = function(x, na.rm=TRUE){
  # exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  # }
  # 
  # geoMeans = apply(counts(dds), 1, gm_mean)
  
  # dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  
  #carry out test
  #the DESeq2 differential expression analysis steps are wrapped into a single function DESEq.
  #we have performed this on all four conditions (Habitats) and from this we can extract results for all comparisosns
  dds = DESeq(dds) #comparison is based on alphabetic order of the condition
  
  #so for each ASV, we can ask DESeq to compare two conditions by:
  #estimating the family-wise dispersions
  #using those values to fit a negative binomial generalized linear model 
   #- a model to look at the relationship of abundance between two conditions
  #performing a Wald significance test on that model
  #and give us the following information:
  #baseMean: the average of the normalised counts taken over all samples
  #log2FoldChange: the log2 fold change between the two groups, e.g. log2 fold change of 2 = a 4-fold increase
  #IfcSE: standard error of the log2FoldChange estimate
  #stat: Wald statistic 
  #p-value: Wald test p-value
  #padj = adjusted p-value of the Wald test
  
  #investigate test results
  #extract results using results function
  #first condition = numerator for fold change
  #second condition = denominator for form change
  res1 = results(dds, contrast = c("Habitat", "Snow", "Spring.Ice"), cooksCutoff = FALSE)
  res1 #info about results
  look.res1 = data.frame(res1)
  res2 = results(dds, contrast = c("Habitat", "Snow", "Summer.Ice"), cooksCutoff = FALSE)
  res2 #info about results
  res3 = results(dds, contrast = c("Habitat", "Snow", "Cryoconite"), cooksCutoff = FALSE)
  res3 #info about results
  res4 = results(dds, contrast = c("Habitat", "Spring.Ice", "Summer.Ice"), cooksCutoff = FALSE)
  res4 #info about results
  res5 = results(dds, contrast = c("Habitat", "Spring.Ice", "Cryoconite"), cooksCutoff = FALSE)
  res5 #info about results
  res6 = results(dds, contrast = c("Habitat", "Summer.Ice", "Cryoconite"), cooksCutoff = FALSE)
  res6 #info about results
  
  #results into a list
  res.list = list(res1, res2, res3, res4, res5, res6)
  
  #only keep results which are statistically and practically significant
  alpha = 0.05 #p-value (significance cutoff)
  effect.size = 1
  test.nam = c("Snow Vs Spring Ice", "Snow Vs Summer Ice", "Snow Vs Cryoconite", 
               "Spring Ice Vs Summer Ice", "Spring Ice Vs Cryoconite",
               "Summer Ice Vs Cryoconite")
  
  for (i in 1:length(res.list)){
    
    res = res.list[[i]]
    sigtab = res[which(res$padj < alpha & abs(res$log2FoldChange) > effect.size), ]
    sigtab = cbind("Test" = test.nam[[i]], as(sigtab, "data.frame"), as(tax_table(ps.trim)[rownames(sigtab), ], "matrix"))
    
    res.list[[i]] = sigtab
  }
  
  
  #put all our results into one table
  data = bind_rows(res.list)
  
  #remove NA data
  data = data[!is.na(data$Family) ,]
  
  #remove uncultured data
  data = data[data$Family != "uncultured",]
  data = data[data$Family != "Incertae_Sedis",]
  
  return(data)
  
}

# Aldex2
# Estimates within sample variation using the Dirichlet distribution.
# Applied centrered-log-ratio transformation.
# If comparing two conditions Aldex2 uses Welch's t-test and Wilcoxon tests.
# If comparing three or more conditions it will perform a one-way ANOVA and Kruskall-Wallis test.
# The Benjamini-Hochberg p-value adjustment is used by default.
# https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html
# phylo = ps.f.pro
# group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,
#           data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID)
# lev = "Family"
# hab1 = "Snow"
# hab2 = "Spring Ice"

 ALDEx2 <- function(phylo, group, lev){
  
  set.seed(666)
   
  #Separate habitat data
  ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(ps, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  #ps.trim <- filter_taxa(ps, function(x) sum(x) > 0, TRUE)
  
  #extract count table of un-normalised data with samples as columns
  counts <- data.frame(otu_table(ps.trim), check.names = FALSE)
  
  #get vector of conditions in the same order as the samples in the counts table
  conds <- c(data.frame(sample_data(ps.trim))$Habitat)

  #aldex is a wrapper function linking the modular elements of the aldex2 package. 
  #aldex function performs a two-sample t-test and calculates effect sizes. 
  #if test="t" then effect should = TRUE
  #"t" option evaluates the data as a two-factor experiment using Welch's t-test and Wilcoxon rank tests 
#Welch's t test = null hypothesis two population have equal means (the populations are assumed to be normal distributed but don't need to have equal variances)
#Wilcoxon rank tests = null hypothesis two populations have equal means (there are no assumptions about the distribution of the data)
#For multiple sample tests using ANOVA (to test difference in means of multiple groups) effect should = FALSE was effect size should not be calculated.
#Effect size = value measuring strength of the relationship between two variables in a population. This is different from p-value which tells you the effect exists. p-value = statistical significance (there is a difference). effect size = practical significance (the difference is big enough to care about).  
#test = "kw" evaluate the data as a one-way ANOVA (Analysis of Variance - how does ONE independent factor, e.g. habitat, impact the dependent variable, e.g., CLR values. A two-way ANOVA evaluated the impact of TWO independent variables on a dependent variable) using the glm and Kruskal-Wallace tests (tests for significant differences between the medians of three or more independent groups consisting of non-parametric data).
#All tests include a Benjamini-Heckberg correction of raw p-values. 
#data can be plotted into Bland-Altman (MA) or effect (MW) plots for two0way tests. 
#Bland-Altmann plots are difference plots. There are used to evaluate the agreement between two measurements.
#as we want to compare pairs of habitats here we will use the following arguments:
#test="t" - compare two conditions
#effect=TRUE - calculate the effect size of these differences
#mc.samples = 1000 - the number of Dirichlet Monte-Carlo Instances/Samples to use when estimating the underlying distributions. These are Monte-Carlo method random sampling of the Dirichlet probability distribution (to estimate the chance of getting our result at random).
#include.sample.summary = FALSE - whether to include median clr values for each sample
#verbose = TRUE - print diagnostic information while running functions
#denom = "all" - which features to retain as the denominator for the Geometric Mean calculation
#paired.test = FALSE - are the samples paired
  
x.all <- aldex(counts, conds, mc.samples=1000, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=TRUE, paired.test=FALSE)

#insignficant features (families) are shown in grey or black. Statistically significant features (those families which are significant different between habitats) are shown in red
#log-ratio abundance axis is the clr value
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")

#get only our significant values and those with effect size > 1
#we want to look at wi.eBH = expected BH corrected p value of Wilcoxon test 
res.trim = rownames_to_column(x.all, "ID") %>%
  filter(wi.eBH <= 0.05, abs(effect) > 1) # here we chose the wilcoxon output rather than ttest #keep the family ID, wilcoxon BH adj p-value, dif.btw (difference between the group median for each significant feature), effect (effect size), overlap (proportion of effect size that overlaps 0))

#add taxonomic info
tax.info = data.frame(tax_table(ps.trim))
tax.info$ID = rownames(tax.info)
res.trim_tax = left_join(res.trim, tax.info, by="ID")

return(res.trim_tax)

}

# Statistically significant p-val < 0.05, practically significant effect size > 1.
# 
# Combine Results/Plot
# Wilcoxon
# I have a list of Families with significant differences identified by Wilcoxon tests (does not say which is more or less - just different)
# 
# DESeq2
# Returns dataframe with statistically and practically significantly different families between habitats, which pair of habitats are being compared, the effect size (a.k.a. log2FoldChange) and the adjusted p-value.
# 
# Aldex2
# Median clr values for each significantly different family, the difference between the two, the effect size and the adjusted p-value.

#Prokaryotes

#Get results of each analysis and identify which Families overlap

#Wilcoxon
pro.wx1 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID), lev = "Family")
pro.wx1 = cbind(Test = "Snow Vs Spring Ice", pro.wx1)

pro.wx2 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID), lev = "Family")
pro.wx2 = cbind(Test = "Snow Vs Summer Ice", pro.wx2)

pro.wx3 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), lev = "Family")
pro.wx3 = cbind(Test = "Snow Vs Cryoconite", pro.wx3)

pro.wx4 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID), lev = "Family") #no differences
#pro.wx4 = cbind(Test = "Spring Ice Vs Summer Ice", pro.wx4) #no significant differences

pro.wx5 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), lev = "Family")
pro.wx5 = cbind(Test = "Spring Ice Vs Cryoconite", pro.wx5)

pro.wx6 = wilcoxon.daa(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), lev = "Family")
pro.wx6 = cbind(Test = "Summer Ice Vs Cryoconite", pro.wx6)

#DESeq2
pro.d2 = DESeq2(ps.f.pro, "Family")

#Aldex2
pro.ax1 = ALDEx2(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID), lev = "Family")
pro.ax1 = cbind(Test = "Snow Vs Spring Ice", pro.ax1)

pro.ax2 = ALDEx2(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID), lev = "Family")
pro.ax2 = cbind(Test = "Snow Vs Summer Ice", pro.ax2)

pro.ax3 = ALDEx2(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), lev = "Family")
pro.ax3 = cbind(Test = "Snow Vs Cryoconite", pro.ax3)

pro.ax4 = ALDEx2(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID), lev = "Family")
#pro.ax4 = cbind(Test = "Spring Ice Vs Summer Ice", pro.ax4) #no significant differences

pro.ax5 = ALDEx2(phylo = ps.f.pro, 
                 group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID,
                           data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), 
                 lev = "Family")

pro.ax5 = cbind(Test = "Spring Ice Vs Cryoconite", pro.ax5)

pro.ax6 = ALDEx2(phylo = ps.f.pro, group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID), lev = "Family")
pro.ax6 = cbind(Test = "Summer Ice Vs Cryoconite", pro.ax6)


#Which Families were found in all three analyses?

#Prokaryotes

#Snow and Spring Ice
pro.snow.sp.ice.tax = Reduce(intersect, list(pro.wx1$Family,
                                             pro.d2[pro.d2$Test == "Snow Vs Spring Ice",]$Family,
                                             pro.ax1$Family))

#Snow and Summer Ice
pro.snow.sum.ice.tax = Reduce(intersect, list(pro.wx2$Family,
                                              pro.d2[pro.d2$Test == "Snow Vs Summer Ice",]$Family,
                                              pro.ax2$Family))

#Snow and Cryoconite
pro.snow.cry.tax = Reduce(intersect, list(pro.wx3$Family,
                                          pro.d2[pro.d2$Test == "Snow Vs Cryoconite",]$Family,
                                          pro.ax3$Family))

#Spring Ice and Summer Ice
# pro.sp.sum.ice.tax = Reduce(intersect, list(pro.wx4$Family,
#                        pro.d2[pro.d2$Test == "Spring Ice Vs Summer Ice",]$Family,
#                        pro.ax4$Family))

#Spring Ice and Cryoconite
pro.sp.ice.cry.tax = Reduce(intersect, list(pro.wx5$Family,
                                            pro.d2[pro.d2$Test == "Spring Ice Vs Cryoconite",]$Family,
                                            pro.ax5$Family))

#Summer Ice and Cryoconite
pro.sum.ice.cry.tax = Reduce(intersect, list(pro.wx6$Family,
                                             pro.d2[pro.d2$Test == "Summer Ice Vs Cryoconite",]$Family,
                                             pro.ax6$Family))

#All families with significant differences
tot.fam = unique(c(pro.snow.sp.ice.tax, pro.snow.sum.ice.tax, pro.snow.cry.tax, pro.sp.ice.cry.tax, pro.sum.ice.cry.tax))

#get full tax data
tax.data = data.frame(tax_table(ps.pro))
tax.data = tax.data[tax.data$Family %in% tot.fam,]
tax.data.fam = unique(data.frame(Domain=tax.data$Domain, Phylum = tax.data$Phylum, Class=tax.data$Class, 
                                 Order=tax.data$Order,
                                 Family = tax.data$Family))
#sort by phylum
tax.data.fam = tax.data.fam[order(tax.data.fam$Phylum), ]

#print taxa with significant differences 
pro.snow.sp.ice.tax
pro.snow.sum.ice.tax
pro.snow.cry.tax
pro.sp.ice.cry.tax
pro.sum.ice.cry.tax

#Plot DESeq2 results only for those families with significant differences in all three analyses.


#SNOW AND SPRING ICE

data2plot1 = pro.d2[pro.d2$Test == "Snow Vs Spring Ice",]
data2plot1 = data2plot1[data2plot1$Family %in% pro.snow.sp.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot1$Family,]
rows2add2 = data.frame(Test = "Snow Vs Spring Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot1 = rbind(data2plot1, rows2add2)

data2plot1 =arrange(data2plot1 ,desc(Phylum),Family)
data2plot1$Family = paste0(data2plot1$Family, " (", data2plot1$Phylum, ")")

data2plot1$Family = factor(data2plot1$Family, levels = unique(data2plot1$Family))

#plot results
p1 = ggplot(data2plot1, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Spring Ice'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Spring Ice/Snow")

print(p1)


#SNOW AND SUMMER ICE

data2plot2 = pro.d2[pro.d2$Test == "Snow Vs Summer Ice",]
data2plot2 = data2plot2[data2plot2$Family %in% pro.snow.sum.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot2$Family,]
rows2add2 = data.frame(Test = "Snow Vs Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot2 = rbind(data2plot2, rows2add2)

data2plot2 =arrange(data2plot2 ,desc(Phylum),Family)
data2plot2$Family = paste0(data2plot2$Family, " (", data2plot2$Phylum, ")")

data2plot2$Family = factor(data2plot2$Family, levels = unique(data2plot2$Family))

#plot results
p2 = ggplot(data2plot2, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Summer Ice'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Summer Ice/Snow")

print(p2)

#SNOW AND CRYOCONITE

data2plot = pro.d2[pro.d2$Test == "Snow Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% pro.snow.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Snow Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p3 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Cryoconite/Snow")

print(p3)

#SPRING ICE AND SUMMER ICE

#   data2plot = pro.d2[pro.d2$Test == "Snow Vs Spring Ice",]
# data2plot = data2plot[data2plot$Family %in% pro.snow.sp.ice.tax,]

#add empty rows for all the taxa we will include
data2plot = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
data2plot = data.frame(Test = "Spring Ice Vs Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=data2plot$Domain, Phylum=data2plot$Phylum, Class=data2plot$Class, Order=data2plot$Order, Family=data2plot$Family,
                       Genus=NA, Species=NA)

# data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

#this is so we see both in the legend
data2plot = data.frame(rbind(data2plot, data2plot[1:2,]))
data2plot$log2FoldChange[nrow(data2plot)-1: nrow(data2plot)] = -0.00000000000000000000000000001
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p4 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Summer Ice'), paste0('Enriched in Spring Ice')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Summer Ice/Spring Ice")

print(p4)

#SPRING ICE AND CRYOCONITE
data2plot = pro.d2[pro.d2$Test == "Spring Ice Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% pro.sp.ice.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Spring Ice Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = -0.00000000000000000000000000001

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p5 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Spring Ice')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Cryoconite/Spring Ice")


print(p5)


#SUMMER ICE AND CRYOCONITE


data2plot = pro.d2[pro.d2$Test == "Summer Ice Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% pro.sum.ice.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Summer Ice Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = -0.00000000000000000000000000001

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p6 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Summer Ice')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-6,10)+
  ggtitle("Cryoconite/Summer Ice")

print(p6)

(p1 | p2 | p3 | p4 | p5 | p6)


#Eukaryotes

#Get results of each analysis and identify which Families overlap

#Wilcoxon
euk.wx1 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID), lev = "Family")
euk.wx1 = cbind(Test = "Snow Vs Spring Ice", euk.wx1)

euk.wx2 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID), lev = "Family")
euk.wx2 = cbind(Test = "Snow Vs Summer Ice", euk.wx2)

euk.wx3 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.wx3 = cbind(Test = "Snow Vs Cryoconite", euk.wx3)

euk.wx4 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID), lev = "Family") #no differences
euk.wx4 = cbind(Test = "Spring Ice Vs Summer Ice", euk.wx4)

euk.wx5 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.wx5 = cbind(Test = "Spring Ice Vs Cryoconite", euk.wx5)

euk.wx6 = wilcoxon.daa(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.wx6 = cbind(Test = "Summer Ice Vs Cryoconite", euk.wx6)

#DESeq2
euk.d2 = DESeq2(ps.f.euk, "Family")

#Aldex2
euk.ax1 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID), lev = "Family")
#euk.ax1 = cbind(Test = "Snow Vs Spring Ice", euk.ax1) #no significant difference

euk.ax2 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID), lev = "Family")
euk.ax2 = cbind(Test = "Snow Vs Summer Ice", euk.ax2)

euk.ax3 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.ax3 = cbind(Test = "Snow Vs Cryoconite", euk.ax3)

euk.ax4 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID), lev = "Family")
euk.ax4 = cbind(Test = "Spring Ice Vs Summer Ice", euk.ax4)

euk.ax5 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.ax5 = cbind(Test = "Spring Ice Vs Cryoconite", euk.ax5)

euk.ax6 = ALDEx2(phylo = ps.f.euk, group = c(data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID,data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID), lev = "Family")
euk.ax6 = cbind(Test = "Summer Ice Vs Cryoconite", euk.ax6)


#Which Families were found in all three analyses?

#eukaryotes

#Snow and Spring Ice
euk.snow.sp.ice.tax = Reduce(intersect, list(euk.wx1$Family,
                       euk.d2[euk.d2$Test == "Snow Vs Spring Ice",]$Family,
                       euk.ax1$Family))

#Snow and Summer Ice
euk.snow.sum.ice.tax = Reduce(intersect, list(euk.wx2$Family,
                                              euk.d2[euk.d2$Test == "Snow Vs Summer Ice",]$Family,
                                              euk.ax2$Family))

#Snow and Cryoconite
euk.snow.cry.tax = Reduce(intersect, list(euk.wx3$Family,
                                          euk.d2[euk.d2$Test == "Snow Vs Cryoconite",]$Family,
                                          euk.ax3$Family))

#Spring Ice and Summer Ice
euk.sp.sum.ice.tax = Reduce(intersect, list(euk.wx4$Family,
                                            euk.d2[euk.d2$Test == "Spring Ice Vs Summer Ice",]$Family,
                                            euk.ax4$Family))

#Spring Ice and Cryoconite
euk.sp.ice.cry.tax = Reduce(intersect, list(euk.wx5$Family,
                                            euk.d2[euk.d2$Test == "Spring Ice Vs Cryoconite",]$Family,
                                            euk.ax5$Family))

#Summer Ice and Cryoconite
euk.sum.ice.cry.tax = Reduce(intersect, list(euk.wx6$Family,
                                             euk.d2[euk.d2$Test == "Summer Ice Vs Cryoconite",]$Family,
                                             euk.ax6$Family))

#All families with significant differences
tot.fam = unique(c(euk.snow.sum.ice.tax, euk.snow.cry.tax, euk.sp.sum.ice.tax, euk.sp.ice.cry.tax, euk.sum.ice.cry.tax))

#get full tax data
tax.data = data.frame(tax_table(ps.euk))
tax.data = tax.data[tax.data$Family %in% tot.fam,]
tax.data.fam = unique(data.frame(Domain=tax.data$Domain, Phylum = tax.data$Phylum, Class=tax.data$Class, 
                                 Order=tax.data$Order,
                                 Family = tax.data$Family))
#sort by phylum
tax.data.fam = tax.data.fam[order(tax.data.fam$Phylum), ]

#Plot DESeq2 results only for those families with significant differences in all three analyses.

#SNOW AND SPRING ICE - no significant differences so making a blank plot

#add empty rows for all the taxa we will include
data2plot = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
data2plot = data.frame(Test = "Snow Vs Spring Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=data2plot$Domain, Phylum=data2plot$Phylum, Class=data2plot$Class, Order=data2plot$Order, Family=data2plot$Family,
                       Genus=NA, Species=NA)

# data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

#this is so we see both in the legend
data2plot = data.frame(rbind(data2plot, data2plot[1:2,]))
data2plot$log2FoldChange[nrow(data2plot)-1: nrow(data2plot)] = -0.00000000000000000000000000001
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p7 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Spring Ice'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)

print(p7)

#SNOW AND SUMMER ICE

data2plot = euk.d2[euk.d2$Test == "Snow Vs Summer Ice",]
data2plot = data2plot[data2plot$Family %in% euk.snow.sum.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Snow Vs Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = 0.00000000000000000000000000001

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p8 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Summer Ice'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust=0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)

print(p8)

#SNOW AND CRYOCONITE

data2plot = euk.d2[euk.d2$Test == "Snow Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% euk.snow.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Snow Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = -0.00000000000000000000000000001

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p9 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Snow')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title=element_text(hjust=0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)

print(p9)

#SPRING ICE AND SUMMER ICE

data2plot = euk.d2[euk.d2$Test == "Spring Ice Vs Summer Ice",]
data2plot = data2plot[data2plot$Family %in% euk.sp.sum.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Spring Ice Vs Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = 0.00000000000000000000000000001


data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p10 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Summer Ice'), paste0('Enriched in Spring Ice')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title=element_text(hjust=0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)

print(p10)


#SPRING ICE AND CRYOCONITE
data2plot = euk.d2[euk.d2$Test == "Spring Ice Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% euk.sp.ice.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Spring Ice Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = 0.00000000000000000000000000001

data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p11 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Spring Ice')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title=element_text(hjust=0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)

print(p11)


#SUMMER ICE AND CRYOCONITE


data2plot = euk.d2[euk.d2$Test == "Summer Ice Vs Cryoconite",]
data2plot = data2plot[data2plot$Family %in% euk.sum.ice.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Summer Ice Vs Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
# data2plot = data2plot[order(data2plot$Phylum, decreasing = TRUE), ]
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")

data2plot = data.frame(rbind(data2plot, data2plot[1,]))
data2plot$log2FoldChange[nrow(data2plot)] = -0.00000000000000000000000000001


data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

#plot results
p12 = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched in Cryoconite'), paste0('Enriched in Summer Ice')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title=element_text(hjust=0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=7),
        #legend.box="vertical", legend.margin=margin()
  )+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-6,10)


print(p12)

all.p = (p1 | p2 | p3 | p4 | p5 | p6) / (p7 | p8 | p9 | p10 | p11 | p12) +
  plot_layout(heights = c(1,0.4))

all.p

pdf("../results/daa-pro-euk.pdf", width=15, height=5)
print(all.p)
dev.off()

###### LOCATION DAA #######################################################

#Differential Abundance Analysis with Wilcoxon rank sum, DESeq2 and Aldex2

# Wilcoxon Rank Sum Test

# phylo = ps.f.pro
# group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
# lev = "Family"
# hab1 = "Snow"
# hab2 = "Spring Ice"

wilcoxon.daa <- function(phylo, group, lev){
  
  #Separate habitat data
  ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(ps, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  
  #clr transformation
  ps.clr <- microbiome::transform(ps.trim, "clr")
  
  #extract the data
  aa_data <- data.frame(t(otu_table(ps.clr)), check.names = FALSE)
  #counts <- data.frame(t(otu_table(ps.trim)))
  
  #add info to the data
  aa_data <- cbind(aa_data, 
                   Location = data.frame(sample_data(ps.clr))$Location)
  
  agg <- names(aa_data[, !names(aa_data) %in% "Location"])
  
  wilcoxon_p <- c() # Initialize empty vector for p-values
  foldchange <- c() # Initiate empty vector for effect size (difference between log means) - see https://microbiome.github.io/tutorials/all.html
  
  for (i in agg) {
    
    result <- wilcox.test(aa_data[, i] ~ Location,
                          data = aa_data)
    
    e <- coef(lm(aa_data[,i] ~ Location, data = aa_data))
    
    # Stores p-value to the vector with this column name
    wilcoxon_p[[i]]  <- result$p.value
    foldchange[[i]] <- e[[2]]
    
  }
  
  wilcoxon_p <- data.frame(taxa =  names(wilcoxon_p),
                           p_raw = unlist(wilcoxon_p))
  
  foldchange <- data.frame(taxa =  names(foldchange),
                           effect = unlist(foldchange))
  test.res = full_join(wilcoxon_p, foldchange, by = "taxa")
  
  #adjust p-values
  test.res$p_adjusted <- p.adjust(test.res$p_raw, method = "fdr")
  
  #only keep rows with significant p-values
  test.res_sig <- test.res[which(test.res$p_adjusted < 0.05 & abs(test.res$effect) > 1),]
  
  #add taxonomic info
  tax.info = data.frame(tax_table(ps.clr))
  tax.info$taxa = rownames(tax.info)
  test.res_sig_tax = left_join(test.res_sig, tax.info, by="taxa")
  
  #remove those that say uncultured or of unknown placement
  test.res_sig_tax_rm = test.res_sig_tax[test.res_sig_tax$Family != "uncultured",]
  test.res_sig_tax_rm = test.res_sig_tax[test.res_sig_tax$Family != "Incertae_Sedis",]
  
  return(test.res_sig_tax_rm)
  
}


# DESeq2
# http://127.0.0.1:19618/library/DESeq2/doc/DESeq2.html
# phylo = ps.f.pro
# group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
# lev = "Family"
# hab1 = "Snow"
# hab2 = "Spring Ice"

DESeq2 <- function(phylo, group, lev){
  
  #Separate habitat data
  ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(ps, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  
  #convert phylo object to deseq2 object
  #no need to normalise the data
  dds = phyloseq_to_deseq2(ps.trim, ~ Location)
  dds #look at features of our deseq2 object 
  
  #side note - look at the distribution of our count data
  #we can see that it is NOT normally distributed 
  #therefore we cannot use methods for normally distributed data to predict the likelihood
  #of events or features like the mean/sd etc
  #poisson distributions are more appropriate to fit for microbiome/RNA-Seq data
  
  # hist(rowSums(data.frame(otu_table(ps.trim))))
  
  #because our count data is sparce DESeq has issues calculating the geometric mean 
  #so we need to calculate the size factors using a zero-tolerant method
  #calculate geometric means prior to estimate size factors 
  #(https://github.com/joey711/phyloseq/issues/387)
  
  # gm_mean = function(x, na.rm=TRUE){
  # exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  # }
  # 
  # geoMeans = apply(counts(dds), 1, gm_mean)
  
  # dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  
  #carry out test
  #the DESeq2 differential expression analysis steps are wrapped into a single function DESEq.
  #we have performed this on all four conditions (Habitats) and from this we can extract results for all comparisosns
  dds = DESeq(dds) #comparison is based on alphabetic order of the condition
  dds
  res = results(dds, cooksCutoff = FALSE)
  alpha = 0.05
  effect.size = 1
  sigtab = res[which(res$padj < alpha & abs(res$log2FoldChange) > effect.size), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.trim)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  
  
  #put all our results into one table
  data = sigtab
  
  #remove NA data
  data = data[!is.na(data$Family) ,]
  
  #remove uncultured data
  data = data[data$Family != "uncultured",]
  data = data[data$Family != "Incertae_Sedis",]
  
  return(data)
  
}

# Aldex2
# Estimates within sample variation using the Dirichlet distribution.
# Applied centrered-log-ratio transformation.
# If comparing two conditions Aldex2 uses Welch's t-test and Wilcoxon tests.
# If comparing three or more conditions it will perform a one-way ANOVA and Kruskall-Wallis test.
# The Benjamini-Hochberg p-value adjustment is used by default.
# https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html

# phylo = ps.f.pro
# group = c(data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
# lev = "Family"
# hab1 = "Snow"
# hab2 = "Spring Ice"

 ALDEx2 <- function(phylo, group, lev){
  
  set.seed(666)
   
  #Separate habitat data
  ps <- prune_samples(group, phylo)
  
  #agglomerate
  ps.glom <- tax_glom(ps, lev)
  
  #remove groups with zero counts
  ps.trim <- filter_taxa(ps.glom, function(x) sum(x) > 0, TRUE)
  
  #extract count table of un-normalised data with samples as columns
  counts <- data.frame(otu_table(ps.trim), check.names = FALSE)
  
  #get vector of conditions in the same order as the samples in the counts table
  conds <- c(data.frame(sample_data(ps.trim))$Location)

  #aldex is a wrapper function linking the modular elements of the aldex2 package. 
  #aldex function performs a two-sample t-test and calculates effect sizes. 
  #if test="t" then effect should = TRUE
  #"t" option evaluates the data as a two-factor experiment using Welch's t-test and Wilcoxon rank tests 
#Welch's t test = null hypothesis two population have equal means (the populations are assumed to be normal distributed but don't need to have equal variances)
#Wilcoxon rank tests = null hypothesis two populations have equal means (there are no assumptions about the distribution of the data)
#For multiple sample tests using ANOVA (to test difference in means of multiple groups) effect should = FALSE was effect size should not be calculated.
#Effect size = value measuring strength of the relationship between two variables in a population. This is different from p-value which tells you the effect exists. p-value = statistical significance (there is a difference). effect size = practical significance (the difference is big enough to care about).  
#test = "kw" evaluate the data as a one-way ANOVA (Analysis of Variance - how does ONE independent factor, e.g. habitat, impact the dependent variable, e.g., CLR values. A two-way ANOVA evaluated the impact of TWO independent variables on a dependent variable) using the glm and Kruskal-Wallace tests (tests for significant differences between the medians of three or more independent groups consisting of non-parametric data).
#All tests include a Benjamini-Heckberg correction of raw p-values. 
#data can be plotted into Bland-Altman (MA) or effect (MW) plots for two0way tests. 
#Bland-Altmann plots are difference plots. There are used to evaluate the agreement between two measurements.
#as we want to compare pairs of habitats here we will use the following arguments:
#test="t" - compare two conditions
#effect=TRUE - calculate the effect size of these differences
#mc.samples = 1000 - the number of Dirichlet Monte-Carlo Instances/Samples to use when estimating the underlying distributions. These are Monte-Carlo method random sampling of the Dirichlet probability distribution (to estimate the chance of getting our result at random).
#include.sample.summary = FALSE - whether to include median clr values for each sample
#verbose = TRUE - print diagnostic information while running functions
#denom = "all" - which features to retain as the denominator for the Geometric Mean calculation
#paired.test = FALSE - are the samples paired
x.all <- aldex(counts, conds, mc.samples=1000, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=TRUE, paired.test=FALSE)

#insignficant features (families) are shown in grey or black. Statistically significant features (those families which are significant different between habitats) are shown in red
#log-ratio abundance axis is the clr value
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="wilcox", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="wilcox", xlab="Dispersion",
           ylab="Difference")

#get only our significant values and those with effect size > 1
#we want to look at wi.eBH = expected BH corrected p value of Wilcoxon test 
res.trim = rownames_to_column(x.all, "ID") %>%
  filter(wi.eBH <= 0.05, abs(effect) > 1) # here we chose the wilcoxon output rather than ttest #keep the family ID, wilcoxon BH adj p-value, dif.btw (difference between the group median for each significant feature), effect (effect size), overlap (proportion of effect size that overlaps 0))

#add taxonomic info
tax.info = data.frame(tax_table(ps.trim))
tax.info$ID = rownames(tax.info)
res.trim_tax = left_join(res.trim, tax.info, by="ID")

return(res.trim_tax)

}


#Prokaryotes

#Get results of each analysis and identify which Families overlap

#Wilcoxon
pro.wx1 = wilcoxon.daa(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID, lev = "Family") #no differences
#pro.wx1 = cbind(Test = "Snow", pro.wx1) #no significant differences

pro.wx2 = wilcoxon.daa(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID, lev = "Family")
pro.wx2 = cbind(Test = "Spring Ice", pro.wx2)

pro.wx3 = wilcoxon.daa(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID, lev = "Family")
pro.wx3 = cbind(Test = "Summer Ice", pro.wx3)

pro.wx4 = wilcoxon.daa(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID, lev = "Family") 
pro.wx4 = cbind(Test = "Cryoconite", pro.wx4)



#DESeq2
pro.dq1 = DESeq2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID, lev = "Family")
pro.dq1 = cbind(Test = "Snow", pro.dq1)

pro.dq2 = DESeq2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID, lev = "Family")
pro.dq2 = cbind(Test = "Spring Ice", pro.dq2)

pro.dq3 = DESeq2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID, lev = "Family")
pro.dq3 = cbind(Test = "Summer Ice", pro.dq3)

pro.dq4 = DESeq2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID, lev = "Family")
pro.dq4 = cbind(Test = "Cryoconite", pro.dq4)


#Aldex2
pro.ax1 = ALDEx2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID, lev = "Family")
#pro.ax1 = cbind(Test = "Snow", pro.ax1) #no significant difference

pro.ax2 = ALDEx2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID, lev = "Family")
pro.ax2 = cbind(Test = "Spring Ice", pro.ax2)

pro.ax3 = ALDEx2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID, lev = "Family")
pro.ax3 = cbind(Test = "Summer Ice", pro.ax3)

pro.ax4 = ALDEx2(phylo = ps.f.pro, group = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID, lev = "Family")
pro.ax4 = cbind(Test = "Cryoconite", pro.ax4)

#Which Families were found in all three analyses?

#Prokaryotes

#Snow
pro.snow.tax = Reduce(intersect, list(pro.wx1$Family,
                                      pro.dq1$Family,
                                      pro.ax1$Family))

#Spring Ice
pro.sp.ice.tax = Reduce(intersect, list(pro.wx2$Family,
                                        pro.dq2$Family,
                                        pro.ax2$Family))

#Summer Ice
pro.sum.ice.tax = Reduce(intersect, list(pro.wx3$Family,
                                         pro.dq3$Family,
                                         pro.ax3$Family))

#Cryoconite
pro.cry.tax = Reduce(intersect, list(pro.wx4$Family,
                                     pro.dq4$Family,
                                     pro.ax4$Family))

#All families with significant differences
tot.fam = unique(c(pro.sp.ice.tax, pro.sum.ice.tax, pro.cry.tax))

#get full tax data
tax.data = data.frame(tax_table(ps.pro))
tax.data = tax.data[tax.data$Family %in% tot.fam,]
tax.data.fam = unique(data.frame(Domain=tax.data$Domain, Phylum = tax.data$Phylum, Class=tax.data$Class, 
                                 Order=tax.data$Order,
                                 Family = tax.data$Family))

#Plot DESeq2 results only for those families with significant differences in all three analyses.

#SPRING ICE

data2plot = pro.dq2[pro.dq2$Test == "Spring Ice",]
data2plot = data2plot[data2plot$Family %in% pro.sp.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Spring Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p1b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-10,12)+
  ggtitle("Spring Ice")

print(p1b)


#SUMMER ICE

data2plot = pro.dq3[pro.dq3$Test == "Summer Ice",]
data2plot = data2plot[data2plot$Family %in% pro.sum.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p2b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-10,12)+
  ggtitle("Summer Ice")

print(p2b)

#CRYOCONITE
data2plot = pro.dq4[pro.dq4$Test == "Cryoconite",]
data2plot = data2plot[data2plot$Family %in% pro.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p3b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=3))+
  ylab("log2FoldChange")+
  ylim(-10,12)+
  ggtitle("Cryoconite")

print(p3b)


p1b | p2b | p3b

#Eukaryotes

#Get results of each analysis and identify which Families overlap

#Wilcoxon
euk.wx1 = wilcoxon.daa(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID, lev = "Family") #no differences
#euk.wx1 = cbind(Test = "Snow", euk.wx1) #no significant difference

euk.wx2 = wilcoxon.daa(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID, lev = "Family")
euk.wx2 = cbind(Test = "Spring Ice", euk.wx2)

euk.wx3 = wilcoxon.daa(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID, lev = "Family")
euk.wx3 = cbind(Test = "Summer Ice", euk.wx3)

euk.wx4 = wilcoxon.daa(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID, lev = "Family") 
euk.wx4 = cbind(Test = "Cryoconite", euk.wx4)



#DESeq2
euk.dq1 = DESeq2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID, lev = "Family")
euk.dq1 = cbind(Test = "Snow", euk.dq1)

euk.dq2 = DESeq2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID, lev = "Family")
euk.dq2 = cbind(Test = "Spring Ice", euk.dq2)

euk.dq3 = DESeq2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID, lev = "Family")
euk.dq3 = cbind(Test = "Summer Ice", euk.dq3)

euk.dq4 = DESeq2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID, lev = "Family")
euk.dq4 = cbind(Test = "Cryoconite", euk.dq4)


#Aldex2
euk.ax1 = ALDEx2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID, lev = "Family")
#euk.ax1 = cbind(Test = "Snow", euk.ax1) #no significant difference

euk.ax2 = ALDEx2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID, lev = "Family")
euk.ax2 = cbind(Test = "Spring Ice", euk.ax2)

euk.ax3 = ALDEx2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID, lev = "Family")
euk.ax3 = cbind(Test = "Summer Ice", euk.ax3)

euk.ax4 = ALDEx2(phylo = ps.f.euk, group = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID, lev = "Family")
euk.ax4 = cbind(Test = "Cryoconite", euk.ax4)

#Which Families were found in all three analyses?

#eukaryotes

#Snow
euk.snow.tax = Reduce(intersect, list(euk.wx1$Family,
                                      euk.dq1$Family,
                                      euk.ax1$Family))

#Spring Ice
euk.sp.ice.tax = Reduce(intersect, list(euk.wx2$Family,
                                        euk.dq2$Family,
                                        euk.ax2$Family))

#Summer Ice
euk.sum.ice.tax = Reduce(intersect, list(euk.wx3$Family,
                                         euk.dq3$Family,
                                         euk.ax3$Family))

#Cryoconite
euk.cry.tax = Reduce(intersect, list(euk.wx4$Family,
                                     euk.dq4$Family,
                                     euk.ax4$Family))

#All families with significant differences
tot.fam = unique(c(euk.sp.ice.tax, euk.sum.ice.tax, euk.cry.tax))

#get full tax data
tax.data = data.frame(tax_table(ps.euk))
tax.data = tax.data[tax.data$Family %in% tot.fam,]
tax.data.fam = unique(data.frame(Domain=tax.data$Domain, Phylum = tax.data$Phylum, Class=tax.data$Class, 
                                 Order=tax.data$Order,
                                 Family = tax.data$Family))


#Plot DESeq2 results only for those families with significant differences in all three analyses.

#SPRING ICE

data2plot = euk.dq2[euk.dq2$Test == "Spring Ice",]
data2plot = data2plot[data2plot$Family %in% euk.sp.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Spring Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p4b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=5))+
  ylab("log2FoldChange")+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylim(-10,12)

print(p4b)


#SUMMER ICE

data2plot = euk.dq3[euk.dq3$Test == "Summer Ice",]
data2plot = data2plot[data2plot$Family %in% euk.sum.ice.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Summer Ice", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p5b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=11))+
  ylab("log2FoldChange")+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylim(-10,12)

print(p5b)

#CRYOCONITE
data2plot = euk.dq4[euk.dq4$Test == "Cryoconite",]
data2plot = data2plot[data2plot$Family %in% euk.cry.tax,]

#add empty rows for all the taxa we will include
rows2add = tax.data.fam[!tax.data.fam$Family %in% data2plot$Family,]
rows2add2 = data.frame(Test = "Cryoconite", baseMean=0, log2FoldChange=0, lfcSE=0, stat=0, pvalue=0, padj=0, Domain=rows2add$Domain, Phylum=rows2add$Phylum, Class=rows2add$Class, Order=rows2add$Order, Family=rows2add$Family,
                       Genus=NA, Species=NA)
data2plot = rbind(data2plot, rows2add2)
data2plot =arrange(data2plot ,desc(Phylum),Family)
data2plot$Family = paste0(data2plot$Family, " (", data2plot$Phylum, ")")
data2plot$Family = factor(data2plot$Family, levels = unique(data2plot$Family))

p6b = ggplot(data2plot, aes(x=Family, y=log2FoldChange)) + 
  geom_bar(aes(fill = log2FoldChange < 0), stat = "identity") + 
  scale_fill_manual(drop = FALSE, breaks = c(TRUE, FALSE), values=c("gray", "red"), labels=c(paste0('Enriched on Lower Foxfonna'), paste0('Enriched on Upper Foxfonna')))+
  coord_flip()+
  theme(legend.position = "none", legend.title = element_blank(), plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=5))+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  ylab("log2FoldChange")+
  ylim(-10,12)

print(p6b)



all.p2 = (p1b | p2b | p3b) / (p4b | p5b | p6b) +
  plot_layout(heights = c(1,0.5))

pdf("../results/daa-pro-euk-location.pdf", width=10, height=10)
print(all.p2)
dev.off()


#combine all plots with titles

#my plots
p1
p2
p3
p4
p5
p6
p1b
p2b
p3b
p4b
p5b
p6b
all.p
all.p2

#get our titles
title1 <- ggdraw() + draw_label("DAA Between Habitats", size=20, hjust=2.8, fontface = "bold")
title2 <- ggdraw() + draw_label("DAA Between Upper/Lower Foxfonna", size=20, hjust=1.7, fontface = "bold")

combi.plot = title1 / plot_spacer() / ((p1 + labs(title="Spring Ice Vs \n Snow") + theme(axis.text = element_text(size=10), axis.title.y = element_text(size=20), plot.title = element_text(size=20))) | (p2 + labs(title="Summer Ice Vs \n Snow") + theme(plot.title = element_text(size=20))) | (p3 + labs(title="Cryoconite Vs \n Snow") + theme(plot.title = element_text(size=20))) | (p4 + labs(title="Summer Ice Vs \n Spring Ice") + theme(plot.title = element_text(size=20))) | (p5 + labs(title="Cryoconite Vs \n Spring Ice") + theme(plot.title = element_text(size=20))) | (p6 + labs(title="Cryoconite Vs \n Summer Ice") + theme(plot.title = element_text(size=20)))) / ((p7 + theme(legend.text = element_text(size=10), axis.text = element_text(size=10), axis.title.y = element_text(size=20)))| (p8 + theme(legend.text = element_text(size=10))) | (p9 + theme(legend.text = element_text(size=10))) | (p10 + theme(legend.text = element_text(size=10))) | (p11 + theme(legend.text = element_text(size=10))) | (p12 + theme(legend.text = element_text(size=10)))) / plot_spacer() /
  title2 / plot_spacer() /((p1b + theme(axis.text = element_text(size=10), axis.title.y = element_text(size=20), plot.title = element_text(size=25))) | (p2b + theme(plot.title = element_text(size=25))) | (p3b + theme(plot.title = element_text(size=25)))) / ((p4b + theme(axis.text = element_text(size=10), axis.title.y = element_text(size=20))) | (p5b + theme(legend.text = element_text(size=20))) | p6b) +
  plot_layout(heights = c(0.1, 0.1, 1, 0.4, 0.1, 0.1, 0.1, 2.7, 1.4))

pdf("../results/daa-plot.pdf", width=15, height=25)
print(combi.plot)
dev.off()
