#This script provides:
#1. Tests for qPCR inhibition due to PMA inteference.

#plot 1: undiluted and normalised diluted samples plot as a scatter plot with a line of equality to show if, after normalising, 
#the diluted samples had similar qPCR values (e.g. the dots were on the line of equality).
#plot 2. Bland-Altman difference plot (a.k.a Tukey mean-difference plot), is a method of plotting data to see how much they agree!
#it plots the mean of two datapoints, by their difference.


#Clear workspace
rm(list=ls())
graphics.off()

#load packages
library(ggplot2)
library(PairedData)
library(dplyr)
library(reshape)
library(scales)

#load data

#results of inhibition test
inhib_test <- read.csv("../data/Solman-Standards-and-Inhibition-Test-090823_results.csv")

#qpcr metadata
metadata = read.csv("../data/qpcr-metadata.csv", sep=",")

##################################################################################################################

#1. Tests for qPCR inhibition due to PMA interference.

#INHIBITION TEST RESULTS

#I included 10-fold dilutions for my samples
#I need to calculate the mean for each true sample and the mean for each dilution and find the percentage difference

inhib.data = inhib_test

#get column names
names(inhib.data) = inhib.data[7,]

#remove unwanted rows
inhib.data = inhib.data[8:(nrow(inhib.data)),]

#remove rows with empty "Task" Column
inhib.data = inhib.data[inhib.data$Task != "",]

#okay keep "UNKNOWN SAMPLES"
inhib.data.unknown = inhib.data[inhib.data$Task == "UNKNOWN",]

#remove rows with undetermined in with Ct
inhib.data.unknown = inhib.data.unknown[inhib.data.unknown$Cт != "Undetermined",]

#get only the data we need for this analysis
data = as.data.frame(unique(cbind(inhib.data.unknown$`Sample Name`, inhib.data.unknown$`Quantity Mean`)))
data$group = c("Undiluted", "Undiluted", "Diluted", "Undiluted", "Undiluted", "Diluted", "Diluted",
               "Diluted", "Undiluted", "Undiluted", "Diluted", "Diluted", "Undiluted", "Undiluted", "Diluted", "Diluted")

#give sample names
df_new <- data %>% mutate(across(c('V1'), substr, 5, nchar(V1)))
nams = gsub("\\.","", df_new$V1)
nams = gsub("D","", nams)
data$Sample.Name = nams

#add metadata
metadata$Sample.Name = gsub(" ", "", metadata$Sample.Name) #remove space from metadata names
df = left_join(data, metadata, by="Sample.Name")

#back calculate gene abundances
df$DNA.Concentration = gsub("Below detection", 0.0025, df$DNA.Concentration) #replace below detection limit with half the detection limit
df$GeneCopiesPerReaction = as.numeric(df$V2)*4 #calculate gene copies per reaction
df$GeneCopiesPerNgDNA = df$GeneCopiesPerReaction/(as.numeric(df$DNA.Concentration)*4) #calculate gene copies per ng DNA
df$GeneCopiesPeruLElute = df$GeneCopiesPerNgDNA * as.numeric(df$DNA.Concentration) #calculate gene copies per uL elute
for (i in 1:nrow(df)){ #multiply by volume of elute to get total gene copies extracted
  if (df$Double.Extraction[i] == TRUE){
    df$TotalGeneCopiesExtracted[i] = df$GeneCopiesPeruLElute[i] * 200
  } else if (df$Double.Extraction[i] == FALSE){
    df$TotalGeneCopiesExtracted[i] = df$GeneCopiesPeruLElute[i] * 100
  } 
}
df$GeneCopiesPerMgMl = df$TotalGeneCopiesExtracted / df$Volume.Weight #find the gene copies per mg or ml

#multiplying by 10 if the sample is a diluted sample
df$NormalisedGeneCopies = df$GeneCopiesPerMgMl
for (i in 1:nrow(df)){
  if (df$group[i] == "Diluted"){
    df$NormalisedGeneCopies[i] = df$NormalisedGeneCopies[i] * 10
  }
}


#add log transformed data
df$log10genecopies = log10(df$NormalisedGeneCopies)

#Plot undiluted mean log10 gene copies mL-1 for each sample against mean log10 diluted gene copies mL-1. 
#make dataset wide for plotting
plot.df = as.data.frame(cbind(df$Sample.Name, df$group, as.numeric(df$log10genecopies)))
data_wide <- reshape(plot.df, idvar = "V1", timevar="V2", direction="wide")
names(data_wide) = c("Sample", "Undiluted", "Diluted")

#plot
p = ggplot(data_wide, aes(x=as.numeric(Undiluted), y=as.numeric(Diluted))) + 
  geom_point(size=6) +
  geom_abline(slope=1, intercept = 0)+
  theme_bw()+
  xlab(expression(paste("Undiluted log"[10]," 16S rRNA gene copies per mL/mg")))+
  ylab(expression(paste("Diluted log"[10]," 16S rRNA gene copies per mL/mg * 10")))
p

pdf("../results/inhibition-test-graph1.pdf", height=6, width=6)
print(p)
dev.off()

#graph shows that diluted and undiluted samples had very similar results.

#plot 2. Bland-Altman difference plot

#let's get our data again but not the log10 data, just the real data
plot.df2 = as.data.frame(cbind(df$Sample.Name, df$group, as.numeric(df$NormalisedGeneCopies)))
data_wide2 <- reshape(plot.df2, idvar = "V1", timevar="V2", direction="wide")
names(data_wide2) = c("Sample", "Undiluted", "Diluted")
#make sure our data is numeric
data_wide2$Undiluted = as.numeric(data_wide2$Undiluted) 
data_wide2$Diluted = as.numeric(data_wide2$Diluted) 

#get the difference between each pair of values
data_wide2$dif = data_wide2$Undiluted - data_wide2$Diluted

#get the mean for each pair of values
data_wide2$mean = (data_wide2$Undiluted + data_wide2$Diluted)/2

#plot the log10 mean (x axis) by the log 10 difference (y axis)

plot(log10(data_wide2$mean), log10(data_wide2$dif))

#alternatively we cn use the log10 values to calculate the mean and the difference
data_wide$log10Dif = data_wide$Undiluted - data_wide$Diluted

#Plot mean log10 gene copies mL-1 for each sample (using diluted and undiluted samples) against the log10 difference between diluted/undiluted samples to see if there is an increase/decrease in difference (i.e. inhibition) with log10 gene copies mL-1.

#Calculate mean log10 gene copies mL-1 for each sample using normalised diluted and undiluted samples. 
#get mean values for each sample
# data.sum2 = data.unknown %>%
#   group_by(R_Names) %>%
#   summarise(Mean = mean(NormDilution, na.rm=TRUE),
#             log10Mean = log10(mean(NormDilution, na.rm=TRUE)),
#             SD = sd(NormDilution, na.rm=TRUE))
# 
# #make dataset wide for calculating differences
# data_wide2 <- spread(data.sum[,-c(3,5,6)], Dilution, Mean)
# data_wide2$Dif = data_wide2$Diluted - data_wide2$Undiluted
# data_wide2$log10Diluted = log10(data_wide2$Diluted)
# data_wide2$log10Undiluted = log10(data_wide2$Undiluted)
# data_wide2$log10Dif = data_wide2$log10Diluted - data_wide2$log10Undiluted

#find the mean value using diluted and undiluted samples 
plot.df2 = as.data.frame(cbind(df$Sample.Name, df$group, as.numeric(df$NormalisedGeneCopies)))
data_wide2 <- reshape(plot.df2, idvar = "V1", timevar="V2", direction="wide")
names(data_wide2) = c("Sample", "Undiluted", "Diluted")
data_wide2$Mean = (as.numeric(data_wide2$Undiluted) + as.numeric(data_wide2$Diluted)) / 2
data_wide2$log10Mean = log10(data_wide2$Mean)

#find the difference between the diluted and undiluted samples
data_wide2$log10Dif = log10(as.numeric(data_wide2$Undiluted) - as.numeric(data_wide2$Diluted))
mean_dif = mean(data_wide2$log10Dif)

#confidence intervals
lower_limit <- mean_dif - 1.91*sd(data_wide2$log10Dif)
upper_limit <- mean_dif + 1.91*sd(data_wide2$log10Dif)

#plot

p2 = ggplot(data_wide2, aes(x=log10TotalMean, y=log10Dif)) + 
  geom_point(size=4) +
  geom_hline(yintercept=0)+
  geom_hline(yintercept=mean_dif, linetype="dashed", color = "red")+
  geom_hline(yintercept=lower_limit, color = "blue")+
  geom_hline(yintercept=upper_limit, color = "blue")+
  theme_bw()+
  xlab(expression(paste("Mean Log"[10]," DNA quantified by undiluted and diluted qPCR")))+
  ylab(expression(paste("Difference in Log"[10]," DNA quantity (16S rRNA gene copies per mL"^-1, "/mg"^-1)))
p2


pdf("../results/inhibition-test-graph2.pdf", height=6, width=6)
print(p2)
dev.off()

#qPCR inhibition will be considered present if there is a one log increase in gene copies mL-1 in normalised diluted samples compared to undiluted samples. For example, if mean undiluted sample quantity is 1 x 10^5 gene copies mL-1 and the mean normalised diluted sample quantity is 1x10^6 gene copies mL-1 we would confirm there is a qPCR inhibition due to copurified substances.

