#Script to analyse PMA-qPCR test data

#Analysis steps
#1. Clear workspace, load packages and data
#2. Total genes in reaction: Multiply mean gene copies per uL by the number of uL added to the reaction (=4uL)
#3. Genes per ng DNA: Divide total genes in the reaction by the total ng DNA added to the reaction
#4. Gene copies per uL: Multiply gene copies per ng DNA by original sample DNA concentration (remember to multiply by the X-fold dilution if necessary)
#5. Gene copies in total elute: Multiply gene copies per uL by the total volume of DNA extraction elute = 100uL
#6. Gene copies per mg: Divide gene copies in total elute by the original sample weight
#7. Calculate basic stats
#8. Plot the data
#9. Test for significant differences 
#Great tutorials
#http://www.sthda.com/english/wiki/one-way-anova-test-in-r
#http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r

#1. Clear workspace, load packages and data
rm(list=ls())
graphics.off()

library(tidyverse)
library(dplyr)
library(dunn.test)
library(ggpubr)

data = read.csv("../data/PMA-qPCR-test-results.csv")

#2. Total genes in reaction: Multiply mean gene copies per uL by the number of uL added to the reaction (=4uL)

#dataframe into tibble
my_data <- as_tibble(data)
#remove controls
my_data = my_data[my_data$Group != "Control",]
#only keep distinct rows with mean cell abundances
my_data2 = my_data[!duplicated(my_data$Quantity.Mean),]
#get total gene in reaction 
my_data2$GenesInReaction = my_data2$Quantity.Mean*4

#3. Genes per ng DNA: Divide total genes in the reaction by the total ng DNA added to the reaction
my_data2$GenesPerNg = my_data2$GenesInReaction / (as.numeric(my_data2$FinalDNAConcentration)*4)

#4. Gene copies per uL: Multiply gene copies per ng DNA by original sample DNA concentration (remember to multiply by the X-fold dilution in necessary)
my_data2$GenesPeruL = my_data2$GenesPerNg * (as.numeric(my_data2$FinalDNAConcentration))
#multiply diluted samples by 20
my_data2 = my_data2 %>%
  mutate(GenesPeruLCorrected = case_when(Diluted == TRUE ~ GenesPeruL*20,
                            Diluted == FALSE ~ GenesPeruL*1))    

#5. Gene copies in total elute: Multiply gene copies per uL by the total volume of DNA extraction elute = 100uL
my_data2$GenesPerElute = my_data2$GenesPeruLCorrected * 100

#6. Gene copies per mg: Divide gene copies in total elute by the original sample weight
my_data2$GenesPerMg = my_data2$GenesPerElute / my_data2$E.coli.Weight..mg.

#7. Calculate basic stats

#when calculating my stats I see that the standard deviations (measure of variation around the mean) is very high for my 
#dataset. So I calculated the coefficient of variation by dividing SD by mean. CV ≥ 1 indicates a relatively high variation.
#CV = the ratio between SD and Mean in the dataset. Essentially we see here that our data is spread over a wide range of values,
#but the variance isn't "high" as most groups have CV < 1.

#summary stats
my_stats = group_by(my_data2, Group, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(GenesPerMg, na.rm = TRUE),
    median = median(GenesPerMg, na.rm = FALSE),
    sd = sd(GenesPerMg, na.rm = TRUE)
  )

#calculate the coefficient of variation
my_stats$CV = my_stats$sd / my_stats$mean

#export data
write.csv(my_stats, "../results/PMA-qPCR-test-res.csv")

#8. Plot the data

#create a combined group and treatment categorical variable
my_data2$CombiGroups = paste(my_data2$Group, my_data2$Treatment)
#put them in the correct order for plotting
my_data2$CombiGroups <- factor(my_data2$CombiGroups , levels=c("Live tDNA", "Live iDNA", "Killed tDNA", "Killed iDNA"))
#add log10 transformed genes
my_data2$log10GenesPerMg = log10(my_data2$GenesPerMg)

p <- ggplot(my_data2, aes(x=CombiGroups, y=log10GenesPerMg)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  ylab("log10 Genes per Mg")
p

#save the plot
pdf("../results/PMA-qPCR-results.pdf", height=3, width=4)
print(p)
dev.off()

#we can see from the high CV value and from the plot that we have an anomolous result with sample KP2. 
#Let's remove this and try again

#Summary stats after removing anomolous result
my_data3 = my_data2[my_data2$Sample != "KP2",]

#summary stats
my_stats2 = group_by(my_data3, Group, Treatment) %>%
  summarise(
    count = n(),
    mean = mean(GenesPerMg, na.rm = TRUE),
    median = median(GenesPerMg, na.rm = FALSE),
    sd = sd(GenesPerMg, na.rm = TRUE)
  )

#calculate the coefficient of variation
my_stats2$CV = my_stats2$sd / my_stats2$mean

#export data
write.csv(my_stats2, "../results/PMA-qPCR-test-res-no-anom.csv")

#8. Plot the data

p2 <- ggplot(my_data3, aes(x=CombiGroups, y=log10GenesPerMg, fill=Treatment)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.position="bottom")+
  ylab("log10 16S rRNA Genes per Mg E.coli Cells")
p2

#save the plot
pdf("../results/PMA-qPCR-results-no-anom.pdf", height=3, width=4)
print(p2)
dev.off()

#9. Test for significant differences
#are there differences in the genes per mg between Live iDNA, Live tDNA, Killed iDNA and Killed tDNA?
#Which groups are significantly different?
# Null hypothesis: the means of the different groups are the same
# Alternative hypothesis: At least one sample mean is not equal to the others

#PARAMETRIC TESTS
#Compute the analysis of variance to see if there are significant differences between the groups
res.aov1 <- aov(GenesPerMg ~ CombiGroups, data = my_data3) 

#but is our data normally distibuted?
# 1. Homogeneity of variances
plot(res.aov1, 1)
# 2. Normality
plot(res.aov1, 2)
# Extract the residuals
aov_residuals1 <- residuals(object = res.aov1)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals1) #p < 0.05 so data not normally distributed

#DATA NOT NORMALLY DISTRIBUTED SO USE NON-PARAMETRIC TEST 

#NON PARAMETRIC TESTS
#this performs kruskal-wallis test for differences between >2 groups as well as post hoc dunn test for pairwise differences
dunn.test(my_data3$GenesPerMg, my_data3$CombiGroups, method = "none", altp=TRUE)

#we have significant differences between live iDNA and Live tDNA 
#and near signficiant differences between killed iDNA and killed tDNA 
#but this is without p value corrections for multiple comparisons