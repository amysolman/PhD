#This script provides:
#1. Removes anomalous results.
#2. Calculate proportion of reads for each sample from iDNA.
#3. Plots gene abundances for iDNA and tDNA for each sample. 
#4. Calculates mean and sd gene abundances per mL or mg for treated (iDNA) and untreated (tDNA) samples for each habitat.
#5. Performs statistical test for differences in abundance of iDNA and tDNA genes for all habitats, spring habitats and summer habitats.


#Clear workspace
rm(list=ls())
graphics.off()

#load packages
library(ggplot2)
library(PairedData)
library(dplyr)
library(reshape)

#load data
df = read.csv("../data/meta_qpcr_res.csv", sep="\t")
df2 = read.csv("../data/qPCR-mean.csv", sep="\t")

#set the habitat order
df$Habitat = factor(df$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite", "Control"))
df2$Habitat = factor(df2$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite", "Control"))

##################################################################################################################
#1. Removes anomalous results.
#remove the anomolous result
df = df[df$Sample != "S21.30",]

##################################################################################################################
#2. Calculate proportion of reads for each sample from iDNA.
###################What proportion of reads from each sample was from tDNA and iDNA??????
#group by sample and treatment
#subtract iDNA reads from tDNA reads
df.p = df %>%
    group_by(Sample, Treatment) %>%
    summarise(reads = sum(Gene.Copies.Per.mL.mg))

#into wide format
#reshape(df.p, idvar = "Sample", timevar = "Treatment", direction = "wide")

dif = vector() #vector for storing the difference in number of reads from sum
pro = vector()

for (i in 1:length(unique(df.p$Sample))){
  samps = unique(df.p$Sample)
  data = df.p[df.p$Sample == samps[i],]
  print(data)
  dif[i] = data$reads[2] - data$reads[1]
  pro[i] = data$reads[1] / data$reads[2]
}

x = cbind(unique(df.p$Sample), dif, pro)

#save these results
write.csv(x, "../results/pma-proportion-idna-tdna.csv")

##################################################################################################################

#3. Plots gene abundances for iDNA and tDNA for each sample. 

#plot
p = ggplot(df) +
  geom_point(aes(Sample.Name.To.Plot, log10(Gene.Copies.Per.mL.mg), color = Treatment, shape = Treatment, stroke=2), size=4) +
  ylab("log10 Gene Copies per mL or mg")+
  xlab("Sample")+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  guides(color=guide_legend(override.aes=list(shape=c(21,21), size=4, stroke=2)))+
  scale_shape_manual(values = c(21,21))+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, size=10),
        legend.title=element_blank(),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_text(size=15))

#save the plot
pdf("../results/pma-qPCR-plot.pdf", height=5, width=15)
print(p)
dev.off()

#second plot
p2 = ggplot(df2) +
  geom_point(aes(Sample, log10(Gene.Copies.Per.mL.mg), color = Treatment, shape = Treatment, stroke=2), size=6) +
  ylab("log10 Gene Copies per mL or mg")+
  xlab("Sample")+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  guides(color=guide_legend(override.aes=list(shape=c(21,24), size=6, stroke=2)))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5, size=10),
        legend.title=element_blank(),
        legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        strip.text = element_text(size=15))

#save the plot
pdf("../results/pma-qPCR-plot2.pdf", height=5, width=15)
print(p2)
dev.off()

#another way to plot paired data

# Subset gene data without treatment
tdna <- subset(df2,  Treatment == "tDNA", Gene.Copies.Per.mL.mg,
                 drop = TRUE)
# subset gene data with treatment
idna <- subset(df2,  Treatment == "iDNA", Gene.Copies.Per.mL.mg,
                drop = TRUE)
# Plot paired data
pd <- paired(tdna, idna)
plot(pd, type = "profile") + theme_bw()

##################################################################################################################

#4. Calculates mean and sd gene abundances per mL or mg for treated (iDNA) and untreated (tDNA) samples for each habitat.

#basic stats
sum_stats = group_by(df2, Treatment, Habitat) %>%
  summarise(
    count = n(),
    mean = mean(Gene.Copies.Per.mL.mg, na.rm = TRUE),
    median = median(Gene.Copies.Per.mL.mg, na.rm = TRUE),
    IQR = IQR(Gene.Copies.Per.mL.mg, na.rm = TRUE),
    sd = sd(Gene.Copies.Per.mL.mg, na.rm = TRUE)
  )

#save results
write.csv(sum_stats, "../results/pma-qpcr-summary-stats.csv")

##################################################################################################################

#5. Performs statistical test for differences in abundance of iDNA and tDNA genes for all habitats, spring habitats and summer habitats.
#statistical difference between gene copies in treated and untreated samples

#remove anomolous result
#df2 = df2[df2$Sample != "S21.30",]

#is the data normally distributed?

# compute the difference
d <- with(df, 
          Gene.Copies.Per.mL.mg[Treatment == "tDNA"] - Gene.Copies.Per.mL.mg[Treatment == "iDNA"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value < 0.05 so the data is not normally distributed

#get the data
x = df[df$Treatment == "iDNA",]$Gene.Copies.Per.mL.mg
y = df[df$Treatment == "tDNA",]$Gene.Copies.Per.mL.mg

#parametric paired sample t test
#t.test(Gene.Copies.Per.mL.mg ~ Treatment, paired = TRUE, data = df)

#nonparametric paired sample Wilcoxon test
wilcox.test(x, y, paired = TRUE) #no significant difference between 

#what does this look like per group?
#this doesn't work for snow and spring ice because we don't have enough datapoints
# sn <- wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Habitat == "Snow",], paired = TRUE)
# sn  
# sp <- wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Habitat == "Spring Ice",], paired = TRUE)
# sp  
sm <- wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Habitat == "Summer Ice",], paired = TRUE)
sm  
cr <- wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Habitat == "Cryoconite",], paired = TRUE)
cr  

#what happens if we test spring and summer samples separately?
#nonparametric paired sample Wilcoxon test
wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Season == "Spring",], paired = TRUE)
wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Season == "Summer",], paired = TRUE)

#there is a significant difference for summer samples overall so which is it?
wilcox.test(Gene.Copies.Per.mL.mg ~ Treatment, data = df[df$Season == "Summer",], paired = TRUE, alternative="less")= 0
