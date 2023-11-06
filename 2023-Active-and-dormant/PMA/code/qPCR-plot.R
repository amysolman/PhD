#Clear workspace
rm(list=ls())
graphics.off()

#load packages
library(ggplot2)
library(PairedData)
library(dplyr)

#load data
df = read.csv("../data/meta_qpcr_res.csv", sep="\t")
df2 = read.csv("../data/qPCR-mean.csv", sep="\t")

#set the habitat order
df$Habitat = factor(df$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite", "Control"))
df2$Habitat = factor(df2$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite", "Control"))

#plot
p = ggplot(df) +
  geom_point(aes(Sample.Name.To.Plot, log(Gene.Copies.Per.mL.g), color = Treatment, shape = Treatment, stroke=2), size=6) +
  ylab("log Gene Copies per mL/g")+
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
pdf("../results/qPCR-plot.pdf", height=5, width=15)
print(p)
dev.off()

#second plot
p2 = ggplot(df2) +
  geom_point(aes(Sample, log(Gene.Copies.Per.mL.g), color = Treatment, shape = Treatment, stroke=2), size=6) +
  ylab("log Gene Copies per mL/g")+
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
pdf("../results/qPCR-plot2.pdf", height=5, width=15)
print(p2)
dev.off()

#another way to plot paired data

# Subset gene data without treatment
tdna <- subset(df2,  Treatment == "tDNA", Gene.Copies.Per.mL.g,
                 drop = TRUE)
# subset gene data with treatment
idna <- subset(df2,  Treatment == "iDNA", Gene.Copies.Per.mL.g,
                drop = TRUE)
# Plot paired data
pd <- paired(tdna, idna)
plot(pd, type = "profile") + theme_bw()


#basic stats
group_by(df2, Treatment, Habitat) %>%
  summarise(
    count = n(),
    mean = mean(Gene.Copies.Per.mL.g, na.rm = TRUE),
    median = median(Gene.Copies.Per.mL.g, na.rm = TRUE),
    IQR = IQR(Gene.Copies.Per.mL.g, na.rm = TRUE)
  )

#statistical difference between gene copies in treated and untreated samples

#is the data normally distributed?

# compute the difference
d <- with(df2, 
          Gene.Copies.Per.mL.g[Treatment == "tDNA"] - Gene.Copies.Per.mL.g[Treatment == "iDNA"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value < 0.05 so the data is not normally distributed

#parametric paired sample t test
#t.test(Gene.Copies.Per.mL.g ~ Treatment, paired = TRUE, data = df)

#nonparametric paired sample Wilcoxon test
wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df2, paired = TRUE) #no significant difference between 


#nonparametic paired sample Wilcoxon test (is the median gene number in samples treated with PMA (iDNA) 
#less than untreated samples (tDNA)?)
# x = df2[df2$Treatment == "iDNA",]
# y = df2[df2$Treatment == "tDNA",]
# wilcox.test(x$Gene.Copies.Per.mL.g, y$Gene.Copies.Per.mL.g, paired = TRUE,
#                    alternative="less")

#nonparametic paired sample Wilcoxon test (is the median gene number in samples treated with PMA (iDNA) 
#greater than untreated samples (tDNA)?)
# wilcox.test(x$Gene.Copies.Per.mL.g, y$Gene.Copies.Per.mL.g, paired = TRUE,
#             alternative="greater")

#these results suggest that overall the PMA treated samples had a significantly higher number of gene copies per mL/g sample than PMA un-treated samples

#what does this look like per group?
sn <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df2[df2$Habitat == "Snow",], paired = TRUE)
sn  
sp <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df2[df2$Habitat == "Spring Ice",], paired = TRUE)
sp  
sm <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df2[df2$Habitat == "Summer Ice",], paired = TRUE)
sm  
cr <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df2[df2$Habitat == "Cryoconite",], paired = TRUE)
cr  

#what happens if we test spring and summer samples separately?
#nonparametric paired sample Wilcoxon test
wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Season == "Spring",], paired = TRUE)
wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Season == "Summer",], paired = TRUE)

#Spring
x = df[df$Season == "Spring" & df$Treatment == "iDNA",]
y = df[df$Season == "Spring" & df$Treatment == "tDNA",]
wilcox.test(x$Gene.Copies.Per.mL.g, y$Gene.Copies.Per.mL.g, paired = TRUE)

#Summer
x = df[df$Season == "Summer" & df$Treatment == "iDNA",]
y = df[df$Season == "Summer" & df$Treatment == "tDNA",]
wilcox.test(x$Gene.Copies.Per.mL.g, y$Gene.Copies.Per.mL.g, paired = TRUE, alternative = "less")

#test for differences between individual pairs
data = df[df$Sample == "S21.29",]
x = data[data$Treatment == "iDNA",]
y = data[data$Treatment == "tDNA",]
wilcox.test(x$Gene.Copies.Per.mL.g, y$Gene.Copies.Per.mL.g, paired = TRUE)