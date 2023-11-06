#Clear workspace
rm(list=ls())
graphics.off()

#load packages
library(ggplot2)
library(PairedData)
library(dplyr)

#load data
df = read.csv("../data/meta_qpcr_res.csv", sep="\t")

#set the habitat order
df$Habitat = factor(df$Habitat, levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite", "Control"))

#plot
p = ggplot(df) +
  geom_point(aes(Sample.Name.To.Plot, log(Gene.Copies.Per.mL.g), color = Treatment, shape = Treatment), size=4) +
  ylab("log Gene Copies per mL/g")+
  xlab("Sample")+
  facet_grid(cols = vars(Habitat), scales = "free", space = "free", labeller = label_wrap_gen(width=10))+
  guides(color=guide_legend(override.aes=list(shape=c(21,24), alpha=.5, size=4)))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank())

#save the plot
pdf("../results/qPCR-plot.pdf", height=5, width=15)
print(p)
dev.off()

#another way to plot paired data

# Subset gene data without treatment
tdna <- subset(df,  Treatment == "tDNA", Gene.Copies.Per.mL.g,
                 drop = TRUE)
# subset gene data with treatment
idna <- subset(df,  Treatment == "iDNA", Gene.Copies.Per.mL.g,
                drop = TRUE)
# Plot paired data
pd <- paired(tdna, idna)
plot(pd, type = "profile") + theme_bw()


#basic stats
group_by(df, Treatment, Habitat) %>%
  summarise(
    count = n(),
    mean = mean(Gene.Copies.Per.mL.g, na.rm = TRUE),
    median = median(Gene.Copies.Per.mL.g, na.rm = TRUE),
    IQR = IQR(Gene.Copies.Per.mL.g, na.rm = TRUE)
  )

#statistical difference between gene copies in treated and untreated samples

#is the data normally distributed?

# compute the difference
d <- with(df, 
          Gene.Copies.Per.mL.g[Treatment == "tDNA"] - Gene.Copies.Per.mL.g[Treatment == "iDNA"])
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value < 0.05 so the data is not normally distributed

#parametric paired sample t test
#t.test(Gene.Copies.Per.mL.g ~ Treatment, paired = TRUE, data = df)

#nonparametric paired sample Wilcoxon test
res <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df, paired = TRUE)
res  

#nonparametic paired sample Wilcoxon test (is the median gene number after treatment less than 
#the median gene number before treatment?)
res <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df, paired = TRUE,
                   alternative="less")
res  

#these results suggest that overallthe PMA treated samples had a significantly higher number of gene copies per mL/g sample than PMA un-treated samples

#what does this look like per group?
sn <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Habitat == "Snow",], paired = TRUE)
sn  
sp <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Habitat == "Spring Ice",], paired = TRUE)
sp  
sm <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Habitat == "Summer Ice",], paired = TRUE)
sm  
cr <- wilcox.test(Gene.Copies.Per.mL.g ~ Treatment, data = df[df$Habitat == "Cryoconite",], paired = TRUE)
cr  

