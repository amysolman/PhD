## Method

# Physical and chemical characteristics of the data were investigated for significant differences between habitats: snow, spring ice, summer ice, cryoconite. Chemical data has been QC-ed according to the following method: <LOD and <LOQ samples were set at LOD/2 as advised by the following sources: 
# https://tidsskriftet.no/en/2022/09/medicine-and-numbers/measurements-below-detection-limit
# https://www.researchgate.net/post/How_should_one_treat_data_with_LOQ_values_during_statistcal_analysis
# https://www.genasis.cz/time-series/index.php?pg=home--treatment-of-values-under-loq
# 
# The mean value of control samples was then subtracted from each analyte. 
# Any resulting negative values were set to zero. 
# The data was visualized and Kruskal-Wallis and Wilcoxon tests were computed 
# to export significant differences between the means of our habitats. 
# These non-parametric test were used as our data was not normally distributed, 
# therefore did not meet the assumptions of a one-way ANOVA.


#clear workspace and load packages
rm(list=ls())
library(dplyr)
library(ggpubr)
library(car) #for levin's test of homogeneity of variances
library(broom)
library(tidyr)
library(stringr) #to split the column
source("00-solman-functions.R")
library(cowplot)
#remotes::install_github("coolbutuseless/ggpattern") #patterns for plots
library(ggpattern)
library(ggplot2)
#install.packages("ggbreak")
library(ggbreak)
library(ggfortify) #for function autoplot
library(patchwork)

#load data
chem.phy.ed = read.csv("../data/metadata.csv")
#remove controls
chem.phy.ed = chem.phy.ed[chem.phy.ed$Habitat != "Control",]

#PCA Plot of chemical characteristics
#https://www.geeksforgeeks.org/how-to-make-pca-plot-with-r/

#prep dataframe
#keep only major ions
keep = c("Habitat", "Location", "pH", "Conductivity_muS", "TN", "Cl", "Fe", "Na", "Mg", "K", "Ca")
pca.df = chem.phy.ed[, names(chem.phy.ed) %in% keep] #get edited chemical data with certain columns removed if they are not appropriate for this analysis 
pca.df$Location = paste(pca.df$Location, "Foxfonna")

#keep complete cases only
pca.df.comp = pca.df[complete.cases(pca.df),]

#principal component analysis
loc.pca = prcomp(pca.df.comp[,-c(1:2)], center=TRUE, scale. = TRUE)

#summary of the pca object
summary(loc.pca) #we can see 50% of the variance is captured in the first axis, while 14% is captured in the second axis.
pca.df.comp$Habitat = factor(pca.df.comp$Habitat,
                             levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
#plot the results
loc.pca.plot <- autoplot(loc.pca,
                         data = pca.df.comp,
                         fill = 'Habitat',
                         #color="Habitat",
                         shape='Location', 
                         size=8, 
                         alpha=0.7)+
  scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
  # #guides(fill=guide_legend(override.aes=list(colour=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme_bw()

loc.pca.plot

pdf("../results/pca-chem-variables.pdf", height=4)
print(loc.pca.plot)
dev.off()

biplot.loc.pca <- biplot(loc.pca)
biplot.loc.pca

#Plot chemical values

df = chem.phy.ed

for (k in 10:ncol(df)){ #for each item to analyse
  
  data = data.frame(Location = c(df$Location), Habitat = c(df$Habitat), Season = c(df$Season), X = df[,k])
  
  data$Habitat <- factor(data$Habitat,
                         levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  
  p = ggplot(data, aes(x=Habitat, y=X, fill=Location))+
    geom_boxplot()+
    facet_wrap(~Season, scales = "free_x")+
    theme_bw()+
    ylab(names(df)[k])
  
  print(p)
  
}


major.ions = c("pH", "Conductivity_muS", "DOC", "Cl", "NO3", "SO4", "Fe", "K",  "Na", "Mg")

plots2save = list()

for (m in 1:length(major.ions)){ #for each item to analyse
  
  data = data.frame(Location = c(df$Location), Habitat = c(df$Habitat), Season = c(df$Season), X = df[,major.ions[m]])
  data$Location = paste(data$Location, "Foxfonna")
  
  data$Habitat <- factor(data$Habitat,
                         levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  data$Habitat_Location = paste(data$Habitat, data$Location)
  data$Habitat_Location <- factor(data$Habitat_Location,
                                  levels=c("Snow Lower Foxfonna", "Snow Upper Foxfonna", "Spring Ice Lower Foxfonna", 
                                           "Spring Ice Upper Foxfonna",
                                           "Summer Ice Lower Foxfonna",
                                           "Summer Ice Upper Foxfonna",
                                           "Cryoconite Lower Foxfonna",
                                           "Cryoconite Upper Foxfonna"))
  
  #get units of concentration
  if (major.ions[m] == "pH"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle("pH")
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "Conductivity_muS"){
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Conducitivy" ~ mu * "S cm" ^-1))
    
    print(plots2save[[m]])
    
    
    
  } else if (major.ions[m] == "DOC"){
  
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("DOC mg L"^-1))+
      ylim(0,4)
    
    print(plots2save[[m]])
    
    
    
  } else if (major.ions[m] == "Cl"){
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Cl mg L"^-1))+
      ylim(0,4)
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "NO3"){
  
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("NO"[3] ~ "mg L"^-1))
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "SO4"){
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("SO"[4] ~ "mg L"^-1))
    
    print(plots2save[[m]])
    
  }else if (major.ions[m] == "Fe"){
    
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Fe" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "K"){
  
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("K" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "Na"){
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Na" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  } else if (major.ions[m] == "Mg"){
    
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Mg" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  }
  
  
}

print(plots2save[[1]])

multi.p = plot_grid(plots2save[[1]],
                    plots2save[[2]],
                    plots2save[[3]],
                    plots2save[[4]],
                    plots2save[[5]],
                    plots2save[[6]],
                    plots2save[[7]],
                    plots2save[[8]],
                    plots2save[[9]],
                    plots2save[[10]],
                    ncol=2, align=c("hv"))
multi.p


pdf("../results/major-ion-plots.pdf", width=8, height=14)
print(multi.p)
dev.off()

#Combine plots
pca.plot = loc.pca.plot + theme(legend.position=c(.75,.85), legend.box="vertical", legend.text = element_text(size=20),
                                legend.title = element_text(size=25), axis.text = element_text(size=15), axis.title = element_text(size=20),
)

# final.p = ggpubr::ggarrange(pca.plot, multi.p, ncol=2, labels=c("A", "B"), font.label = list(size=20))
# final.p

# final.p = ggpubr::ggarrange(pca.plot, multi.p, ncol=2, labels=c("Principal Components Analysis", "Values of Key Analytes"), font.label = list(size=20))
# final.p

#use patchwork to arrange the plots
#make titles and legend
title1 <- ggdraw() + draw_label("Principal Components Analysis", size=16, hjust=1.84, fontface = "bold")
title2 <- ggdraw() + draw_label("Key Analytes", size=16, hjust=4.45,fontface = "bold")
legend <- get_legend(loc.pca.plot + theme(legend.box="horizontal", legend.title = element_blank(), legend.text=element_text(size=16)) + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1), fill = guide_legend(nrow = 1, override.aes=list(shape=22))))

final.p = title1 / (loc.pca.plot + theme(legend.position = "none") )/ plot_spacer() / title2 / (plots2save[[1]] | plots2save[[2]] | plots2save[[3]] | plots2save[[4]]) / (plots2save[[5]] | plots2save[[6]] | plots2save[[7]] | plots2save[[8]]) / plot_spacer() /legend + plot_layout(heights = c(0.1, 1, 0.1, 0.1, 0.8, 0.8, 0.1, 0.1))
final.p


pdf("../results/chemical-analysis.pdf", width=13, height=13)
print(final.p)
dev.off()


##############################################################################################################
#MINOR IONS PLOT
minor.ions = c("Al", "Mn", "Li", "Ti", "V", "Co")

plots2save2 = list()

for (m in 1:length(minor.ions)){ #for each item to analyse
  
  data = data.frame(Location = c(df$Location), Habitat = c(df$Habitat), Season = c(df$Season), X = df[,minor.ions[m]])
  data$Location = paste(data$Location, "Foxfonna")
  
  data$Habitat <- factor(data$Habitat,
                         levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
  data$Habitat_Location = paste(data$Habitat, data$Location)
  data$Habitat_Location <- factor(data$Habitat_Location,
                                  levels=c("Snow Lower Foxfonna", "Snow Upper Foxfonna", "Spring Ice Lower Foxfonna", 
                                           "Spring Ice Upper Foxfonna",
                                           "Summer Ice Lower Foxfonna",
                                           "Summer Ice Upper Foxfonna",
                                           "Cryoconite Lower Foxfonna",
                                           "Cryoconite Upper Foxfonna"))
  
  #remove rows with NA in the X column
  data = data[!is.na(data$X),]
  
  
  
  if (minor.ions[m] == "Al"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Al" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  } else if (minor.ions[m] == "Mn"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Mn" ~ mu*"g L"^-1))
    
    print(plots2save[[m]])
    
  } else if (minor.ions[m] == "Li"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Li" ~ mu*"g L"^-1))+
      ylim(0, 3)
    
    print(plots2save[[m]])
    
  } else if (minor.ions[m] == "Ti"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9),
                  position=position_jitter(width = 0.2,
                                           height = 0.2))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Ti" ~ mu*"g L"^-1))+
      ylim(0, 3)
    
    print(plots2save[[m]])
    
  } else if (minor.ions[m] == "V"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("V" ~ mu*"g L"^-1)) +ylim(0, 0.15)
    
    print(plots2save[[m]])
    
  } else if (minor.ions[m] == "Co"){
    
    #plot
    plots2save[[m]] = ggplot(data, aes(x=Habitat_Location, y=X, fill=Habitat_Location))+
      #geom_boxplot_pattern(aes(pattern= Location, pattern_fill = Location),pattern_spacing = 0.03)+
      geom_boxplot()+
      theme_bw()+
      geom_jitter(aes(shape=Location, fill=Habitat_Location, size=3, alpha=.9))+
      scale_fill_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_colour_manual(values=c("#fa9f99", "#fa9f99", "#a4c64d", "#a4c64d", "#4dd2d6", "#4dd2d6", "#d8a4ff", "#d8a4ff"))+
      scale_shape_manual(values=c(21,24,21,24,21,24,21,24))+
      theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      ggtitle(expression("Co" ~ mu*"g L"^-1))+
      ylim(0, 0.15)
    
    print(plots2save[[m]])
    
  }
  
}


multi.p2 = plot_grid(plots2save[[1]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     plots2save[[2]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     plots2save[[3]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     plots2save[[4]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     plots2save[[5]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     plots2save[[6]] + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y=element_blank(),
                                             axis.text.x = element_blank(), axis.ticks.x=element_blank()),
                     ncol=2, align=c("hv"))
multi.p2

#get legend
legend <- get_legend(
  loc.pca.plot + theme(legend.box="verticle", legend.margin=margin(), legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=16)) 
  + guides(fill = guide_legend(override.aes=list(shape=22)))
  
)

minor.ion.p = plot_grid(multi.p2, legend, ncol=1, rel_heights = c(1, 0.1))

pdf("../results/minor-ion-plots.pdf", width=10, height=10)
print(minor.ion.p)
dev.off()


##############################################################################################################


ree = c("Nd", "Ce", "Pr", "Sc", "Sm",
        "Eu", "Tm", "Pm", "La", "Gd", 
        "Dy")

#elements on the x axis, concentrations on y axis, line graph  grouped/coloured by habitat

#get ree data
ree.data = df[,names(df) %in% c("Habitat",  "Location", ree)]

#turn into long data
ree.long = gather(ree.data, Analyte, Concentration, Sc:Tm)

# #data summary function
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# #summarise data
# df3 <- data_summary(ree.long, varname="Concentration", 
#                     groupnames=c("Habitat", "Location", "Analyte"))
# head(df3)
# 
# df3$Habitat <- factor(df3$Habitat,
#                       levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
# 
# df3[df3 == "Lower"] <- "Lower Foxfonna"
# df3[df3 == "Upper"] <- "Upper Foxfonna"
# 
# ree.plot = ggplot(df3, aes(x=Analyte, y=Concentration, group=Habitat, color=Habitat)) + 
#   geom_errorbar(aes(ymin=Concentration-sd, ymax=Concentration+sd), width=.1) +
#   geom_line() + geom_point()+
#   # scale_color_brewer(palette="Paired")+
#   theme_bw()+
#   facet_grid(rows=vars(Location))+
#   scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
#   theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
#   ylab("REE Concentration (ng L-1)")
# ree.plot
# 
# pdf("../results/ree-plot.pdf", width=8, height=4)
# print(ree.plot)
# dev.off()
# 
# 
# #ree.boxplot
# ree.long$Habitat <- factor(ree.long$Habitat,
#                            levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
# 
# ree.box = ggplot(ree.long, aes(x=Analyte, y=Concentration, fill=Habitat))+
#   geom_boxplot()+
#   theme_bw()+
#   scale_colour_manual(values=c("#fa9f99", "#a4c64d", "#4dd2d6", "#d8a4ff"))+
#   theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size=25))+
#   ylab(expression("REE Concentration (ng L"^-1*")"))+
#   scale_y_break(c(55,70))+
#   facet_grid(cols=vars(Location))
# 
# ree.box
# 
# #save
# pdf.options(reset = TRUE, onefile = FALSE)
# pdf("../results/ree-boxplot.pdf", width=20, height=8)
# print(ree.box)
# dev.off()

#ree.boxplot
ree.long$Habitat <- factor(ree.long$Habitat,
                           levels=c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))

ree.box = ggplot(ree.long, aes(x=Analyte, y=Concentration, fill=Habitat))+
  geom_boxplot()+
  facet_grid(rows=("Location"))+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank())+
  ylab(expression("REE Concentration (ng L"^-1*")"))
#ylim(0,40)
ree.box

#save
pdf.options(reset = TRUE, onefile = FALSE)
pdf("../results/ree-boxplot.pdf", width=20, height=8)
print(ree.box)
dev.off()


#Combine REE and minor elements plots

#make titles and legend
title1 <- ggdraw() + draw_label("Minor Ions", size=16, hjust=6.2, fontface = "bold")
title2 <- ggdraw() + draw_label("Rare Earth Elements", size=16, hjust=3.2,fontface = "bold")

legend <- get_legend(loc.pca.plot + theme(legend.box="horizontal", legend.title = element_blank(), legend.text=element_text(size=16)) + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1), fill = guide_legend(nrow = 1, override.aes=list(shape=22))))

combi.plot = title1 / (plots2save[[1]] | plots2save[[2]]) / (plots2save[[3]] | plots2save[[4]] ) / (plots2save[[5]] | plots2save[[6]]) / plot_spacer() / title2 / (ree.box + theme(legend.position = "none", axis.text = element_text(size=15), strip.text = element_text(size=15))) / plot_spacer() / legend + plot_layout(heights = c(0.1, 1, 1, 1, 0.1, 0.1, 3, 0.1, 0.1))
combi.plot

pdf("../results/chem-combi-plot.pdf", width=15, height=25)
print(combi.plot)
dev.off()
