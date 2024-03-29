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

#load data

chem = read.csv("../data/meta_chem.csv")
phy = read.csv("../data/meta_phy.csv")
cont = read.csv("../data/meta_cont.csv")

#Prep Data for Analysis
#make sure values are set to <LOQ and <LOD

chem.mod.trim = chem

chem.mod.trim[is.na(chem.mod.trim)] <- "NA"

for (i in 3:ncol(chem.mod.trim)){ #for each column
  
  for (j in 6:length(chem.mod.trim[,i])){ #for each row
    
    if (chem.mod.trim$Season[j] == "Spring"){ #if the season is spring
      
      if (chem.mod.trim[j,i] == "<LOD" | chem.mod.trim[j,i] == "<LOQ"){ #if the value is <LOD or <LOQ
        
        chem.mod.trim[j,i] = (as.numeric(chem.mod.trim[2,i])/2)  #replace with corresponding LOD/2
        
      }
      
    } else if (chem.mod.trim$Season[j] == "Summer"){ #if the season is summer
      
      if (chem.mod.trim[j,i] == "<LOD" | chem.mod.trim[j,i] == "<LOQ"){ #if the value is <LOD or <LOQ
        
        chem.mod.trim[j,i] = (as.numeric(chem.mod.trim[4,i])/2)  #replace with corresponding LOD/2
        
        
        
      }
    }
  }
  
}

#make sure everything is numeric
chem.mod.trim.num = chem.mod.trim
chem.mod.trim.num[,c(3:ncol(chem.mod.trim.num))] <- mutate_all(chem.mod.trim.num[,c(3:ncol(chem.mod.trim.num))], function(x) as.numeric(as.character(x)))
is.numeric(chem.mod.trim.num[3,3]) #should equal true

#separate sample data from blanks/QC/unit values
chem.samp = chem.mod.trim.num[-c(1:5),]
names(chem.samp) = c("SampleID", names(chem.samp)[2:ncol(chem.samp)])

#separate spring and summer data
sum.chem = chem.samp[chem.samp$Season == "Summer",]
spr.chem = chem.samp[chem.samp$Season == "Spring",]

#get our control data
sum.cont = sum.chem[c(1:4),]
spr.cont = spr.chem[c(1:7),]

#subtract the mean value of control samples from each analyte
sum.chem.no.b = sum.chem[-c(1:4),]

for (j in 3:ncol(sum.chem.no.b)){ #for each column in summer chemical data
  c = sum.chem.no.b[,j] #get the column data
  c2 = c-mean(sum.cont[,j], na.rm = TRUE)
  sum.chem.no.b[,j] <- c2
}

spr.chem.no.b = spr.chem[-c(1:7),]

for (j in 3:ncol(spr.chem.no.b)){
  c = spr.chem.no.b[,j] #get the column data
  c2 = c-mean(spr.cont[,j], na.rm = TRUE)
  spr.chem.no.b[,j] <- c2
}

#recombin spring and summer data
sp.sum = rbind(spr.chem.no.b, sum.chem.no.b)

#any resulting negative values were set to zero.
sp.sum.0 = sp.sum
sp.sum.0[sp.sum.0 < 0] <- 0 

#combine physical and chemical analytes
chem.phy = cbind(phy, sp.sum.0[,-c(1:2)])

#columns to remove
rm <- c("DiamNS", "DiamEW", "SnowDepth1", "SnowDepth2", "SnowDepth3")
chem.phy.ed <- chem.phy[,!names(chem.phy) %in% rm]

#one of our results is anomalous - we have a conductivity result of 183, while the next lowest is 42 so I'm going to remove this value
chem.phy.ed = chem.phy.ed %>% 
  mutate(Conductivity_muS = replace(Conductivity_muS, Conductivity_muS > 150, NA))

#add in control data to edited chemical data
chem.phy.ed.cont = rbind(cont, chem.phy.ed)

#export the final dataframe
write.csv(chem.phy.ed.cont, "../data/metadata.csv", row.names = FALSE)

#Mean and standard deviation for each analyte for each group

df = chem.phy.ed

#set zero to NA so it won't be included in calculating minimums
df[df == 0] <- NA

df.sum = df[,c(4,10:ncol(df))] %>%
  dplyr::group_by(Habitat) %>%
  #filter() #filter out zeros??
  dplyr::summarise(across(everything(), .f = list(mean = mean,  sd = sd, max = max, min=min), na.rm = TRUE)) 

#get minimums that aren't zero!!!!!!!!!!

#into long format
data_long <- gather(df.sum, Var, Measurement, pH_mean:Yb_min, factor_key=TRUE)

#split the var column
data_long[c('Var', 'Metric')] <- str_split_fixed(data_long$Var, '_', 2)

final.df = data.frame(group1=data_long$Habitat, var=data_long$Var, Metric=data_long$Metric, Value=round(data_long$Measurement, 5))

#make sure our units are correct
final.df[final.df == "muS_max"] <- "max"
final.df[final.df == "muS_mean"] <- "mean"
final.df[final.df == "muS_min"] <- "min"
final.df[final.df == "muS_sd"] <- "sd"
final.df[final.df == "cm_max"] <- "max"
final.df[final.df == "cm_mean"] <- "mean"
final.df[final.df == "cm_min"] <- "min"
final.df[final.df == "cm_sd"] <- "sd"
# df.save[df.save == "Conductivity_muS"] <- "Conductivity"

#respread the data
data_wide <- spread(final.df, Metric, Value)

final.df2 = final.df
final.df2$Metric2 = paste(final.df2$group1, final.df2$Metric)

#respread the data
final.df3 = data.frame(Metric=final.df2$Metric2, Var=final.df2$var, Value=round(final.df2$Value, 2))

#spread the data, make Var a factor first
data_wide2 = final.df3 %>% 
  mutate(Var = factor(Var, levels = unique(Var))) %>%
  spread(Metric, Value)

#change inf + nan values to NA
data_wide2[data_wide2 == -Inf] <- NA
data_wide2[data_wide2 == "Inf"] <- NA
data_wide2[data_wide2 == NaN] <- NA
data_wide2[data_wide2 == "NaN"] <- NA

# data_wide2 <- spread(final.df3, Metric, Value)
data_wide3 = data.frame(cbind(data_wide2$`Snow max`,
                              data_wide2$`Snow min`,
                              data_wide2$`Snow mean`,
                              data_wide2$`Snow sd`,
                              data_wide2$`Spring Ice max`,
                              data_wide2$`Spring Ice min`,
                              data_wide2$`Spring Ice mean`,
                              data_wide2$`Spring Ice sd`,
                              data_wide2$`Summer Ice max`,
                              data_wide2$`Summer Ice min`,
                              data_wide2$`Summer Ice mean`,
                              data_wide2$`Summer Ice sd`,
                              data_wide2$`Cryoconite max`,
                              data_wide2$`Cryoconite min`,
                              data_wide2$`Cryoconite mean`,
                              data_wide2$`Cryoconite sd`))
names(data_wide3) = c("Snow Max", "Snow Min", "Snow Mean", "Snow SD",
                      "Spring Ice Max", "Spring Ice Min", "Spring Ice Mean", "Spring Ice SD",
                      "Summer Ice Max", "Summer Ice Min", "Summer Ice Mean", "Summer Ice SD",
                      "Cryoconite Max", "Cryoconite Min", "Cryoconite Mean", "Cryoconite SD")
#add analytes as row names
rownames(data_wide3) = data_wide2$Var

#Add LOD values and units
lod = chem.mod.trim[c(1,2,4),-c(1,2)] #get data
# lod = data.frame(cbind(paste(lod$X, lod$Season), lod[,3:ncol(lod)])) #add LOD Spring, LOD Summer columns
# names(lod) = c("")
# lod.l = gather(lod, Analyte, Value, TC:Yb) #make data long

lod2 = data.frame(Unit = as.character(lod[1,]), "LOD Spring" = as.numeric(lod[2,]), "LOD Summer" = as.numeric(lod[3,]), Analyte = names(lod))
#round our LOD values
lod2$LOD.Spring = round(lod2$LOD.Spring, 4)
lod2$LOD.Summer = round(lod2$LOD.Summer, 4)

#combine with data_wide3
#add an analyte column
data_wide4 = data.frame(cbind(rownames(data_wide3), data_wide3))
#name that column
names(data_wide4) = c("Analyte", names(data_wide4)[2:length(names(data_wide4))])
#merge dataframes
data = left_join(data_wide4, lod2, by="Analyte")

#rearrange columns
data2save = data.frame(cbind(data$Analyte, data$Unit, data$LOD.Spring, data$LOD.Summer, data[,2:17]))

#replace values with <LOD where appropriate
for (i in 6:nrow(data2save)){
  
  if (is.na(data2save[i,5]) == FALSE & data2save[i,5] < data2save[i,3]){ #if snow max is less than spring LOD
    data2save[i,5] <- "<LOD"
    data2save[i,8] <- "NA" #make SD NA
  } 
  
  if (is.na(data2save[i,6]) == FALSE & data2save[i,6] < data2save[i,3]){ #if snow min is less than spring LOD
    data2save[i,6] <- "<LOD"
  }
  
  if (is.na(data2save[i,7]) == FALSE & data2save[i,7] < data2save[i,3]){ #if snow mean is less than spring LOD
    data2save[i,7] <- "<LOD"
  }
  
  if (is.na(data2save[i,9]) == FALSE & data2save[i,9] < data2save[i,3]){ #if spring max is less than spring LOD
    data2save[i,9] <- "<LOD"
    data2save[i,12] <- "NA" #make SD NA
  } 
  
  if (is.na(data2save[i,10]) == FALSE & data2save[i,10] < data2save[i,3]){ #if spring min is less than spring LOD
    data2save[i,10] <- "<LOD"
  }
  
  if (is.na(data2save[i,11]) == FALSE & data2save[i,11] < data2save[i,3]){ #if spring mean is less than spring LOD
    data2save[i,11] <- "<LOD"
  }
  
  if (is.na(data2save[i,13]) == FALSE & data2save[i,13] < data2save[i,4]){ #if summer max is less than summer LOD
    data2save[i,13] <- "<LOD"
    data2save[i,16] <- "NA" #make SD NA
  } 
  
  if (is.na(data2save[i,14]) == FALSE & data2save[i,14] < data2save[i,4]){ #if summer min is less than summer LOD
    data2save[i,14] <- "<LOD"
  }
  
  if (is.na(data2save[i,15]) == FALSE & data2save[i,15] < data2save[i,4]){ #if summer mean is less than summer LOD
    data2save[i,15] <- "<LOD"
  }
  
  if (is.na(data2save[i,17]) == FALSE & data2save[i,17] < data2save[i,4]){ #if cryo max is less than summer LOD
    data2save[i,17] <- "<LOD"
    data2save[i,20] <- "NA" #make SD NA
  } 
  
  if (is.na(data2save[i,18]) == FALSE & data2save[i,18] < data2save[i,4]){ #if cryo min is less than summer LOD
    data2save[i,18] <- "<LOD"
  }
  if (is.na(data2save[i,19]) == FALSE & data2save[i,19] < data2save[i,4]){ #if cryo mean is less than summer LOD
    data2save[i,19] <- "<LOD"
  }
  
}

#order for table rows
#this is the order we want our analytes in
order = c("pH", "Conductivity", "Temp", "Area", "Depth", "TC", "DOC", "TN", "NO2", "NO3", "PO4", "SO4", "Cl", "Br", 
          "Fl", "Fe", "Na", "Mg", "K", "Ca", "Li", "Al", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn", "Y", "Rb", 
          "Sr", "Mo", "Ag", "Ba", "Lu", "Pb", "Zr", "Cd", "As", "Se", "Nb", "Sn", "Cs", "Hf", "Ta", "Re", "U", "La", "Ce", 
          "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Yb", "Ru", "Sc", "Tm")
#get the index of the order we want against the analyte rows we want to reorder
idx <- match(order, data2save$data.Analyte)
#reorder based on those indexes
data2save_ordered  <- data2save[idx,]


#save dataframe
write.csv(data2save_ordered, "../results/physical-data-stats.csv")
