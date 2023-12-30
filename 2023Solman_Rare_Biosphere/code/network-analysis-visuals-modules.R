#Import Gephi Data Table

rm(list=ls())

library(tibble)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(phyloseq)
#install.packages("ggtern")
library(ggtern)
source("00-solman-functions.R")
source("00-abundance-functions.R")
library(RColorBrewer) #for our pie chart colours

#gephi table
gephi.tab.ant = read.csv("../results/antarctic-gephi-network-table.csv")
gephi.tab.arc = read.csv("../results/arctic-gephi-network-table.csv")

#keystone taxa
keystone.ant = read.csv("../results/Antarctic-keystone-taxa.csv")
keystone.arc = read.csv("../results/Arctic-keystone-taxa.csv")

#Load data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) >= 1, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun-var-trans.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int-var-trans.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare-var-trans.rds")
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun-var-trans.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int-var-trans.rds")
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare-var-trans.rds")

pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds")
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds")
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds")
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds")
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds")
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds")

#Get abundance dataframes

#Abundance df function
abun_dataframe_func <- function(phylo.abun, phylo.int, phylo.rare){
  
  #Get data frame with classifications
  AT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.abun))), 
                       Abundance = rep("AT", length(rownames(data.frame(otu_table(phylo.abun))))), 
                       Colour = "#ef7a76")
  
  MT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.int))), Abundance = rep("MT", length(rownames(data.frame(otu_table(phylo.int))))), Colour = "#74a9d8")
  RT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.rare))), Abundance = rep("RT", length(rownames(data.frame(otu_table(phylo.rare))))), Colour="#82ca81")
  
  
  abun.df = rbind(AT_DATA, MT_DATA, RT_DATA)
  
  return(abun.df)
  
}

pro.ant.df = abun_dataframe_func(pro.ant.abun, pro.ant.int, pro.ant.rare)
pro.arc.df = abun_dataframe_func(pro.arc.abun, pro.arc.int, pro.arc.rare)
euk.ant.df = abun_dataframe_func(euk.ant.abun, euk.ant.int, euk.ant.rare)
euk.arc.df = abun_dataframe_func(euk.arc.abun, euk.arc.int, euk.arc.rare)

#merge
ant.df = rbind(pro.ant.df, euk.ant.df)
arc.df = rbind(pro.arc.df, euk.arc.df)


#Add Abundance data to ASV ID in gephi data frame
abun = ant.df[ant.df$ID %in% gephi.tab.ant$name, ]
#this gives the index order in which to put our M1 data frame to match the gephi.tab order
order = match(gephi.tab.ant$name, abun$ID)
#so now we just reorder it!
reordered = abun[order,]

#now append to our gephi data table
gephi.tab.ant$Abundance = reordered$Abundance
gephi.tab.ant$Colour = reordered$Colour

#Add Abundance data to ASV ID in gephi data frame
abun = arc.df[arc.df$ID %in% gephi.tab.arc$name, ]
#this gives the index order in which to put our M1 data frame to match the gephi.tab order
order = match(gephi.tab.arc$name, abun$ID)
#so now we just reorder it!
reordered = abun[order,]
#now append to our gephi data table
gephi.tab.arc$Abundance = reordered$Abundance
gephi.tab.arc$Colour = reordered$Colour

#get taxonomic data

pro.ant.tax <- data.frame(tax_table(pro.ant))
pro.tax.tab <- tibble::rownames_to_column(pro.ant.tax, "ID")
euk.ant.tax <- data.frame(tax_table(euk.ant))
euk.tax.tab <- tibble::rownames_to_column(euk.ant.tax, "ID")
tax.tab = rbind(pro.tax.tab, euk.tax.tab)
#Append to gephi table
#subset tax data by network IDs
tax.trim = tax.tab[tax.tab$ID %in% gephi.tab.ant$name, ]
#this gives the index order in which to put our data frame to match the gephi.tab order
tax.order = match(gephi.tab.ant$name, tax.trim$ID)
#so now we just reorder it!
tax.reordered = tax.trim[tax.order,]
#now append to our gephi data table
final.df.ant = cbind(gephi.tab.ant, tax.reordered)
#Save CSV file without row names otherwise Gephi won't recognise the data 
#export
write.csv(final.df.ant, "../results/antarctic-gephi-data-table-edit.csv", row.names = FALSE)

pro.arc.tax <- data.frame(tax_table(pro.arc))
pro.tax.tab <- tibble::rownames_to_column(pro.arc.tax, "ID")
euk.arc.tax <- data.frame(tax_table(euk.arc))
euk.tax.tab <- tibble::rownames_to_column(euk.arc.tax, "ID")
tax.tab = rbind(pro.tax.tab, euk.tax.tab)
#Append to gephi table
#subset tax data by network IDs
tax.trim = tax.tab[tax.tab$ID %in% gephi.tab.arc$name, ]
#this gives the index order in which to put our data frame to match the gephi.tab order
tax.order = match(gephi.tab.arc$name, tax.trim$ID)
#so now we just reorder it!
tax.reordered = tax.trim[tax.order,]
#now append to our gephi data table
final.df.arc = cbind(gephi.tab.arc, tax.reordered)
#Save CSV file without row names otherwise Gephi won't recognise the data 
#export
write.csv(final.df.arc, "../results/arctic-gephi-data-table-edit.csv", row.names = FALSE)

#Modular Analysis

###Pie Charts

#for each module I need the number of abundant, intermediate and rare taxa
#sort(table(final.df.arc$modularity_class)) 

#here specify the numbers of the 6 largest modules
ant.order = c(as.numeric(names(sort(table(gephi.tab.ant$modularity_class), decreasing = TRUE)[1:6])))
arc.order = c(as.numeric(names(sort(table(gephi.tab.arc$modularity_class), decreasing = TRUE)[1:6])))

# final.df = final.df.ant
# mod = ant.order[[1]]
# name = "Module I"

pie_data <- function(final.df, mod, name, pole){
  
  counts = final.df[final.df$modularity_class == mod,]
  
  x = data.frame(table(counts$Abundance))
  levels(x$Var1) <- c(levels(x$Var1), "Abundant")
  levels(x$Var1) <- c(levels(x$Var1), "Intermediate")
  levels(x$Var1) <- c(levels(x$Var1), "Rare")
  x$Var1[x$Var1 == 'AT'] <- 'Abundant'
  x$Var1[x$Var1 == 'MT'] <- 'Intermediate'
  x$Var1[x$Var1 == 'RT'] <- 'Rare'
  x$Module = name
  x$Pole = pole
  x$Freq = x$Freq/sum(x$Freq)
  
  return(x)
}

ant.pie1 = pie_data(final.df.ant, ant.order[[1]], "Module I", "Antarctic")
ant.pie2 = pie_data(final.df.ant, ant.order[[2]], "Module II", "Antarctic")
ant.pie3 = pie_data(final.df.ant, ant.order[[3]], "Module III", "Antarctic")
ant.pie4 = pie_data(final.df.ant, ant.order[[4]], "Module IV", "Antarctic")
ant.pie5 = pie_data(final.df.ant, ant.order[[5]], "Module V", "Antarctic")
ant.pie6 = pie_data(final.df.ant, ant.order[[6]], "Module VI", "Antarctic")
arc.pie1 = pie_data(final.df.arc, arc.order[[1]], "Module I", "Arctic")
arc.pie2 = pie_data(final.df.arc, arc.order[[2]], "Module II", "Arctic")
arc.pie3 = pie_data(final.df.arc, arc.order[[3]], "Module III", "Arctic")
arc.pie4 = pie_data(final.df.arc, arc.order[[4]], "Module IV", "Arctic")
arc.pie5 = pie_data(final.df.arc, arc.order[[5]], "Module V", "Arctic")
arc.pie6 = pie_data(final.df.arc, arc.order[[6]], "Module VI", "Arctic")

pie.df = rbind(ant.pie1, ant.pie2, ant.pie3, ant.pie4, ant.pie5, ant.pie6, arc.pie1, arc.pie2, arc.pie3, arc.pie4, arc.pie5, arc.pie6)

#specify the subcommunity order
pie.df$Var1 = factor(pie.df$Var1, levels = c("Rare", "Intermediate", "Abundant"))

pie.p = ggplot(pie.df, aes(x="", y=Freq, fill=Var1)) +geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +theme_void()+
  facet_wrap(~Pole + Module, nrow = 2)+
  scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
  theme(legend.position = "bottom", legend.title = element_blank())

pdf("../results/network-module-pies.pdf")
print(pie.p)
dev.off()

#pie chart of module taxonomy
# 
# final.df = final.df.ant
# phylo1 = pro.ant
# phylo2 = euk.ant
# mod = ant.order[[1]]
# name = "Module I"

tax_pie_func <- function(phylo1, phylo2, final.df, mod, name, pole){
  
  #add the total counts for each ASV to the dataframe
  ps = merge_my_phylo_same_pole(phylo1, phylo2)
  counts = data.frame(otu_table(ps))
  counts = counts[rownames(counts) %in% final.df$name,]
  data = data.frame(name=rownames(counts), num = rowSums(counts))
  #join with module data
  data.to.plot = left_join(data, final.df, by="name")
  #subset by module
  data.to.plot2 = data.to.plot[data.to.plot$modularity_class == mod,]
  richness = length(unique(data.to.plot2$Class))
  #convert num into relative abundance
  data.to.plot2$rel = data.to.plot2$num/sum(data.to.plot2$num)
  plot.df = data.frame(Tax=data.to.plot2$Class, Abun=data.to.plot2$rel)
  #get values for each group
  out = plot.df %>%
    dplyr::group_by(Tax) %>%
    dplyr::summarise(sum = sum(Abun))
  
  #combine anything less than 0.1 - 10%
  out.keep = out[out$sum > 0.1,]
  out.keep = rbind(out.keep, data.frame(Tax="Other", sum=1-sum(out.keep$sum)))
  out.keep$Module = name
  out.keep$Pole = pole
  
  return(res.list = list(out.keep, richness))
  
}

ant.pie1 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[1]], "Module I", "Antarctic")
ant.pie2 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[2]], "Module II", "Antarctic")
ant.pie3 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[3]], "Module III", "Antarctic")
ant.pie4 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[4]], "Module IV", "Antarctic")
ant.pie5 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[5]], "Module V", "Antarctic")
ant.pie6 = tax_pie_func(pro.ant, euk.ant, final.df.ant, ant.order[[6]], "Module VI", "Antarctic")
arc.pie1 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[1]], "Module I", "Arctic")
arc.pie2 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[2]], "Module II", "Arctic")
arc.pie3 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[3]], "Module III", "Arctic")
arc.pie4 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[4]], "Module IV", "Arctic")
arc.pie5 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[5]], "Module V", "Arctic")
arc.pie6 = tax_pie_func(pro.arc, euk.arc, final.df.arc, arc.order[[6]], "Module VI", "Arctic")

pie.df = rbind(ant.pie1[[1]], ant.pie2[[1]], ant.pie3[[1]], ant.pie4[[1]], ant.pie5[[1]], ant.pie6[[1]], arc.pie1[[1]], arc.pie2[[1]], arc.pie3[[1]], arc.pie4[[1]], arc.pie5[[1]], arc.pie6[[1]])

rich.df = data.frame(Module=rep(c("Module I", "Module II", "Module III", "Module IV", "Module V", "Module VI"), 2), Pole = c(rep("Antarctic", 6), rep("Arctic", 6)), Richness = c(ant.pie1[[2]], ant.pie2[[2]], ant.pie3[[2]], ant.pie4[[2]], ant.pie5[[2]], ant.pie6[[2]], arc.pie1[[2]], arc.pie2[[2]], arc.pie3[[2]], arc.pie4[[2]], arc.pie5[[2]], arc.pie6[[2]]))
length(unique(pie.df$Tax))

getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(pie.df, aes(x="", y=sum, fill=Tax)) +geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +theme_void()+
  facet_wrap(~Module + Pole)+
  #  scale_fill_brewer(palette="Set1")+
  scale_fill_manual(values = getPalette(length(unique(pie.df$Tax))))+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size=4))

#Below is additional module analysis code that is not necessary for this paper

#In which modules do keystone species occur and what is their taxonomy
# key = keystone.ant
# final.df = final.df.ant

# ant.order = c(4, 15, 39, 34, 1, 37)
# arc.order = c(4, 9, 11, 1, 3, 13)

# old.ant <- ant.order
# new.ant <- c(1,2,3,4,5,6)
# old.arc <- arc.order
# new.arc <- c(1,2,3,4,5,6)
# 
# final.df.ant$modularity_class[final.df.ant$modularity_class %in% old.ant] <- new.ant[match(final.df.ant$modularity_class, old.ant, nomatch = 0)]
# final.df.arc$modularity_class[final.df.arc$modularity_class %in% old.arc] <- new.arc[match(final.df.arc$modularity_class, old.arc, nomatch = 0)]
# 
# df1 = final.df.ant[final.df.ant$name %in% keystone.ant$X,]
# df2 = final.df.arc[final.df.arc$name %in% keystone.arc$X,]
# 
# key.mod = data.frame(Pole = c(rep("Antarctic", nrow(df1)), rep("Arctic", nrow(df2))), Module=c(df1$modularity_class, df2$modularity_class), Abundance = c(df1$Abundance, df2$Abundance), Phylum = c(df1$Phylum, df2$Phylum), Class=c(df1$Class, df2$Class), Order=c(df1$Order, df2$Order), Family=c(df1$Family, df2$Family))
# 
# 
# # df = final.df.ant
# # phylo1 =  pro.ant
# # phylo2 = euk.ant
# # key = keystone.ant
# 
# module_analysis <- function(df, phylo1, phylo2, key){
#   
#   
#   #merge phyloseq objects
#   merge.ps = merge_my_phylo_same_pole(phylo1, phylo2)
#   
#   #get ASV table
#   counts = data.frame(otu_table(merge.ps))
#   
#   #get the samples each ASV is found in
#   out <- simplify2array(apply(counts[1:ncol(counts)], 1, function(x) paste(names(counts[1:ncol(counts)])[x != 0], collapse = " ")))
#   y = data.frame(out)
#   z = data.frame(str_split_fixed(y$out, " ", ncol(counts)))
#   z[z==""] = NA
#   z1 = cbind(rownames(counts), z)
#   names(z1) = c("ASV", seq(from = 1, to =ncol(counts)))
#   
#   #reshape the data
#   melt.me = melt(z1, "ASV")
#   #remove rows with missing values
#   melt.me.keep = melt.me[complete.cases(melt.me), ]
#   names(melt.me.keep) = c("ASV", "variable", "Name")
#   
#   #match sample names with glacier and region
#   meta = data.frame(sample_data(merge.ps))
#   res = melt.me.keep %>%
#     left_join(meta, by='Name')
#   
#   #How many modules are there and what percentage of the network do they account for?
#   #number of modules
#   length(unique(df$modularity_class))
#   
#   #get dataframe of modules
#   mod.df = data.frame(table(df$modularity_class))
#   mod.df$perc = round(mod.df$Freq/sum(mod.df$Freq)*100, 2)
#   mod.df = mod.df[order(-mod.df$Freq),]
#   
#   names(mod.df) = c("ModuleID", "Nodes", "Percentage")
#   
#   #number of modules with over 50 nodes
#   x = nrow(mod.df[mod.df$Nodes > 50,])
#   
#   # What is the taxonomic composition of each major module? 
#   
#   #get the numbers of our major modules
#   major = mod.df[mod.df$Nodes > 50,]
#   major.mods = as.numeric(major$ModuleID)
#   
#   #empty list for saving our results
#   #looking at major mods only 
#   abun.list = list()
#   phy.list = list()
#   order.list = list()
#   class.list = list()
#   glac.list = list()
#   reg.list=list()
#   key.list=list()
#   key.tab.list = list()
#   
#   for (i in 1:length(unique(df$modularity_class))){
#     
#     df.mod1 = df[df$modularity_class == unique(df$modularity_class)[i],]
#     
#     #abundance
#     abun.tab = as.data.frame(table(df.mod1$Abundance, useNA = c("ifany")))
#     abun.tab$perc = round(abun.tab$Freq/sum(abun.tab$Freq)*100, 2)
#     abun.tab = abun.tab[order(-abun.tab$Freq),]
#     abun.tab$module = mod.df$Module[i]
#     abun.list[[i]] = abun.tab
#     
#     #phylum
#     phy.tab = as.data.frame(table(df.mod1$Phylum, useNA = c("ifany")))
#     phy.tab$perc = round(phy.tab$Freq/sum(phy.tab$Freq)*100, 2)
#     phy.tab = phy.tab[order(-phy.tab$Freq),]
#     phy.tab$module = mod.df$Module[i]
#     phy.list[[i]] = phy.tab
#     
#     #class
#     class.tab = as.data.frame(table(df.mod1$Class, useNA = c("ifany")))
#     class.tab$perc = round(class.tab$Freq/sum(class.tab$Freq)*100, 2)
#     class.tab = class.tab[order(-class.tab$Freq),]
#     class.tab$module = mod.df$Module[i]
#     class.list[[i]] = class.tab
#     
#     #order
#     order.tab = as.data.frame(table(df.mod1$Order, useNA = c("ifany")))
#     order.tab$perc = round(order.tab$Freq/sum(order.tab$Freq)*100, 2)
#     order.tab = order.tab[order(-order.tab$Freq),]
#     order.tab$module = mod.df$Module[i]
#     order.list[[i]] = order.tab
#     
#     #And in which sites do those taxa occur? 
#     #reduce our res dataframe by those in this model
#     res.mod1 = res[res$ASV %in% df.mod1$name,]
#     unique(res.mod1$ASV) == unique(df.mod1$name) #should be TRUE
#     
#     #glacier
#     glac.tab = data.frame(table(res.mod1$Glacier, useNA = c("ifany")))
#     glac.tab$perc = round(glac.tab$Freq/sum(glac.tab$Freq)*100, 2)
#     glac.tab = glac.tab[order(-glac.tab$Freq),]
#     glac.tab$module = mod.df$Module[i]
#     glac.list[[i]] = glac.tab
#     
#     #region
#     reg.tab = data.frame(table(res.mod1$Region, useNA = c("ifany")))
#     reg.tab$perc = round(reg.tab$Freq/sum(reg.tab$Freq)*100, 2)
#     reg.tab = reg.tab[order(-reg.tab$Freq),]
#     reg.tab$module = mod.df$Module[i]
#     reg.list[[i]] = reg.tab
#     
#     #In which modules do the keystone taxa occur and what is there taxonomy?
#     #subset the modular data by keystone taxa
#     key.mod = df.mod1[df.mod1$name %in% key$X,]
#     #subset res by our keystone taxa so we can see which glaciers they're associated with
#     key.res = res[res$ASV %in% key.mod$name,]
#     key.tab = data.frame(table(key.res$Glacier))
#     key.list[[i]] = key.mod
#     key.tab.list[[i]] = key.tab
#     
#   }
#   
#   
#   abun.out = bind_rows(abun.list, .id = "Module")
#   phy.out = bind_rows(phy.list, .id = "Module")
#   class.out = bind_rows(class.list, .id = "Module")
#   order.out = bind_rows(order.list, .id = "Module")
#   glac.out = bind_rows(glac.list, .id = "Module")
#   reg.out = bind_rows(reg.list, .id = "Module")
#   key.out = bind_rows(key.list, .id = "Module")
#   key.tab.out = bind_rows(key.tab.list, .id = "Module")
#   
#   return(final.res.list = list(major, abun.out, phy.out, class.out, order.out, glac.out, reg.out, key.out, key.tab.out))
#   
# }
# 
# ant.mod.res = module_analysis(final.df.ant, pro.ant, euk.ant, keystone.ant)
# arc.mod.res = module_analysis(final.df.arc, pro.arc, euk.arc, keystone.arc)
# 
# ## Modular Analysis of the Antarctic and Arctic Networks
# 
# #Major modules of the Antarctic network, identified as those with over 50 nodes
# ant.mod.res[[1]]
# 
# #Major modules of the Arctic network, identified as those with over 50 nodes
# arc.mod.res[[1]]
# 
# #Subcommunity breakdown per Antarctic module
# ant.mod.res[[2]]
# 
# #Subcommunity breakdown per Arctic module
# arc.mod.res[[2]]
# 
# #Phylum composition per Antarctic module
# ant.mod.res[[3]]
# 
# #Phylum composition per Arctic module
# arc.mod.res[[3]]
# 
# #Class composition per Antarctic module
# ant.mod.res[[4]]
# 
# #Class composition per Arctic module
# arc.mod.res[[4]]
# 
# #Order composition per Antarctic module
# ant.mod.res[[5]]
# 
# #Order composition per Arctic module
# arc.mod.res[[5]]
# 
# #Glaciers associated with each Antarctic module
# ant.mod.res[[6]]
# 
# #Glaciers associated with each Arctic module
# arc.mod.res[[6]]
# 
# #Regions associated with each Antarctic module
# ant.mod.res[[7]]
# 
# #Regions associated with each Arctic module
# arc.mod.res[[7]]
# 
# key.tab.ant = cbind(ant.mod.res[[8]]$Module, ant.mod.res[[8]]$modularity_class, ant.mod.res[[8]]$Abundance, ant.mod.res[[8]]$Kingdom, ant.mod.res[[8]]$Phylum, ant.mod.res[[8]]$Class, ant.mod.res[[8]]$Order)
# 
# #Keystone taxa associated with each Antarctic module
# key.tab.ant
# 
# key.tab.arc = cbind(arc.mod.res[[8]]$Module, arc.mod.res[[8]]$modularity_class, arc.mod.res[[8]]$Abundance, arc.mod.res[[8]]$Kingdom, arc.mod.res[[8]]$Phylum, arc.mod.res[[8]]$Class, arc.mod.res[[8]]$Order)
# 
# #Keystone taxa associated with each Arctic module
# key.tab.arc
# 
# #The glaciers to which each keystone taxa are associated in each Antarctic module.
# ant.mod.res[[9]]
# 
# #The glaciers to which each keystone taxa are associated in each Arctic module.
# arc.mod.res[[9]]
# 
# #ternary diagrams
# 
# #at some point use ggtern to make tertiary diagrams
# 
# #ternary diagrams of modular distribution between abundant, intermediate and rare ASVs
# #ternary diagrams of modular distribution between greenland, svalbard and sweden for Arctic dataset
# 
# #Dataframe showing the regions each ASV within the network is associated to
# 
# #For each region in module I add the relative abundance of the ASVs 
# 
# # module_data = gephi.tab.arc
# # phylo1 = pro.arc
# # phylo2 = euk.arc
# # module = 9
# # name = "Module II"
# # 
# # #subset to module 1 (=4)
# # sort(table(module_data$modularity_class)) #module 4
# 
# ternary_plot_data_prep <- function(module_data, phylo1, phylo2, module, name){
#   
#   mod = module_data[module_data$modularity_class == module,] #these are the ASVs within this module
#   #subset the count table by these ASVs and agglomerate the samples into regions to test
#   ps = merge_my_phylo_same_pole(phylo1, phylo2)
#   mergedPS = merge_samples(ps, "Region")
#   counts = data.frame(t(otu_table(mergedPS)))
#   sum(counts) == sum(data.frame(otu_table(phylo1)))+sum(data.frame(otu_table(phylo2))) #should be TRUE
#   counts = counts[row.names(counts) %in% mod$name,]
#   
#   #transform to relative abundances
#   abundance <- counts
#   for (i in 1:ncol(abundance)) {
#     abundance[i] <- abundance[i] / sum(abundance[i])
#   }
#   abundance
#   abundance$Module = name
#   
#   return(abundance)
#   
# }
# 
# arc.mod1 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 1, "Module I")
# arc.mod2 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 2, "Module II")
# arc.mod3 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 3, "Module III")
# arc.mod4 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 4, "Module IV")
# arc.mod5 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 5, "Module V")
# arc.mod6 = ternary_plot_data_prep(gephi.tab.arc, pro.arc, euk.arc, 6, "Module VI")
# 
# arc.mod.final = rbind(arc.mod1, arc.mod2, arc.mod3, arc.mod4, arc.mod5, arc.mod6)
# 
# #plot
# ggtern(data=arc.mod.final,aes(x=Greenland,y=Svalbard, z=Sweden)) +
#   geom_point()+
#   theme_rgbw()+
#   stat_density_tern(geom = 'polygon',
#                     n         = 200,
#                     aes(fill  = ..level..,
#                         alpha = ..level..)) +
#   scale_fill_gradient(low = "blue",high = "red")  +
#   guides(color = "none", fill = "none", alpha = "none")+
#   facet_wrap(~ Module)
# 
# 
# #For the Antarctic divide up into Diamond, Ustein scoop and Other sites
# 
# # module_data = gephi.tab.ant
# # phylo1 = pro.ant
# # phylo2 = euk.ant
# # module = 4
# # name = "Module I"
# 
# # #subset to module 1 (=4)
# # sort(table(module_data$modularity_class)) #module 4
# 
# ternary_plot_data_prep2 <- function(module_data, phylo1, phylo2, module, name){
#   
#   mod = module_data[module_data$modularity_class == module,] #these are the ASVs within this module
#   #subset the count table by these ASVs and agglomerate the samples into regions to test
#   ps = merge_my_phylo_same_pole(phylo1, phylo2)
#   
#   meta = data.frame(sample_data(ps))
#   Region2 = meta$Glacier
#   out = sub("Miers", "Other", Region2)
#   out = sub("Lower Wright", "Other", out)
#   out = sub("Upper Wright", "Other", out)
#   out = sub("Canada", "Other", out)
#   out = sub("Commonwealth", "Other", out)
#   out = sub("Taylor", "Other", out)
#   out = sub("Upper Koettlitz", "Other", out)
#   out = sub("Lower Koettlitz", "Other", out)
#   out = sub("Usteinen Scoop", "Queen Maud Land", out)
#   out = sub("Duboisbreen", "Queen Maud Land", out)
#   meta$Region2 = out
#   merge.ps = phyloseq(otu_table(ps), tax_table(ps),sample_data(meta))
#   mergedPS = merge_samples(merge.ps, "Region2")
#   counts = data.frame(t(otu_table(mergedPS)))
#   sum(counts) == sum(data.frame(otu_table(phylo1)))+sum(data.frame(otu_table(phylo2))) #should be TRUE
#   counts = counts[row.names(counts) %in% mod$name,]
#   
#   #transform to relative abundances
#   abundance <- counts
#   for (i in 1:ncol(abundance)) {
#     abundance[i] <- abundance[i] / sum(abundance[i])
#   }
#   abundance
#   abundance$Module = name
#   
#   return(abundance)
#   
# }
# 
# ant.mod1 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 1, "Module I")
# ant.mod2 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 2, "Module II")
# ant.mod3 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 3, "Module III")
# ant.mod4 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 4, "Module IV")
# ant.mod5 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 5, "Module V")
# ant.mod6 = ternary_plot_data_prep2(gephi.tab.ant, pro.ant, euk.ant, 6, "Module VI")
# 
# ant.mod.final = rbind(ant.mod1, ant.mod2, ant.mod3, ant.mod4, ant.mod5, ant.mod6)
# 
# 
# #plot
# ggtern(data=ant.mod.final,aes(x=Diamond,y=Other, z=Queen.Maud.Land)) +
#   geom_point()+
#   theme_rgbw()+
#   stat_density_tern(geom = 'polygon',
#                     n         = 200,
#                     aes(fill  = ..level..,
#                         alpha = ..level..)) +
#   scale_fill_gradient(low = "blue",high = "red")  +
#   guides(color = "none", fill = "none", alpha = "none")+
#   facet_wrap(~ Module)
# 
# #By abundance
# 
# # I want a dataframe with ASVs as rows, abundance classifications as columns and the ASVs mean relative abundance in the dataset as data
# ternary_plot_data_prep3 <- function(module_data, phylo1, phylo2, module, name, df1, df2){
#   
#   mod = module_data[module_data$modularity_class == module,] #these are the ASVs within this module
#   
#   #get count data
#   ps = merge_my_phylo_same_pole(phylo1, phylo2)
#   #transform into relative abundance 
#   counts1 = data.frame(otu_table(transform_sample_counts(ps, function(x) x / sum(x) )))
#   #get row means for each ASV
#   abun.res = data.frame(rowMeans(counts1))
#   #classify as abundant, intermediate or rare 
#   df = rbind(df1, df2)
#   abun.res$ID = rownames(abun.res)
#   #match relative abundance to ASV 
#   out = left_join(abun.res, df, by="ID")
#   #remove ASVs that aren't in this module
#   out = out[out$ID %in% mod$name,]
#   #divide dataframe by abundance
#   x = reshape(out, idvar = "ID", timevar = "Abundance", direction = "wide")
#   rownames(x) <- x[,1]
#   x[,1] <- NULL
#   #remove colour columns
#   x <- x %>% select(-contains("Colour"))
#   names = str_sub(names(x),-2,-1)
#   names(x) = names
#   x[is.na(x)] <- 0
#   x$Module = name
#   
#   return(x)
#   
# }
# 
# arc.mod1a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 4, "Module I", pro.arc.df, euk.arc.df)
# arc.mod2a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 9, "Module II", pro.arc.df, euk.arc.df)
# arc.mod3a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 11, "Module III", pro.arc.df, euk.arc.df)
# arc.mod4a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 1, "Module IV", pro.arc.df, euk.arc.df)
# arc.mod5a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 3, "Module V", pro.arc.df, euk.arc.df)
# arc.mod6a = ternary_plot_data_prep3(gephi.tab.arc, pro.arc, euk.arc, 13, "Module VI", pro.arc.df, euk.arc.df)
# 
# arc.mod.final2 = dplyr::bind_rows(arc.mod1a, arc.mod2a, arc.mod3a, arc.mod4a, arc.mod5a, arc.mod6a)
# 
# arc.mod.final2[is.na(arc.mod.final2)] <- 0
# names(arc.mod.final2) <- c("Intermediate", "Abundant", "Module", "Rare")
# 
# arc.mod.final2$Intermediate = arc.mod.final2$Intermediate+0.002
# arc.mod.final2$Abundant = arc.mod.final2$Abundant+0.001
# arc.mod.final2$Rare = arc.mod.final2$Rare+0.003
# 
# #plot
# ggtern(data=arc.mod.final2, aes(x=Abundant,y=Intermediate, z=Rare)) +
#   geom_point()+
#   theme_rgbw()+
#   stat_density_tern(geom = 'polygon',
#                     n         = 200,
#                     aes(fill  = ..level..,
#                         alpha = ..level..)) +
#   scale_fill_gradient(low = "blue",high = "red")  +
#   guides(color = "none", fill = "none", alpha = "none")+
#   facet_wrap(~ Module)