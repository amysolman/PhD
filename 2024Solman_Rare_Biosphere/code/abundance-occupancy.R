# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

# library(ggpubr)
#install.packages("ggstatsplot")
# library(ggstatsplot)
source("00-solman-functions.R")
# library(ggcorrplot)
#library(fdrtool)
# library(psych)
# library(cowplot)
# #devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)
library(stringr) #for subsetting strings in lat long function
library(fossil) #for earth.dist function
library(vegan) #for community dissimilarity matrix
library(tibble) #for rownames to column
library(dplyr) #coalesce function + rows update function
library(reshape2) #for melt function
library(funrar)
library(phyloseq)
library(ggplot2)
#install.packages("ggpmisc")
library(ggpmisc)
library(picante) #for faith's pd
library(cowplot)

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds")
euk <- readRDS("../results/18S-phylo-object-rarefied.rds")

# pro <- readRDS("../results/16S-phylo-object-sub-coms-merged.rds") 
# euk <- readRDS("../results/18S-phylo-object-sub-coms-merged.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
pro.arc <- subset_samples(pro, Pole=="Arctic")
pro.arc = filter_taxa(pro.arc, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)
euk.arc <- subset_samples(euk, Pole=="Arctic")
euk.arc = filter_taxa(euk.arc, function(x) sum(x) >= 1, TRUE)

#Arctic prokaryotes
pro.arc.abun <- readRDS("../results/16S-phylo-object-arc-abun.rds") 
pro.arc.int <- readRDS("../results/16S-phylo-object-arc-int.rds") 
pro.arc.rare <- readRDS("../results/16S-phylo-object-arc-rare.rds") 

#Antarctic prokaryotes
pro.ant.abun <- readRDS("../results/16S-phylo-object-ant-abun.rds") 
pro.ant.int <- readRDS("../results/16S-phylo-object-ant-int.rds") 
pro.ant.rare <- readRDS("../results/16S-phylo-object-ant-rare.rds") 

#Arctic eukaryotes
euk.arc.abun <- readRDS("../results/18S-phylo-object-arc-abun.rds") 
euk.arc.int <- readRDS("../results/18S-phylo-object-arc-int.rds") 
euk.arc.rare <- readRDS("../results/18S-phylo-object-arc-rare.rds") 

#Antarctic eukaryotes
euk.ant.abun <- readRDS("../results/18S-phylo-object-ant-abun.rds") 
euk.ant.int <- readRDS("../results/18S-phylo-object-ant-int.rds") 
euk.ant.rare <- readRDS("../results/18S-phylo-object-ant-rare.rds") 


#CREATE ABUNDANCE-OCCUPANCY PLOT FUNCTION
# my_phylo = pro.ant
# pole = "Antarctic"
# # title = "Blah"
# abun = readRDS("../results/16S-phylo-object-ant-abun.rds") 
# int =  readRDS("../results/16S-phylo-object-ant-int.rds") 
# rare =  readRDS("../results/16S-phylo-object-ant-rare.rds") 

abundance_occupancy_plot <- function(my_phylo, abun, int, rare){
  
  ASV.table = data.frame(t(otu_table(my_phylo)), check.names = FALSE)
  #Get the mean number of reads per sample (mean sum of each row) e.g. individuals per community 
  N <- mean(apply(ASV.table, 1, sum))
  #Get the mean number of reads for each ASV across all samples (mean of each column)
  p.m <- apply(ASV.table, 2, mean)
  #Remove any zeros
  p.m <- p.m[p.m != 0] #remove any ASVs with zero counts
  #divide the number of mean reads by the total number of reads per sample (mean relative abundance of each ASV globally)
  p <- p.m/N #mean relative abundance of each ASV. 
  #Make ASV.table into presence/absence table
  ASV.table.bi <- 1*(ASV.table>0)
  #find the mean of frequence of each column (ASV) e.g. the mean number of samples each ASV is found in
  freq.table <- apply(ASV.table.bi, 2, sum) 
  #only keep ASVs with a frequency other than 0
  freq.table <- freq.table[freq.table != 0] #only keep data that isn't zero
  #Put the average relative abundance of each taxa into a dataframe
  p.df = data.frame(p) %>%
    rownames_to_column(var="ASV") 
  #Make into dataframe with ASV name and frequence of occurence (percentage of samples the ASV is present in )
  freq.df = data.frame(ASV=names(freq.table), freq=freq.table) 
  #Combine dataframes and arrange by relative abundance 
  C <- inner_join(p.df,freq.df, by="ASV") %>%
    arrange(p)
  # Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
  C.no0 <- C %>%
    filter(freq != 0, p != 0)
  
  #get subcommunity classifications
  asv_subs = data.frame(ASV = c(rownames(data.frame(tax_table(abun))),
                               rownames(data.frame(tax_table(int))),
                               rownames(data.frame(tax_table(rare)))), 
                        Subcommunity=c(rep("Abundant", length(rownames(data.frame(tax_table(abun))))),
                                       rep("Intermediate", length(rownames(data.frame(tax_table(int))))),
                                       rep("Rare", length(rownames(data.frame(tax_table(rare)))))))
  
  #combine with abundant/occupancy data
  # tax = rbind(abun_tax_names, int_tax_names, rare_tax_names)
  # tax = as.data.frame(tax)
  # colnames(tax) = c("ID", "Subcommunity")
  #df$ID = rownames(df)
  df2 = full_join(C.no0, asv_subs, by="ASV")
  colnames(df2) = c("ASV", "rel_abundance", "occupancy", "Subcommunity")
  
  #Calculate spearman's rank correlation
  #Abundant community
  test_df1 = df2[df2$Subcommunity == "Abundant",]
  res1 <- cor.test(test_df1$occupancy, test_df1$rel_abundance, method="spearman", exact=FALSE) #exact=FALSE so we don't get the "cannot compute exact p-value with ties" error
  test_df2 = df2[df2$Subcommunity == "Intermediate",]
  res2 <- cor.test(test_df2$occupancy, test_df2$rel_abundance, method="spearman", exact=FALSE)
  test_df3 = df2[df2$Subcommunity == "Rare",]
  res3 <- cor.test(test_df3$occupancy, test_df3$rel_abundance, method="spearman", exact=FALSE)
  #look at our results
  res1
  res2
  res3
  
  #pull out the rho value
  rho_val1 <- as.numeric(res1$estimate)
  rho_val2 <- as.numeric(res2$estimate)
  rho_val3 <- as.numeric(res3$estimate)
  
  #pull out p-value
  p_val1 <- mod_my_p_val(res1$p.value)
  p_val2 <- mod_my_p_val(res2$p.value)
  p_val3 <- mod_my_p_val(res3$p.value)
  
  df2$Subcommunity <- factor(df2$Subcommunity, levels=c("Rare", "Intermediate", "Abundant"))
  
  #add rho value, p value and N = number of asvs to the plot
  my_plot <- ggplot(df2, aes(y=occupancy, x=log10(rel_abundance), color=Subcommunity)) +
    geom_point()+
    scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    #annotate("text", x=4, y=-5, label=paste0("rho=", round(rho_val, 2)))+
    xlim(-6, -0.5)+
    ylim(0, 55)+
    annotate("text", y=49, x=-5, label=paste0("N=", nrow(test_df1)), size=3, color="#ef7a76")+
    annotate("text", y=49, x=-4.3, label=paste0("rho=", round(rho_val1, 2)), size=3, color="#ef7a76")+
    annotate("text", y=49, x=-3.6, label= p_val1, size=3, color="#ef7a76")+
    annotate("text", y=52, x=-5, label=paste0("N=", nrow(test_df2)), size=3, color="#74a9d8")+
    annotate("text", y=52, x=-4.3, label=paste0("rho=", round(rho_val2, 2)), size=3, color="#74a9d8")+
    annotate("text", y=52, x=-3.6, label= p_val2, size=3, color="#74a9d8")+
    annotate("text", y=55, x=-5, label=paste0("N=", nrow(test_df3)), size=3, color="#82ca81")+
    annotate("text", y=55, x=-4.3, label=paste0("rho=", round(rho_val3, 2)), size=3, color="#82ca81")+
    annotate("text", y=55, x=-3.6, label= p_val3, size=3, color="#82ca81")+
    ylab("Sites Occupied")+
    xlab("log10(Mean relative abundance)")+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    #ggtitle(title)+
    theme_bw()+
    theme(text = element_text(size = 15), legend.text=element_text(size=15))
  
  my_plot
  
  return(my_plot)
  
}


out1 = abundance_occupancy_plot(pro.arc, pro.arc.abun, pro.arc.int, pro.arc.rare)
out2 = abundance_occupancy_plot(pro.ant, pro.ant.abun, pro.ant.int, pro.ant.rare)
out3 = abundance_occupancy_plot(euk.arc, euk.arc.abun, euk.arc.int, euk.arc.rare)
out4 = abundance_occupancy_plot(euk.ant, euk.ant.abun, euk.ant.int, euk.ant.rare)

#combine plots
multi.p = plot_grid(
  out2 + theme(legend.position = "none"),
  out1 + theme(legend.position = "none"),
  out4 + theme(legend.position = "none"),
  out3 + theme(legend.position = "none"),
  labels = c('A', 'B', 'C', 'D'),
  ncol = 2)
#rel_widths = c(1, .6))

legend <- get_legend(
  out1 + 
    #guides(color = guide_legend(nrow = 1, override.aes = list(size = 10))) +
    theme(legend.position = "bottom",
          axis.text=element_text(size=25),
          axis.title=element_text(size=25,face="bold"),
          legend.text = element_text(size=15), 
          legend.title = element_blank())
)

final.p = plot_grid(multi.p, legend, ncol = 1, rel_heights = c(1, .1))
final.p

pdf("../results/abundance-occupancy-plot.pdf", width=10, height=8)
print(final.p)
dev.off()