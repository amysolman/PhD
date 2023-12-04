# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Define abundance cut offs to be tested
# 4. Get truncated datasets (removing the rarest ASVs and keeping abundant ASVs) and calculate differences to plot
# 5. Plot the results to find appropriate abundant community cut off
# 6. Get truncated datasets (adding the rarest ASVs) and calculate differences to plot
# 7. Plot the results to find appropriate rare community cut off
# 8. Combine plots together

# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(vegan)
library(MASS)
library(ggplot2)
#library(tidyverse)
#install.packages("funrar")
library(funrar)
library(phyloseq)
#library(DESeq2)
library(cowplot)

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds") 

# 3. Define abundance cut offs to be tested and percentages to remove/retain

#Abundant thresholds to plot
abun_cut_1 = 0.001 #0.1%
abun_cut_2 = 0.0005 #0.05%
abun_cut_3 = 0.0001 #0.01%

#Rare thresholds to plot
rare_cut_1 = 0.0005 #0.05%
rare_cut_2 = 0.0001 #0.01%
rare_cut_3 = 0.00005 #0.005%

# 4. Get truncated datasets (removing the rarest ASVs and keeping abundant ASVs) and calculate differences to plot

#Function for truncating data by removing rare taxa
remove_rare_plot_df <- function(phylo, dataset){
  
  #percentages to remove
  perc = c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 98)
  
  #get count table
  df = data.frame(t(otu_table(phylo)))
  
  #Sort data frame according the abundance
  order_df = df[,order(colSums(-df))]
  
  #how many ASVs/columns in 1%
  perc_1 = ncol(order_df)/100
  
  #list for saving truncated data sets
  trunc_data <- list()
  
  for (i in 1:length(perc)){
    #remove perc[i]*perc_1 number of columns of rarest types and store the truncated data frame
    num_rem = round(perc[i]*perc_1)
    trunc_data[[i]] <- data.frame(order_df[1:(length(order_df)-num_rem)])
  }
  
  cor_res <- vector()
  pro_res <- vector()
  
  # x = trunc_data[[22]]
  # colSums(x)
  
  for (j in 1:length(trunc_data)){
    
    #Calculate pairwise distance matrices between original and truncated matrices using the Bray–Curtis dissimilarity index
    
    #distance matrix of the original dataset
    orig_dist <-vegan::vegdist(order_df,method="bray")	
    
    #distance matrix of truncated dataset
    trunc_dist <-vegan::vegdist(trunc_data[[j]],distance="bray")
    
    #Spearman correlation between matrices
    mat_cor <-cor.test(orig_dist,trunc_dist,method="spearman") #correlation between matrices
    
    #replace 0 with 10e-20 for NMDS
    orig_dist2<-orig_dist
    orig_dist2[orig_dist2==0]<-10e-20	
    trunc_dist2<-trunc_dist
    trunc_dist2[trunc_dist2==0]<-10e-20
    
    #Caluclate NMDS
    orig_dist2NMDS<-MASS::isoMDS(orig_dist2,trace=0)	#NMDS for original dataset
    trunc_dist2NMDS<-MASS::isoMDS(trunc_dist2,trace=0)		#NMDS for truncated dataset
    
    #Get procrustes value from NMDS
    pro_val<-vegan::protest(orig_dist2NMDS,trunc_dist2NMDS)
    cor_res[j] = mat_cor$estimate
    pro_res[j] = pro_val$t0
    
  }
  
  #Calculate the percentage of ASVs above these thresholds
  #Find relative abundances of all ASVs
  #total number of counts 
  tot_counts = sum(colSums(order_df))
  #relative abundance of each ASV
  rel_abun = colSums(order_df)/tot_counts
  
  sum(rel_abun) #should equal 1
  
  #how many ASVs are ABOVE our thresholds
  num1 = rel_abun[rel_abun > abun_cut_1]
  num2 = rel_abun[rel_abun > abun_cut_2]
  num3 = rel_abun[rel_abun > abun_cut_3]
  
  #What are these values as percentages to plot
  #number of ASVs above the abundance threshold e.g. 187 / number of ASVs in 1% of the data set (64) = the percentage of ASVs   the number of ASVs ABOVE the abundance threshold represents = 3
  perc_to_plot1 = round(length(num1)/perc_1) #about X% of ASVs are above this threshold
  perc_to_plot2 = round(length(num2)/perc_1) #about Y% of ASVs are above this threshold
  perc_to_plot3 = round(length(num3)/perc_1) #about Z% of ASVs are above this threshold
  
  #Save the results of our analysis
  value = c(cor_res, pro_res)
  Statistic = c(rep("Spearman's Correlation", length(cor_res)), rep("Procruste's Value", length(pro_res)))
  
  plot_df <- data.frame(perc, value, Statistic)
  write.csv(plot_df, paste0("../results/", dataset, "-remove-rare-df.csv"))
  
  
  res.list = list(plot_df, perc_to_plot1, perc_to_plot2, perc_to_plot3)
  
  return(res.list)

  
}

pro_res = remove_rare_plot_df(pro, "16S")
euk_res = remove_rare_plot_df(euk, "18S")

# 5. Plot the results to find appropriate abundant community cut off


MultiCoLA_plot_abun <- function(multi_res, dataset){
  
  plot_df = multi_res[[1]]
  
  #find the first time the value drops below 0.99
  # perc_to_plot = plot_df[plot_df$value < 0.99,]
  # perc_to_plot = min(perc_to_plot$perc)
  
  #find where both values drop below 0.99
  df1 = plot_df[plot_df$Statistic == "Spearman's Correlation",]
  df2 = plot_df[plot_df$Statistic == "Procruste's Value",]
  
  perc_to_plot1 = df1[df1$value < 0.99,]
  perc_to_plot2 = df2[df2$value < 0.99,]
  
  #which values do they share?
  percs = intersect(perc_to_plot1$perc, perc_to_plot2$perc)
  
  # #if they don't both go below 0.99 then set the cut off to 98%
  # if (length(percs) > 0){
  perc_to_plot = min(percs,na.rm=TRUE)
  # } else if (length(percs) == 0){
  #   perc_to_plot = 98
  # }
  
  base_plot <- ggplot(plot_df,aes(x=perc, y=value, group=Statistic))+
    geom_point()+
    geom_line(aes(linetype = Statistic), size=1)+
    xlab("% rare ASVs removed")+
    ylim(0, 1)+
    ylab("Coefficient")+
    #geom_vline(xintercept = perc_to_plot, linetype="dotted", 
    #color = "orange", size=2)+
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    theme(legend.text = element_text(size=5), 
          legend.title = element_blank())
  #ggtitle("Prokaryote Communities")
  base_plot
  
  
  verticals <- data.frame(
    intercepts=c(100-multi_res[[2]], 100-multi_res[[3]], 100-multi_res[[4]], perc_to_plot),
    Cutoffs=c("0.1%", "0.05%", "0.01%", "Coefficient Below 0.99"), color=c("#E69F00", "#56B4E9", "#009E73", "#999999"), Type=c("solid", "solid", "solid", "solid")
  )
  
  # per_dif <- data.frame(
  #   intercepts=c(perc_to_plot),
  #   'Below 0.9'=c(" ")
  # )
  
  second_plot <- base_plot +
    geom_vline(
      data=verticals,
      mapping=aes(xintercept=intercepts, color=Cutoffs),
      linetype=verticals$Type,
      size=1,
      key_glyph="path"   # this makes the legend key horizontal lines, not vertical
    )+
    scale_colour_manual(values=verticals$color)
  
  second_plot
  
  ggsave(paste0("../results/", dataset, "-dataset-based-remove-rare-plot.pdf"), width=10, height=7)
  
  return(second_plot)
  
}


pro_plot1 = MultiCoLA_plot_abun(pro_res, "16S")
euk_plot1 = MultiCoLA_plot_abun(euk_res, "18S")

# 6. Get truncated datasets (adding the rarest ASVs) and calculate differences to plot

#Function for truncating data by removing rare taxa
keep_rare_plot_df <- function(phylo, dataset){
  
  #percentages to remove
  perc = c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 98)
  
  #get count table
  df = data.frame(t(otu_table(phylo)))
  
  #Sort data frame according the abundance
  order_df = df[,order(colSums(-df))]
  
  #how many ASVs/columns in 1%
  perc_1 = ncol(order_df)/100
  
  trunc_data_abun <- list()
  
  for (i in 1:length(perc)){
    num_keep = round(perc[i]*perc_1)
    if (num_keep != 0){
      trunc_data_abun[[i]] <- data.frame(order_df[,(ncol(order_df)-num_keep+1):ncol(order_df)])
    } else {
      trunc_data_abun[[i]] <- data.frame(rep(NA, length(order_df[,1])))
    }
    if (perc[i] == 100){
      trunc_data_abun[[i]] <- data.frame(order_df[1:(length(order_df))])
    }
    
  }
  
  cor_res_abun <- vector()
  pro_res_abun <- vector()
  for (j in 1:(length(trunc_data_abun)-1)){
    
    #remove samples from truncated dataset with zero counts
    trunc_df = trunc_data_abun[[j+1]]
    trunc_df = trunc_df[rowSums(trunc_df) > 0, ]
    
    #only keep samples from original dataframe that match the truncated dataset
    original_df = order_df[row.names(order_df) %in% row.names(trunc_df),]
    
    #dissimilarity matrix of original dataset
    orig_dist <-vegdist(original_df,method="bray", na.rm=TRUE)
    
    #Calculate pairwise distance matrices between original and truncated matrices using the Bray–Curtis dissimilarity index
    #distance matrix of truncated dataset
    trunc_dist_abun <-vegdist(trunc_df,method="bray", na.rm=TRUE)
    #Spearman correlation between matrices
    mat_cor_abun <-cor.test(orig_dist,trunc_dist_abun,method="spearman") #correlation between matrices
    
    #replace 0 with 10e-20 for NMDS
    orig_dist2<-orig_dist
    orig_dist2[orig_dist2==0]<-10e-20
    trunc_dist_abun2<-trunc_dist_abun
    trunc_dist_abun2[trunc_dist_abun2==0]<-10e-20
    
    #Caluclate NMDS
    orig_dist_NMDS<-isoMDS(orig_dist2,trace=0)		#NMDS for original dataset
    trunc_dist_abun2NMDS<-isoMDS(trunc_dist_abun2,trace=0)		#NMDS for truncated dataset
    
    #Get procrustes value from NMDS
    pro_val_abun<-protest(orig_dist_NMDS,trunc_dist_abun2NMDS)
    cor_res_abun[j+1] = mat_cor_abun$estimate
    pro_res_abun[j+1] = pro_val_abun$t0
    
  }
  
  #Calculate the percentage of ASVs above these thresholds
  #Find relative abundances of all ASVs
  #total number of counts 
  tot_counts = sum(colSums(order_df))
  #relative abundance of each ASV
  rel_abun = colSums(order_df)/tot_counts
  
  sum(rel_abun) #should equal 1
  
  #how many ASVs are ABOVE our thresholds
  rare_num1 = rel_abun[rel_abun < rare_cut_1]
  rare_num2 = rel_abun[rel_abun < rare_cut_2]
  rare_num3 = rel_abun[rel_abun < rare_cut_3]
  
  #What are these values as percentages to plot
  rare_perc_to_plot1 = round(length(rare_num1)/perc_1)
  rare_perc_to_plot2 = round(length(rare_num2)/perc_1)
  rare_perc_to_plot3 = round(length(rare_num3)/perc_1)
  
  #Save the results of our analysis
  value = c(cor_res_abun, pro_res_abun)
  Statistic = c(rep("Spearman's Correlation", length(cor_res_abun)), rep("Procruste's Value", length(pro_res_abun)))
  
  rare_plot_df <- data.frame(perc, value, Statistic)
  write.csv(rare_plot_df, paste0("../results/", dataset, "-dataset-based-add-rare-df.csv"))
  
  
  res.list = list(rare_plot_df, rare_perc_to_plot1, rare_perc_to_plot2, rare_perc_to_plot3)
  
  return(res.list)
  
}

pro_res_rare = keep_rare_plot_df(pro, "16S")
euk_res_rare = keep_rare_plot_df(euk, "18S")

# 7. Plot the results to find appropriate rare community cut off


MultiCoLA_plot_rare <- function(multi_res, dataset){
  
  plot_df = multi_res[[1]]
  
  # rare_perc_to_plot = plot_df[plot_df$value > 0.8,]
  # rare_perc_to_plot = rare_perc_to_plot[complete.cases(rare_perc_to_plot), ]
  # rare_perc_to_plot = min(rare_perc_to_plot$perc)
  
  #find where both values goes above 0.8
  df1 = plot_df[plot_df$Statistic == "Spearman's Correlation",]
  df2 = plot_df[plot_df$Statistic == "Procruste's Value",]
  
  perc_to_plot1 = df1[df1$value > 0.8,]
  perc_to_plot2 = df2[df2$value > 0.8,]
  
  #which values do they share?
  percs = intersect(perc_to_plot1$perc, perc_to_plot2$perc)
  rare_perc_to_plot = min(percs,na.rm=TRUE)
  
  base_plot <- ggplot(plot_df,aes(x=perc, y=value, group=Statistic))+
    geom_point()+
    geom_line(aes(linetype = Statistic), size=1)+
    xlab("% rare ASVs retained")+
    ylim(0, 1)+
    ylab("Coefficient")+
    #geom_vline(xintercept = perc_to_plot, linetype="dotted", 
    #color = "orange", size=2)+
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10,face="bold"))+
    theme(legend.text = element_text(size=5), 
          legend.title = element_blank())
  #ggtitle("Prokaryote Communities")
  base_plot
  
  
  verticals <- data.frame(
    intercepts=c(multi_res[[2]], multi_res[[3]], multi_res[[4]], rare_perc_to_plot),
    Cutoffs=c("0.05%", "0.01%", "0.005%", "Coefficient Above 0.8"), color=c("#E69F00", "#56B4E9", "#009E73", "#999999"), Type=c("solid", "solid", "solid", "solid")
  )
  
  # per_dif <- data.frame(
  #   intercepts=c(perc_to_plot),
  #   'Below 0.9'=c(" ")
  # )
  
  second_plot <- base_plot +
    geom_vline(
      data=verticals,
      mapping=aes(xintercept=intercepts, color=Cutoffs),
      linetype=verticals$Type,
      size=1,
      key_glyph="path"   # this makes the legend key horizontal lines, not vertical
    )+
    scale_colour_manual(values=verticals$color)
  
  second_plot
  
  ggsave(paste0("../results/", dataset, "-dataset-based-keep-rare-plot.pdf"), width=10, height=7)
  
  return(second_plot)
  
}

pro_plot1_rare = MultiCoLA_plot_rare(pro_res_rare, "16S")

euk_plot1_rare = MultiCoLA_plot_rare(euk_res_rare, "18S")

#8. Combine plots together

prow1 <- plot_grid(
  pro_plot1 + theme(legend.position="none"),
  euk_plot1 + theme(legend.position="none"),
  #align = 'vh',
  labels = c("A", "B")
)

legend_b1 <- get_legend(
  pro_plot1 + 
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3),title=NULL)) +
    theme(legend.position = "bottom",
          axis.text=element_text(size=5),
          axis.title=element_text(size=5,face="bold"),
          legend.text = element_text(size=10), 
          legend.title = element_blank())
)

p1 = plot_grid(prow1, legend_b1, ncol = 1, rel_heights = c(1, .1))


prow2 <- plot_grid(
  pro_plot1_rare + theme(legend.position="none"),
  euk_plot1_rare + theme(legend.position="none"),
  #align = 'vh',
  labels = c("C", "D")
)


legend_b2 <- get_legend(
  pro_plot1_rare + 
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3),title=NULL)) +
    theme(legend.position = "bottom",
          axis.text=element_text(size=5),
          axis.title=element_text(size=5,face="bold"),
          legend.text = element_text(size=10), 
          legend.title = element_blank())
)


p2 = plot_grid(prow2, legend_b2, ncol = 1, rel_heights = c(1, .1))

p3 = plot_grid(p1, p2, ncol=1)

pdf("../results/multicola-plot.pdf", width = 10, height = 5)
print(p3)
dev.off()
