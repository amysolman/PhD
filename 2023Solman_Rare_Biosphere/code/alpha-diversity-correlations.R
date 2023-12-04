# **Script breakdown**
# 
# 1. Clear workspace and load packages
# 2. Import data
# 3. Define metadata editing function
# 4. Calculate alpha diversity for analysis, prepare metadata for analysis and plot correlations


## Method

#Spearman’s rank correlation was used to test the impact of environmental factors on alpha diversity within the Antarctic samples, with p-values adjusted using the false discovery rate method and considered significant at p < 0.05. Correlations were calculated using the corr.test function in the psych package and visualised using the ggstatsplot and ggcorrplot packages. 

# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(ggpubr)
library(phyloseq)
library(dplyr) #coalesce function
#install.packages("ggstatsplot")
library(ggstatsplot)
source("00-solman-functions.R")
library(ggcorrplot)
#library(fdrtool)
library(psych)
library(cowplot)
library(naniar) #replacing values with NA

# #devtools::install_github("caijun/ggcorrplot2")
# library(ggcorrplot2)

## Results 

# 2. Import data
pro <- readRDS("../results/16S-phylo-object-rarefied.rds") 
euk <- readRDS("../results/18S-phylo-object-rarefied.rds")

pro.ant <- subset_samples(pro, Pole=="Antarctic")
pro.ant = filter_taxa(pro.ant, function(x) sum(x) >= 1, TRUE)
euk.ant <- subset_samples(euk, Pole=="Antarctic")
euk.ant = filter_taxa(euk.ant, function(x) sum(x) >= 1, TRUE)

pro.ant.abun <- readRDS("../results/16S-phylo-object-rarefied-ant-abun.rds")
pro.ant.int <- readRDS("../results/16S-phylo-object-rarefied-ant-int.rds")
pro.ant.rare <- readRDS("../results/16S-phylo-object-rarefied-ant-rare.rds")

euk.ant.abun <- readRDS("../results/18S-phylo-object-rarefied-ant-abun.rds")
euk.ant.int <- readRDS("../results/18S-phylo-object-rarefied-ant-int.rds")
euk.ant.rare <- readRDS("../results/18S-phylo-object-rarefied-ant-rare.rds")

#PROBLEM WITH EUKARYOTE DATA
#SEEMS LIKE VERY FEW ASVS ARE OBSERVED IN THE SAME SAMPLES

#lets look at the antarctic count table for eukaryotes
counts = data.frame(otu_table(euk.ant))
#turn into presence/absence dataframe
count.pa = counts
count.pa[count.pa > 1] = 1
#now get a histogram of the row sums
hist(rowSums(count.pa))
#we can see that almost all ASVs are found in only one sample, why? This is different to how my data was before.

#correlation plot with ggcorplot + adjust p values + report rho, CIs and p-vale

# phylo1 = euk.ant.abun
# phylo2 = euk.ant.int
# phylo3 = euk.ant.rare

get_my_cor_plots <- function(phylo1, phylo2, phylo3){
  
  
  #put phyloseq objects into a list
  phy_list <- list(phylo1, phylo2, phylo3)
  
  #list for saving out correlations
  out.corr = list()
  
  for (i in 1:length(phy_list)){
    
    #get alpha diversity
    df = estimate_richness(phy_list[[i]], measures=c("Observed", "Shannon"))
    rownames(df) = rownames(data.frame(sample_data(phy_list[[i]])))
    
    #get meta data
    meta.out = data.frame(sample_data(phy_list[[i]]))
    meta.out = meta.out[,c(9:ncol(meta.out))]
    #meta.out = alpha_meta_df2(phy_list[[i]])
    
    #metadata to test
    # keep = c("Latitude", "DistanceToSea", "Altitude", "WaterDepth", "TotalDepth", "SedimentDepth", "Area", "Ice.lid","Temp","EC","pH","pCO2", "DOC","DIC","C.N","Cl.age","DO","N","C","H","HCO3",  "NH4",            "NO2","NO3","TN","TIN","TON","DON","TDN","TDP","DOP","PO4","SiO2","Cl","SO4","Na","K","Mg","Ca","F","DRP")
    # meta.out = meta.out[,names(meta.out) %in% keep]
    
    #correlation test that gives rho and adjusted p-values
    corr.res = corr.test(meta.out, df, use = "pairwise",method="spearman",adjust="fdr",
                         alpha=.05,ci=TRUE,minlength=5,normal=TRUE)
    
    #plot
    #ggcorrplot(corr.res[[1]], method = "circle", lab=TRUE, sig.level = 0.05, p.mat = corr.res[[5]])
    
    #add the correlation to a list
    out.corr[[i]] = corr.res
    
  }
  
  #make results into joined dataframes
  #rho values
  rho1 = as.data.frame(out.corr[[1]][[1]])
  rho2 = as.data.frame(out.corr[[2]][[1]])
  rho3 = as.data.frame(out.corr[[3]][[1]])
  
  #adjusted p-values
  pval1 = as.data.frame(out.corr[[1]][[5]])
  pval2 = as.data.frame(out.corr[[2]][[5]])
  pval3 = as.data.frame(out.corr[[3]][[5]])
  
  df = as.data.frame(cbind(rho1$Observed, rho2$Observed, rho3$Observed, rho1$Shannon, rho2$Shannon, rho3$Shannon))
  names(df) = c("Abundant (O)", "Intermediate (O)", "Rare (O)", "Abundant (S)", "Intermediate (S)", "Rare (S)")
  rownames(df) = rownames(rho1)
  df = as.matrix(df)
  
  pval.df = as.data.frame(cbind(pval1$Observed, pval2$Observed, pval3$Observed, pval1$Shannon, pval2$Shannon, pval3$Shannon))
  names(pval.df) = c("Abundant (O)", "Intermediate (O)", "Rare (O)", "Abundant (S)", "Intermediate (S)", "Rare (S)")
  rownames(pval.df) = rownames(pval1)
  pval.df = as.matrix(pval.df)
  
  
  p = ggcorrplot::ggcorrplot(df, method = "circle", lab=TRUE, sig.level = 0.05, p.mat = pval.df, ggtheme = ggplot2::theme_bw, show.legend = TRUE,
                             insig="blank", lab_size = 4)
  #+theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  #print(p)
  
  res.list= list(p, out.corr)
  return(res.list)
  
}


pro.cor.plot <- get_my_cor_plots(pro.ant.abun, pro.ant.int, pro.ant.rare)
pro.cor.plot[[1]]

euk.cor.plot <- get_my_cor_plots(euk.ant.abun, euk.ant.int, euk.ant.rare)
euk.cor.plot[[1]]


multip1 <- plot_grid(
  pro.cor.plot[[1]],
  euk.cor.plot[[1]],
  #align = 'v',
  labels = c("Prokaryote", "Eukaryote"),
  #rel_widths = c(1,1),
  align = 'v',
  ncol = 1
)
multip1

pdf("../results/alpha-correlation-results.pdf")
print(multip1)
dev.off()


print_my_results <- function(data, name){
  
  #output to print
  data.rho = as.data.frame(data$r)
  data.p = as.data.frame(data$p.adj)
  data.ci = as.data.frame(data$ci.adj)
  
  full_obs_ans = paste0("Spearman's correlation analysis was employed to explore significant correlations between environmental variables and observed ASVs within the ", name, ". ")
  full_shan_ans = paste0("Spearman's correlation analysis was employed to explore significant correlations between environmental variables and Shannon diversity within the ", name, ". ")
  
  for (i in 1:nrow(data.p)){
    
    if (data.p$Observed[[i]] < 0.05){
      
      ans = paste0("There is a signficant correlation between ", rownames(data.p)[[i]], " and number of observed ASVs [ρ=", round(data.rho$Observed[[i]], 2),
                   ", (", round(data.ci$lower.adj[[i]], 2), " to ", round(data.ci$upper.adj[[i]], 2),  " 95% CI), ", mod_my_p_val(data.p$Observed[[i]]), "]. ")
      
      full_obs_ans = paste0(full_obs_ans, ans)
      
    } 
    
    if (data.p$Shannon[[i]] < 0.05){
      
      ans = paste0("There is a signficant correlation between ", rownames(data.p)[[i]], " and Shannon diversity [ρ=", round(data.rho$Shannon[[i]], 2),
                   ", (", round(data.ci$lower.adj[[i+nrow(data.rho)]], 2), " to ", round(data.ci$upper.adj[[i+nrow(data.rho)]], 2),  " 95% CI), ", mod_my_p_val(data.p$Shannon[[i]]), "]. ")
      
      full_shan_ans = paste0(full_shan_ans, ans)
      
    }
  }
  
  all_res = paste0(full_obs_ans, full_shan_ans)
  
  return(all_res)
}


### Abundant Prokaryote Results

print.pro.abun = print_my_results(pro.cor.plot[[2]][[1]], "abundant prokaryote subcommunity")
`r print.pro.abun`

### Intermediate Prokaryote Results

print.pro.int = print_my_results(pro.cor.plot[[2]][[2]], "intermediate prokaryote subcommunity")
`r print.pro.int`

### Rare Prokaryote Results
print.pro.rare = print_my_results(pro.cor.plot[[2]][[3]], "rare prokaryote subcommunity")

`r print.pro.rare`

### Abundant Eukaryote Results
print.euk.abun = print_my_results(euk.cor.plot[[2]][[1]], "abundant eukaryote subcommunity")

`r print.euk.abun`

### Intermediate Eukaryote Results
print.euk.int = print_my_results(euk.cor.plot[[2]][[2]], "intermediate eukaryote subcommunity")

`r print.euk.int`

### Rare Eukaryote Results
print.euk.rare = print_my_results(euk.cor.plot[[2]][[3]], "rare eukaryote subcommunity")

`r print.euk.rare`



