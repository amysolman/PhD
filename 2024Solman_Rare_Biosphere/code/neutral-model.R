# Step 1: Clear workspace and install packages
rm(list=ls())
library(Hmisc) #binconf function
library(minpack.lm) #for non-linear model fitting
# library(stats4) #mle function
library(dplyr) #for %>% function
library(tibble) #function rownames_to_column
#install.packages("GUniFrac") #for Rarefy funtion
#library(GUniFrac)
library(phyloseq)
library(cowplot)
library(scales) #for percentages
library(stats4)
library(stringr)
library(ggplot2)

## Method

# Code developed from https://github.com/mrebolleda/OrganFilters_MimulusMicrobiome/blob/master/Neutral_pol_nopol.md and Burns et al., 2016
# 
# To determine the potential contribution of stochastic processes on the assemblage of abundant, intermediate and rare taxa, Sloan’s neutral community model (Sloan et al., 2006, 2007) was applied to predict the association between ASV observation frequency (i.e. the number of samples each ASV is observed in) and their relative abundance in the metacommunity. Code developed from Gómez & Ashman (2019) (https://github.com/mrebolleda/OrganFilters_MimulusMicrobiome/blob/master/Neutral_pol_nopol.md) was used for fitting the model. In a system shaped by stochastic processes of drift and dispersal, the model predicts taxa with greater relative abundance in the metacommunity will be observed in a greater number of sites, as a large population increases the likelihood of random dispersal between sites. Conversely, rare taxa are less likely to disperse between sites and are more prone to local elimination through ecological drift. With a single free parameter m (migration rate), the model predicts the frequency with which an ASV will be observed, considering it's mean relative abundance in the metacommunity. Akaike informaition criterion (AIC) was used to compare the fit of the neutral model to a binomial distribution, serving as a null model of random sampling from the metacommunity without drift or dispersal limitation. 95% confidence intervals were calculated around model predictions. ASVs within the 95% CI were considered to be distributed according to neutral expectations.
# 
# Step One: Extract the ASV table from the phyloseq object
# Step Two: get the mean number of reads per sample (mean sum of each row) e.g. individuals per community 
# Step Three: Get the mean relative abundance of each ASV
# Step Four: Calculate occurrence frequency of each ASV
# Step Five: Combine data into dataframe
# Step Six: Calculate the limit of detection
# Step Seven: Fit model parameter m (immigration rate) using Non-linear least squares (NLS)
# Step Eight: Fit model parameter m (immigration rate) using Maximum likelihood estimation (MLE)
# Step Nine: Calculate AIC and BIC for model fit
# Step Ten: Make N an integer for binomial model fitting
# Step Eleven: Calculate goodness of fit of model
# Step Twelve: Calculate AIC for binomial model
# Step Thirteen: Goodness of fit for binomial model
# Step Fourteen: Calculate AIC for poisson model
# Step Fifteen: Goodness of fit for poisson model
# Step Sixteen: Get all fitting statistics
# Step Seventeen: Get table of observed and prediced values


# Step 2: Read in data
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

# In the neutral model, we calculate the mean relative abundance of each ASV by:
# 1. Calculating the mean number of reads per sample
# 2. Calculating the mean number of reads in each sample for each ASV (total number of reads/total number of samples)
# 3. Dividing the mean number of reads for each ASV by the mean number of reads per sample
# 
# However, when dividing ASVs into our subcommunities we:
# 1. Divide the number of reads in each sample by the total number of reads, giving the proportion of that sample accounted for by each ASV
# 2. Calculate the mean of those proportions for each ASV across all samples
# 
# Subcommunity method
# Find mean relative abundance of each ASV
# Neutral method
# Get the mean number of reads per ASV in each sample (e.g. 5) and divide that by the mean number of reads per sample



# data = data.frame(otu_table(pro.ant))
# 
# #Method 1 - Subcommunity Method
# abundance <- data
# for (i in 1:ncol(data)) {
#   abundance[,i] <- abundance[,i] / sum(abundance[,i])
# }
# x = rowMeans(abundance) #mean relative abundance of each ASV
# 
# #Method 2 - Neutral Model Method
# N <- mean(apply(data, 2, sum)) #find the mean total number of reads per sample
# p.m <- apply(data, 1, mean) #get the mean number of reads per ASV
# p <- p.m/N #mean relative abundance of each ASV


neutral_model <- function(phylo){
  
#Step One: Extract the ASV table from the phyloseq object
ASV.table = data.frame(t(otu_table(phylo)), check.names = FALSE)

#Step Two: get the mean number of reads per sample (mean sum of each row) e.g. individuals per community 
N <- mean(apply(ASV.table, 1, sum))

#Step Three: Get the mean relative abundance of each ASV.
#Get the mean number of reads for each ASV across all samples (mean of each column)
p.m <- apply(ASV.table, 2, mean)
#Remove any zeros
p.m <- p.m[p.m != 0] #remove any ASVs with zero counts
#divide the number of mean reads by the total number of reads per sample (mean relative abundance of each ASV globally)
p <- p.m/N #mean relative abundance of each ASV


#Step Four: Calculate occurrence frequency of each ASV. 
#Make ASV.table into presence/absence table
ASV.table.bi <- 1*(ASV.table>0)
#find the mean of frequence of each column (ASV) e.g. the mean number of samples each ASV is found in
freq.table <- apply(ASV.table.bi, 2, mean) 
#only keep ASVs with a frequency other than 0
freq.table <- freq.table[freq.table != 0] #only keep data that isn't zero


#Step Five: Combine data into dataframe.
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
#get vectors of our mean relative abundances and frequencies
p <- C.no0$p #this creates a vector of the mean relative abundances of each ASV
freq <- C.no0$freq #this creates a list of the mean frequency of each ASV


#Step Six: Calculate the limit of detection
d <- 1/N


#Step Seven: Fit model parameter m (immigration rate) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001)) 
#get the confidence intervals for m
m.ci <- confint(m.fit, 'm', level=0.95)


#Step Eight: Fit model parameter m (immigration rate) using Maximum likelihood estimation (MLE)
# sncm.LL <- function(m, sigma){
# 		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
# 		R = dnorm(R, 0, sigma)
# 		-sum(log(R))
# }
# m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))


#Step Nine: Calculate AIC and BIC for model fit
# aic.fit = AIC(m.mle, k=2)
aic.fit = AIC(m.fit)
# bic.fit = BIC(m.mle)
bic.fit = BIC(m.fit)


#Step Ten: Make N an integer for binomial model fitting
N = as.integer(N)


#Step Eleven: Calculate goodness of fit of model
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE) #get the predicted ASV frequencies under stochastic processes
#Get R2 value
RSS <- sum(residuals(m.fit)^2) #Residual sum of squares of our NLLS model
TSS <- sum((freq - mean(freq))^2) #Total sum of squares 
Rsqr <- 1 - (RSS/TSS) #R-squared value
RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1)) #root mean squares errors
# Get confidence interval for predictions
freq.pred.ci <- binconf(freq.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)


#Step Twelve: Calculate AIC for binomial model
##Calculate AIC for binomial model
bino.LL <- function(mu, sigma){
  R = freq - pbinom(d, N, p, lower.tail=FALSE)
  R = dnorm(R, mu, sigma)
  -sum(log(R))
}
bino.mle <- mle(bino.LL, start=list(mu=0.1, sigma=0.1), nobs=length(p))

aic.bino <- AIC(bino.mle, k=2)
bic.bino <- BIC(bino.mle)


#Step Thirteen: Goodness of fit for binomial model
bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

bino.pred.ci <- binconf(bino.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)


#Step Fourteen: Calculate AIC for poisson model
pois.LL <- function(mu, sigma){
  R = freq - ppois(d, N*p, lower.tail=FALSE)
  R = dnorm(R, mu, sigma)
  -sum(log(R))
}
pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))

aic.pois <- AIC(pois.mle, k=2)
bic.pois <- BIC(pois.mle)


#Step Fifteen: Goodness of fit for poisson model
pois.pred <- ppois(d, N*p, lower.tail=FALSE)
Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

pois.pred.ci <- binconf(pois.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)


#Step Sixteen: Get all fitting statistics
fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(),
                       poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())

fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1],NA, NA, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(ASV.table), length(p), d)


#Step Seventeen: Get table of observed and prediced values
A <- cbind(p, freq, freq.pred, freq.pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
A <- as.data.frame(A)
colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
rownames(A) = C.no0$ASV			#add ASV IDS
pred.df <- A[order(A[,1]),]


res.list <- list(pred.df, fitstats)
}

#fit the models
pro.ant.out = neutral_model(pro.ant)
pro.arc.out = neutral_model(pro.arc)
euk.ant.out = neutral_model(euk.ant)
euk.arc.out = neutral_model(euk.arc)


#Find the percentage of total community that falls within the neutral model predictions


perc_within_neutral_prediction <- function(data){
  
  perc_check = data
  out_res = vector()
  for (i in 1:nrow(perc_check)){
    row_check = perc_check[i,]
    out_res[[i]] = row_check$freq > row_check$pred.lwr & row_check$freq < row_check$pred.upr
  }
  
  tab_res = table(out_res)
  my_res = tab_res[[2]]/nrow(perc_check)*100
  return(my_res)
  
}

pro.ant.perc = perc_within_neutral_prediction(pro.ant.out[[1]])
pro.arc.perc = perc_within_neutral_prediction(pro.arc.out[[1]])
euk.ant.perc = perc_within_neutral_prediction(euk.ant.out[[1]])
euk.arc.perc = perc_within_neutral_prediction(euk.arc.out[[1]])

###################################################################################################
###################################################################################################

#Report NEUTRAL MODEL RESULTS
sink("../results/neutral-model-results.txt", type="output")
writeLines("===============================================================
NEUTRAL MODEL RESULTS
===============================================================")
writeLines("The percentage of Antarctic prokaryote ASVs within the predictions of the neutral model:")
pro.ant.perc
writeLines("The percentage of Arctic prokaryote ASVs within the predictions of the neutral model:") 
pro.arc.perc
writeLines("The percentage of Antarctic eukaryote ASVs within the predictions of the neutral model:")
euk.ant.perc
writeLines("The percentage of Arctic eukaryote ASVs within the predictions of the neutral model:") 
euk.arc.perc
sink()

###################################################################################################
###################################################################################################

# res = pro.ant.out
# pole = "Ant"

one_results <- function(res, pole){
  
  #combine the fit results dataframes
  x = res[[2]]
  
  neutral_model = c(x$Rsqr, x$RMSE, x$AIC, x$BIC)
  binomial_model = c(x$Rsqr.bino, x$RMSE.bino, x$AIC.bino, x$BIC.bino)
  poisson_model = c(x$Rsqr.pois, x$RMSE.pois, x$AIC.pois, x$BIC.pois)
  
  x.df = data.frame(cbind(neutral_model, binomial_model, poisson_model))
  
  rownames(x.df) = c("R2", "RMSE", "AIC", "BIC")
  names(x.df) = c(paste(pole, " Neutral"), paste(pole, " Binomial"), paste(pole, " Poisson"))
  x.df = round(x.df, 3)
  x.df <- tibble::rownames_to_column(x.df, "Metric")
  
  return(x.df)
}

# res1 = pro.ant.out
# res2 = pro.arc.out

two_results <- function(res1, res2){
  
  
  out1 = one_results(res1, "Ant")
  out2 = one_results(res2, "Arc")
  
  fin = cbind(out1, out2)
  
  return(fin)
  
}

pro.res = two_results(pro.ant.out, pro.arc.out)
euk.res = two_results(euk.ant.out, euk.arc.out)

export_res = rbind(pro.res, euk.res)
export_res = cbind(Group=c(rep("Prokaryote", nrow(pro.res)), rep("Eukaryote", nrow(euk.res))), export_res)

write.csv(export_res, "../results/neutral-model-results.csv")



get_abun_df <- function(phylo1, phylo2, phylo3){
  
  abun_dt = data.frame(ASV=rownames(data.frame(otu_table(phylo1))), Abundance = rep("Abundant", length(rownames(data.frame(otu_table(phylo1))))))
  int_dt = data.frame(ASV=rownames(data.frame(otu_table(phylo2))), Abundance = rep("Intermediate", length(rownames(data.frame(otu_table(phylo2))))))
  rare_dt = data.frame(ASV=rownames(data.frame(otu_table(phylo3))), Abundance = rep("Rare", length(rownames(data.frame(otu_table(phylo3))))))
  df = rbind(abun_dt, int_dt, rare_dt)
  
  return(df)
}

  
  #Antarctic prokaryotes
  df = get_abun_df(pro.ant.abun, pro.ant.int, pro.ant.rare)
  res = neutral_model(pro.ant)
  df2plot = res[[1]]
  num2plot1 = res[[2]]
  #add abundance classifications
  df2plot$ASV <- rownames(df2plot)
  df2plot1 <- merge(df2plot,df, by=c("ASV"))
  #add dataset info
  df2plot1$Group = "Prokaryote"
  df2plot1$Pole =  "Antarctic"
  num2plot1$Group = "Prokaryote"
  num2plot1$Pole =  "Antarctic"
  
  #Arctic prokaryotes
  df = get_abun_df(pro.arc.abun, pro.arc.int, pro.arc.rare)
  res = neutral_model(pro.arc)
  df2plot = res[[1]]
  num2plot2 = res[[2]]
  #add abundance classifications
  df2plot$ASV <- rownames(df2plot)
  df2plot2 <- merge(df2plot,df, by=c("ASV"))
  #add dataset info
  df2plot2$Group = "Prokaryote"
  df2plot2$Pole =  "Arctic"
  num2plot2$Group = "Prokaryote"
  num2plot2$Pole =  "Arctic"
  
  #Antarctic eukaryotes
  df = get_abun_df(euk.ant.abun, euk.ant.int, euk.ant.rare)
  res = neutral_model(euk.ant)
  df2plot = res[[1]]
  num2plot3 = res[[2]]
  #add abundance classifications
  df2plot$ASV <- rownames(df2plot)
  df2plot3 <- merge(df2plot,df, by=c("ASV"))
  #add dataset info
  df2plot3$Group = "Eukaryote"
  df2plot3$Pole =  "Antarctic"
  num2plot3$Group = "Eukaryote"
  num2plot3$Pole =  "Antarctic"
  
  #Arctic prokaryotes
  df = get_abun_df(euk.arc.abun, euk.arc.int, euk.arc.rare)
  res = neutral_model(euk.arc)
  df2plot = res[[1]]
  num2plot4 = res[[2]]
  #add abundance classifications
  df2plot$ASV <- rownames(df2plot)
  df2plot4 <- merge(df2plot,df, by=c("ASV"))
  #add dataset info
  df2plot4$Group = "Eukaryote"
  df2plot4$Pole =  "Arctic"
  num2plot4$Group = "Eukaryote"
  num2plot4$Pole =  "Arctic"
  
  #bind the dataframes together
  final.df = rbind(df2plot1, df2plot2, df2plot3, df2plot4)
  final.num = rbind(num2plot1, num2plot2, num2plot3, num2plot4)
  
  
  final.df$Abundance = factor(final.df$Abundance, levels = c("Rare", "Intermediate", "Abundant"))
  final.df$Group = factor(final.df$Group, levels = c("Prokaryote", "Eukaryote"))
  #final.num$Abundance = factor(final.num$Abundance, levels = c("Rare", "Intermediate", "Abundant"))
  final.num$Group = factor(final.num$Group, levels = c("Prokaryote", "Eukaryote"))
  
  p1 = ggplot(data=final.df) +
    geom_point(data=final.df, aes(x=log10(p), y=freq, fill=Abundance, shape=Abundance),
               alpha=.5, size=5) +
    #geom_line(data=df2plot2, aes(x=log10(p), y=freq), color="black") +
    geom_line(data=final.df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black", size=1) +
    geom_line(data=final.df, aes(x=log10(p), y=pred.upr), linetype=2, color="black", size=1) +
    geom_text(data=final.num, aes(label = paste("R^2 == ", round(Rsqr, 3))), x=-4.9, y=0.8, size=7, parse=TRUE) +
    geom_text(data=final.num, aes(label = paste("italic(m) ==", round(m, 4))), x=-4.9, y=0.68, size=7, parse=TRUE) +
    scale_shape_manual(values = c(21, 21, 21))+
    scale_fill_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    scale_color_manual(values=c("#82ca81", "#74a9d8", "#ef7a76"))+
    facet_wrap(~Group + Pole)+
    labs(x="log10(mean relative abundance)", y="Percentage occurence frequency") +
    xlim(-6, 4)+
    theme_bw() +
    theme(axis.line = element_line(color="black"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size=18),
          axis.title = element_text(size=16),
          axis.text = element_text(size=16),
          strip.text = element_text(size=16))+
    guides(fill = guide_legend(override.aes = list(size = 10, alpha=.8) ) )

  print(p1)
  

pdf("../results/neutral-model.pdf", width=16, height=14)
print(p1)
dev.off()



above_below <- function(phylo1, phylo2, phylo3, phylo4){
  
  df = get_abun_df(phylo2, phylo3, phylo4)
  res = neutral_model(phylo1)
  prediction_results <- res[[1]]
  prediction_results$ASV <- rownames(prediction_results)
  out.df = merge(prediction_results,df, by=c("ASV"))
  
  #divide dataframe by abundance class
  pred_res_abun <- out.df %>% filter(Abundance == "Abundant")
  pred_res_int <- out.df %>% filter(Abundance == "Intermediate")
  pred_res_rare <- out.df %>% filter(Abundance == "Rare")
  
  
  #What percentage of frequency are between the CI?
  abun <- vector()
  for (k in 1:nrow(pred_res_abun)){
    if (pred_res_abun$freq[k] > pred_res_abun$pred.lwr[k] && pred_res_abun$freq[k] < pred_res_abun$pred.upr[k]){
      abun <- c(abun, "Neutral distribution")
    } else if (pred_res_abun$freq[k] > pred_res_abun$pred.lwr[k] && pred_res_abun$freq[k] > pred_res_abun$pred.upr[k]) {
      abun <- c(abun, "Above prediction")
    } else {
      abun <- c(abun, "Below prediction")
    }
  }
  abun_tab <- data.frame(table(abun))
  abun_tab$Perc = round(abun_tab$Freq/sum(abun_tab$Freq)*100, 2)
  abun_tab$Abundance <- "Abundant"
  names(abun_tab) <- c("Prediction","Num", "Perc", "Abundance")
  
  pred_res_abun_position = cbind(pred_res_abun, abun)
  names(pred_res_abun_position) = c("ASV", "p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr", "Abundance", "Prediction")
  
  int <- vector()
  for (l in 1:nrow(pred_res_int)){
    if (pred_res_int$freq[l] > pred_res_int$pred.lwr[l] && pred_res_int$freq[l] < pred_res_int$pred.upr[l]){
      int <- c(int, "Neutral distribution")
    } else if (pred_res_int$freq[l] > pred_res_int$pred.lwr[l] && pred_res_int$freq[l] > pred_res_int$pred.upr[l]) {
      int <- c(int, "Above prediction")
    } else {
      int <- c(int, "Below prediction")
    }
  }
  int_tab <- data.frame(table(int))
  int_tab$Perc = round(int_tab$Freq/sum(int_tab$Freq)*100, 2)
  int_tab$Abundance <- "Intermediate"
  names(int_tab) <- c("Prediction","Num",  "Perc", "Abundance")
  
  pred_res_int_position = cbind(pred_res_int, int)
  names(pred_res_int_position) = c("ASV", "p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr", "Abundance", "Prediction")
  
  rare <- vector()
  for (m in 1:nrow(pred_res_rare)){
    if (pred_res_rare$freq[m] > pred_res_rare$pred.lwr[m] && pred_res_rare$freq[m] < pred_res_rare$pred.upr[m]){
      rare <- c(rare, "Neutral distribution")
    } else if (pred_res_rare$freq[m] > pred_res_rare$pred.lwr[m] && pred_res_rare$freq[m] > pred_res_rare$pred.upr[m]) {
      rare <- c(rare, "Above prediction")
    } else {
      rare <- c(rare, "Below prediction")
    }
  }
  rare_tab <- data.frame(table(rare))
  rare_tab$Perc = round(rare_tab$Freq/sum(rare_tab$Freq)*100, 2)
  rare_tab$Abundance <- "Rare"
  names(rare_tab) <- c("Prediction","Num", "Perc","Abundance")
  
  pred_res_rare_position = cbind(pred_res_rare, rare)
  names(pred_res_rare_position) = c("ASV", "p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr", "Abundance", "Prediction")
  
  #combin dataframes
  neutral_res <- rbind(abun_tab, int_tab, rare_tab)
  plot_res <- rbind(pred_res_abun_position, 
                    pred_res_int_position,
                    pred_res_rare_position)
  
  res.list = list(neutral_res, plot_res)
  
  return(res.list)
  
}



  
  #Antarctic prokaryotes
  pred1 = above_below(pro.ant, pro.ant.abun, pro.ant.int, pro.ant.rare)[[2]]
  res1 = neutral_model(pro.ant)
  num2plot1 = res1[[2]]
  #add dataset info
  pred1$Group = "Prokaryote"
  pred1$Pole =  "Antarctic"
  num2plot1$Group = "Prokaryote"
  num2plot1$Pole =  "Antarctic"
  
  #Arctic prokaryotes
  pred2 = above_below(pro.arc, pro.arc.abun, pro.arc.int, pro.arc.rare)[[2]]
  res2 = neutral_model(pro.arc)
  num2plot2 = res2[[2]]
  #add dataset info
  pred2$Group = "Prokaryote"
  pred2$Pole =  "Arctic"
  num2plot2$Group = "Prokaryote"
  num2plot2$Pole =  "Arctic"
  
  #Antarctic eukaryotes
  pred3 = above_below(euk.ant, euk.ant.abun, euk.ant.int, euk.ant.rare)[[2]]
  res3 = neutral_model(euk.ant)
  num2plot3 = res3[[2]]
  #add dataset info
  pred3$Group = "Eukaryote"
  pred3$Pole =  "Antarctic"
  num2plot3$Group = "Eukaryote"
  num2plot3$Pole =  "Antarctic"
  
  #Arctic eukaryotes
  pred4 = above_below(euk.arc, euk.arc.abun, euk.arc.int, euk.arc.rare)[[2]]
  res4 = neutral_model(euk.arc)
  num2plot4 = res4[[2]]
  #add dataset info
  pred4$Group = "Eukaryote"
  pred4$Pole =  "Arctic"
  num2plot4$Group = "Eukaryote"
  num2plot4$Pole =  "Arctic"
  
  #bind the dataframes together
  final.pred = rbind(pred1, pred2, pred3, pred4)
  final.num = rbind(num2plot1, num2plot2, num2plot3, num2plot4)
  
  
  final.pred$Abundance = factor(final.pred$Abundance, levels = c("Rare", "Intermediate", "Abundant"))
  final.pred$Group = factor(final.pred$Group, levels = c("Prokaryote", "Eukaryote"))
  final.pred$Prediction = factor(final.pred$Prediction, levels = c("Above prediction", "Neutral distribution", "Below prediction"))
  #final.num$Abundance = factor(final.num$Abundance, levels = c("Rare", "Intermediate", "Abundant"))
  final.num$Group = factor(final.num$Group, levels = c("Prokaryote", "Eukaryote"))
  
  
  p2 = ggplot(data=final.pred) +
    geom_point(data=final.pred, aes(x=log10(p), y=freq, fill=Prediction, shape=Prediction),
               alpha=.5, size=3) +
    #geom_line(data=df2plot2, aes(x=log10(p), y=freq), color="black") +
    geom_line(data=final.pred, aes(x=log10(p), y=pred.lwr), linetype=2, color="black", size=1) +
    geom_line(data=final.pred, aes(x=log10(p), y=pred.upr), linetype=2, color="black", size=1) +
    geom_text(data=final.num, aes(label = paste("R^2 == ", round(Rsqr, 3))), x=-4.9, y=0.8, size=4, parse=TRUE) +
    geom_text(data=final.num, aes(label = paste("italic(m) ==", round(m, 4))), x=-4.9, y=0.68, size=4, parse=TRUE) +
    labs(x="log10(mean relative abundance)", y="Percentage occurence frequency") +
    xlim(-6, -1)+
    facet_wrap(~Group + Pole)+
    theme_bw() +
    theme(axis.line = element_line(color="black"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          axis.title = element_text(size=12),
          axis.text = element_text(size=10))+
    #guides(colour=guide_legend(title="Abundance Classification"))+
    scale_color_manual(values=c("#f8c4be", "#88d2d6", "#f18579"))+
    scale_fill_manual(values=c("#f8c4be", "#88d2d6", "#f18579"))+
  scale_shape_manual(values = c(21, 21, 21))+
    guides(fill = guide_legend(override.aes = list(size = 7, alpha=.8) ) )
  
  p2
  
  pdf("../results/neutral-model-2.pdf")
  print(p2)
  dev.off()

pro.ant.pred = above_below(pro.ant, pro.ant.abun, pro.ant.int, pro.ant.rare)[[1]]
pro.ant.pred$Group = "Prokaryote"
pro.ant.pred$Pole =  "Antarctic"
pro.arc.pred = above_below(pro.arc, pro.arc.abun, pro.arc.int, pro.arc.rare)[[1]]
pro.arc.pred$Group = "Prokaryote"
pro.arc.pred$Pole =  "Arctic"
euk.ant.pred = above_below(euk.ant, euk.ant.abun, euk.ant.int, euk.ant.rare)[[1]]
euk.ant.pred$Group = "Eukaryote"
euk.ant.pred$Pole =  "Antarctic"
euk.arc.pred = above_below(euk.arc, euk.arc.abun, euk.arc.int, euk.arc.rare)[[1]]
euk.arc.pred$Group = "Eukaryote"
euk.arc.pred$Pole =  "Arctic"

#get percentages table
#pro.ant.pred

df.2.plot = rbind(pro.arc.pred, pro.ant.pred, euk.arc.pred, euk.ant.pred)
df.2.plot$Prediction = as.factor(df.2.plot$Prediction)

df.2.plot$Abundance = factor(df.2.plot$Abundance, levels = c("Rare", "Intermediate", "Abundant"))
df.2.plot$Group = factor(df.2.plot$Group, levels = c("Prokaryote", "Eukaryote"))

p3 = ggplot(df.2.plot,                         
            aes(x = Abundance,
                y = Perc,
                fill = factor(Prediction, levels=c("Above prediction", "Neutral distribution", "Below prediction")))) + 
  geom_bar(stat = "identity",
           position = "fill") +
  scale_y_continuous(labels=scales::percent)+
  geom_text(
    aes(label = paste(Perc, "%")),
    position = position_fill(0.5), size=5
  )+
  facet_grid(~ Group + Pole)+
  theme_bw()+
  ylab("% ASVs")+
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        axis.title.x=element_blank(), 
        axis.text=element_text(size=10), 
        legend.text = element_text(size=14),
        #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
        )+
  scale_fill_manual(values=c("#f8c4be", "#88d2d6","#f18579"))


p3

pdf("../results/neutral-model-barplot.pdf", width=12, height=8)
print(p3)
dev.off()