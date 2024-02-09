# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

source("00-solman-functions.R")
library(phyloseq)
library(ggplot2)
library(vegan) #for vegdist
library(dplyr) #for %>% function
library(cowplot) #plot_grid
library(Hmisc) #binconf function
library(minpack.lm) #for non-linear model fitting
library(tibble) #function rownames_to_column
library(stats4) #mle function
library(patchwork)
library(egg)
library(ggrepel)
library(tidyverse)
library(scales)


#prokaryotes
ps.pro <- readRDS("../results/16S-ps-controls-removed.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results/18S-ps-no-mm-controls-removed.rds") 
#micrometazoans only
ps.mm <- readRDS("../results/18S-ps-mm-only-controls-removed.rds") 

#remove the sample low read count samples as we did previously

#prokaryotes
ps.p.pro = prune_samples(sample_sums(ps.pro)>= 1900, ps.pro) #retain samples with >= num counts
ps.f.pro = filter_taxa(ps.p.pro, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts
#eukaryotes
ps.p.euk = prune_samples(sample_sums(ps.euk)>= 1900, ps.euk) #retain samples with >= num counts
ps.f.euk = filter_taxa(ps.p.euk, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts
#micrometazoans
ps.p.mm = prune_samples(sample_sums(ps.mm)>= 500, ps.mm) #retain samples with >= num counts
ps.f.mm = filter_taxa(ps.p.mm, function(x) sum(x) >= 1, TRUE) #remove ASVs with zero counts


neutral_model <- function(ps, hab){
  
  #Separate habitat data
  h <- prune_samples(hab, ps)
  
  #remove empty samples or ASVs
  h.f = phyloseq::filter_taxa(h, function(x) sum(x) >0, TRUE)
  h.f = phyloseq::prune_samples(sample_sums(h.f)>0, h.f)
  
  #Step One: Extract the ASV table from the phyloseq object
  ASV.table = data.frame(t(otu_table(h.f)), check.names = FALSE)
  
  #Step Two: get the mean number of reads per sample (mean sum of each row) e.g. individuals per community 
  N <- mean(apply(ASV.table, 1, sum))
  
  #Step Three: Get the mean abundance of each ASV.
  #Get the mean number of reads for each ASV across all samples (mean of each column)
  p.m <- apply(ASV.table, 2, mean)
  
  #Remove any zeros
  p.m <- p.m[p.m != 0] #remove any ASVs with zero counts
  #divide the number of mean reads by the total number of reads per sample (mean relative abundance of each ASV globally)
  
  p <- p.m/N #mean relative abundance of each ASV
  sum(p) == 1 #should be true
  
  
  #Step Four: Calculate occurrence frequency of each ASV. 
  #Make ASV.table into presence/absence table
  ASV.table.bi <- 1*(ASV.table>0)
  
  #find the mean of frequency of each column (ASV) e.g. the mean number of samples each ASV is found in
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
  
  pbeta(1, 2, 3, lower.tail=FALSE)
  
  #get the confidence intervals for m
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  
  #Step Eight: Fit model parameter m (immigration rate) using Maximum likelihood estimation (MLE)
  # sncm.LL <- function(m, sigma){
  # 		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
  # 		R = dnorm(R, 0, sigma)
  # 		-sum(log(R))
  # }
  # 
  # m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  
  #Step Nine: Calculate AIC and BIC for model fit
  #  aic.fit = AIC(m.mle, k=2)
  aic.fit = AIC(m.fit)
  # bic.fit = BIC(m.mle)
  bic.fit = BIC(m.fit)
  
  
  #Step Ten: Make N an integer for binomial model fitting
  N = as.integer(N)
  
  
  #Step Eleven: Calculate goodness of fit of model
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE) #get the predicted ASV frequencies under stochastic processes
  
  #looking at the model
  #the model takes two vectors of values, that means a p parameter and q parameter is given for each ASV
  #using these parameters for each ASV the model predicts the proportion of sites it will be observed in
  freq.pred # the percentage of sites each ASV will be observed in
  a = N*coef(m.fit)*p #the probability of each ASV appearing in a particular site
  b = N*coef(m.fit)*(1-p) #the probability of each ASV not appearing in a particular site
  
  pbeta(0.0041, c(1,2,1,6,7), c(3,4,2,38))
  
  
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
  
  return(res.list)
}




#fit the models
pro.snow <- neutral_model(ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID)
pro.sp <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
pro.sum <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID)
pro.cry <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID)

euk.snow <- neutral_model(ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID)
euk.sp <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID)
euk.sum <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID)
euk.cry <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID)

#Find the percentage of total community that falls within the neutral model predictions


# data = mm.sum[[1]]

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


pro.snow.perc = perc_within_neutral_prediction(pro.snow[[1]])
pro.sp.perc = perc_within_neutral_prediction(pro.sp[[1]])
pro.sum.perc = perc_within_neutral_prediction(pro.sum[[1]])
pro.cry.perc = perc_within_neutral_prediction(pro.cry[[1]])

euk.snow.perc = perc_within_neutral_prediction(euk.snow[[1]])
euk.sp.perc = perc_within_neutral_prediction(euk.sp[[1]])
euk.sum.perc = perc_within_neutral_prediction(euk.sum[[1]])
euk.cry.perc = perc_within_neutral_prediction(euk.cry[[1]])

#get the results for our model fittings for a single habitat

one_results <- function(res, hab){
  
  #combine the fit results dataframes
  x = res[[2]]
  
  neutral_model = c(x$Rsqr, x$RMSE, x$AIC, x$BIC)
  binomial_model = c(x$Rsqr.bino, x$RMSE.bino, x$AIC.bino, x$BIC.bino)
  poisson_model = c(x$Rsqr.pois, x$RMSE.pois, x$AIC.pois, x$BIC.pois)
  
  x.df = data.frame(cbind(neutral_model, binomial_model, poisson_model))
  
  rownames(x.df) = c("R2", "RMSE", "AIC", "BIC")
  names(x.df) = c(paste(hab, " Neutral"), paste(hab, " Binomial"), paste(hab, " Poisson"))
  x.df = round(x.df, 3)
  x.df <- tibble::rownames_to_column(x.df, "Metric")
  
  return(x.df)
}


#get the model fitting results for all habitats and combine

four_results <- function(res1, res2, res3, res4){
  
  
  out1 = one_results(res1, "Snow")
  out2 = one_results(res2, "Spring Ice")
  out3 = one_results(res3, "Summer Ice")
  out4 = one_results(res4, "Cryoconite")
  
  fin = cbind(out1, out2[,2:4], out3[,2:4], out4[,2:4])
  
  return(fin)
  
}

pro.df = four_results(pro.snow, pro.sp, pro.sum, pro.cry)
euk.df = four_results(euk.snow, euk.sp, euk.sum, euk.cry)


#save and export the dataframe
mod.fit.res = rbind(pro.df, euk.df)
write.csv(mod.fit.res, "../results/neutral-model-fit-results.csv") 


###NEUTRAL COLOUR PLOT#####

neutral_colour_plot <- function(ps, hab, leg, x.axis, y.axis, habitat, inset, group){

  res = neutral_model(ps, hab)
  df <- res[[1]]
  #prediction_results$ASV <- rownames(prediction_results)

  #What percentage of frequency are between the CI?
  x <- vector()
  for (k in 1:nrow(df)){
    if (df$freq[k] > df$pred.upr[k]) {
      x <- c(x, "Above prediction")
    } else if (df$freq[k] < df$pred.lwr[k]){
      x <- c(x, "Below prediction")
    } else {
      x <- c(x, "Neutral prediction")
    }
  }

  table(x)
  # asv_tab <- data.frame(table(asv))
  # asv_tab$Perc = round(asv_tab$Freq/sum(asv_tab$Freq)*100, 2)
  # names(asv_tab) <- c("Prediction","Num", "Perc")
  #
  # pred_res_position = cbind(prediction_results, asv)
  # names(pred_res_position) = c("p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr",  "ASV", "Prediction")

  df$Prediction = x

  #label for plotting
  num2plot = res[[2]]
  df.annotations <- data.frame(
    label = c(
      paste0("~R^{2} == ", round(num2plot$Rsqr, 2), "~italic(m) ==", round(num2plot$m,4))
    )
  )
  vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)




  if (leg == "Yes" && x.axis == "Yes" && y.axis == "Yes"){
    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            #legend.position = "none",
            axis.title = element_text(size=10),
            axis.text = element_text(size=10),
            legend.position="bottom")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)
  } else if (leg == "Yes" && x.axis == "No" && y.axis == "Yes"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title.x = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position="bottom")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)
  } else if (leg == "Yes" && x.axis == "Yes" && y.axis == "No"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position="bottom")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)




  }else if (leg == "Yes" && x.axis == "No" && y.axis == "No"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title = element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            legend.position="bottom")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)


  }else if (leg == "No"  && x.axis == "Yes" && y.axis == "Yes"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title = element_text(size=10),
            axis.text = element_text(size=10),
            legend.position="none")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)

  } else if (leg == "No"  && x.axis == "No" && y.axis == "Yes"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title.x = element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.position = "none")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)
    print(p1)

  } else if (leg == "No"  && x.axis == "Yes" && y.axis == "No"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)
    print(p1)



  }else if (leg == "No"  && x.axis == "No" && y.axis == "No"){

    p1 = ggplot(data=df) +
      geom_point(data=df, aes(x=log10(p), y=freq, color=Prediction),
                 alpha=.5, size=1) +
      geom_line(data=df, aes(x=log10(p), y=freq.pred), color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
      geom_line(data=df, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
      geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
                hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
      labs(x="log10(mean relative abundance)", y="Occurence frequency") +
      xlim(-6.3, 1)+
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            axis.title = element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            legend.position = "none")+
      # guides(colour=guide_legend(title="Abundance Classification"))+
      #scale_color_manual(values=c("#46C010","#808080","#d4AA00"))
      scale_colour_brewer(palette="Set1")+
      ggtitle(habitat)

    print(p1)

  }

  if (inset == "Yes"){

    #create inset plot
    pie.df = data.frame(table(x))
    pie = ggplot(pie.df, aes(x="", y=Freq, fill=x))+
      geom_bar(width = 1, stat = "identity")+
      coord_polar("y", start=0)+
      theme_minimal()+
      scale_fill_brewer(palette="Set1")+
      theme(legend.position = "none", axis.title = element_blank(),
            axis.text = element_blank())
    pie


    if (habitat == "Summer Ice" && group == "Eukaryotes" | habitat == "Cryoconite" && group == "Eukaryotes"){

      p1 <-
        ggdraw() +
        draw_plot(p1) +
        draw_plot(pie, x = 0.6, y = .15, width = .38, height = .7)

    } else if (habitat == "Snow" && group == "Eukaryotes" | habitat == "Spring Ice" && group == "Eukaryotes"){

      p1 <-
        ggdraw() +
        draw_plot(p1) +
        draw_plot(pie, x = 0.6, y = .1, width = .38, height = .8)

    } else{
      #combine plots
      p1 <-
        ggdraw() +
        draw_plot(p1) +
        draw_plot(pie, x = 0.6, y = .12, width = .38, height = .7)
    }

  }


  res.list = list(p1, table(x))

  return(res.list)

}




#get our plots
pro.snow.plot <- neutral_colour_plot(ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID, leg="No", x.axis = "No", y.axis="Yes", habitat = "Snow", inset = "Yes", group="Prokaryotes")
pro.snow.plot

pro.sp.plot <- neutral_colour_plot(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID, leg="No", x.axis = "No", y.axis="No", habitat = "Spring Ice", inset = "Yes", group="Prokaryotes")
pro.sp.plot

pro.sum.plot <- neutral_colour_plot(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID, leg="No", x.axis = "No", y.axis="Yes", habitat = "Summer Ice", inset = "Yes", group="Prokaryotes")
pro.sum.plot

pro.cry.plot <- neutral_colour_plot(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID, leg="No", x.axis = "No", y.axis="No", habitat = "Cryoconite", inset = "Yes", group="Prokaryotes")
pro.cry.plot

euk.snow.plot <- neutral_colour_plot(ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID, leg="No", x.axis = "No", y.axis="Yes", habitat = "Snow", inset = "Yes", group="Eukaryotes")
euk.snow.plot

euk.sp.plot <- neutral_colour_plot(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID, leg="No", x.axis = "No", y.axis="No", habitat = "Spring Ice", inset = "Yes", group="Eukaryotes")
euk.sp.plot

euk.sum.plot <- neutral_colour_plot(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID, leg="No", x.axis = "Yes", y.axis="Yes", habitat = "Summer Ice", inset = "Yes", group="Eukaryotes")
euk.sum.plot

euk.cry.plot <- neutral_colour_plot(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID, leg="No", x.axis = "Yes", y.axis="No", habitat = "Cryoconite", inset = "Yes", group="Eukaryotes")
euk.cry.plot

####PREP FINAL NEUTRAL PLOT###

#Snow prokryotes
res1 = neutral_model(ps.f.pro, data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID)
df1 <- res1[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df1)){
  if (df1$freq[k] > df1$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df1$freq[k] < df1$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df1$Prediction = x
#label for plotting
num2plot = res1[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num2plot$Rsqr, 2), "~italic(m)==", round(num2plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p1 = ggplot(data=df1) +
  geom_point(data=df1, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df1, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df1, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df1, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Snow")

print(p1)


####get points for putting labels
pie.df1 = data.frame(table(x))
pie.df1$value = round(pie.df1$Freq/sum(pie.df1$Freq)*100, 1)
pie.df1$group = pie.df1$x

pie.df1b <- pie.df1 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df1b <- pie.df1b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie1 = ggplot(pie.df1, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  #scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df1b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie1

p1p <- p1 + annotation_custom(ggplotGrob(pie1), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p1p)

#Spring Ice prokryotes
res2 = neutral_model(ps.f.pro, data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
df2 <- res2[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df2)){
  if (df2$freq[k] > df2$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df2$freq[k] < df2$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df2$Prediction = x
#label for plotting
num2plot = res2[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num2plot$Rsqr, 2), "~italic(m)==", round(num2plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p2 = ggplot(data=df2) +
  geom_point(data=df2, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df2, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df2, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df2, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Spring Ice")

print(p2)

####get points for putting labels
pie.df2 = data.frame(table(x))
pie.df2$value = round(pie.df2$Freq/sum(pie.df2$Freq)*100, 1)
pie.df2$group = pie.df2$x

pie.df2b <- pie.df2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df2b <- pie.df2b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie2 = ggplot(pie.df2, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df2b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie2

p2p <- p2 + annotation_custom(ggplotGrob(pie2), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p2p)

#Summer Ice prokryotes
res3 = neutral_model(ps.f.pro, data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID)
df3 <- res3[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df3)){
  if (df3$freq[k] > df3$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df3$freq[k] < df3$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df3$Prediction = x
#label for plotting
num3plot = res3[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num3plot$Rsqr, 2), "~italic(m)==", round(num3plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p3 = ggplot(data=df3) +
  geom_point(data=df3, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df3, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df3, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df3, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size=7),
        axis.text=element_text(size=7),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Summer Ice")

print(p3)

####get points for putting labels
pie.df3 = data.frame(table(x))
pie.df3$value = round(pie.df3$Freq/sum(pie.df3$Freq)*100, 1)
pie.df3$group = pie.df3$x

pie.df3b <- pie.df3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df3b <- pie.df3b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie3 = ggplot(pie.df3, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df3b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie3

p3p <- p3 + annotation_custom(ggplotGrob(pie3), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p3p)


#Cryoconite prokryotes
res4 = neutral_model(ps.f.pro, data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID)
df4 <- res4[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df4)){
  if (df4$freq[k] > df4$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df4$freq[k] < df4$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df4$Prediction = x
#label for plotting
num4plot = res4[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num4plot$Rsqr, 2), "~italic(m)==", round(num4plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p4 = ggplot(data=df4) +
  geom_point(data=df4, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df4, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df4, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df4, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Cryoconite")

print(p4)

####get points for putting labels
pie.df4 = data.frame(table(x))
pie.df4$value = round(pie.df4$Freq/sum(pie.df4$Freq)*100, 1)
pie.df4$group = pie.df4$x

pie.df4b <- pie.df4 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df4b <- pie.df4b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie4 = ggplot(pie.df4, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df4b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie4

p4p <- p4 + annotation_custom(ggplotGrob(pie4), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p4p)

#Snow eukaryotes
res5 = neutral_model(ps.f.euk, data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID)
df5 <- res5[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df5)){
  if (df5$freq[k] > df5$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df5$freq[k] < df5$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df5$Prediction = x
#label for plotting
num2plot = res5[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num2plot$Rsqr, 2), "~italic(m)==", round(num2plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p5 = ggplot(data=df5) +
  geom_point(data=df5, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df5, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df5, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df5, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        axis.title = element_text(size=7),
        axis.text = element_text(size=7),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Snow")

print(p5)

#add pie chart
pie.df5 = data.frame(table(x))
pie.df5$value = round(pie.df5$Freq/sum(pie.df5$Freq)*100, 1)
pie.df5$group = pie.df5$x

pie.df5b <- pie.df5 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df5b <- pie.df5b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie5 = ggplot(pie.df5, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df5b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie5

p5p <- p5 + annotation_custom(ggplotGrob(pie5), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p5p)


#Spring Ice eukkryotes
res6 = neutral_model(ps.f.euk, data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID)
df6 <- res6[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df6)){
  if (df6$freq[k] > df6$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df6$freq[k] < df6$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df6$Prediction = x
#label for plotting
num6plot = res6[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num6plot$Rsqr, 2), "~italic(m)==", round(num6plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p6 = ggplot(data=df6) +
  geom_point(data=df6, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df6, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df6, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df6, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Spring Ice")

print(p6)

#add pie chart
pie.df6 = data.frame(table(x))
pie.df6$value = round(pie.df6$Freq/sum(pie.df6$Freq)*100, 1)
pie.df6$group = pie.df6$x

pie.df6b <- pie.df6 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df6b <- pie.df6b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie6 = ggplot(pie.df6, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df6b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie6

p6p <- p6 + annotation_custom(ggplotGrob(pie6), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p6p)

#Summer Ice eukryotes
res7 = neutral_model(ps.f.euk, data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID)
df7 <- res7[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df7)){
  if (df7$freq[k] > df7$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df7$freq[k] < df7$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df7$Prediction = x
#label for plotting
num7plot = res7[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num7plot$Rsqr, 2), "~italic(m)==", round(num7plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p7 = ggplot(data=df7) +
  geom_point(data=df7, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df7, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df7, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df7, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Summer Ice")

print(p7)

#add pie chart
pie.df7 = data.frame(table(x))
pie.df7$value = round(pie.df7$Freq/sum(pie.df7$Freq)*100, 1)
pie.df7$group = pie.df7$x

pie.df7b <- pie.df7 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df7b <- pie.df7b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie7 = ggplot(pie.df7, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df7b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie7

p7p <- p7 + annotation_custom(ggplotGrob(pie7), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p7p)


#Cryoconite eukkryotes
res8 = neutral_model(ps.f.euk, data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID)
df8 <- res8[[1]]
#What percentage of frequency are between the CI?
x <- vector()
for (k in 1:nrow(df8)){
  if (df8$freq[k] > df8$pred.upr[k]) {
    x <- c(x, "Above prediction")
  } else if (df8$freq[k] < df8$pred.lwr[k]){
    x <- c(x, "Below prediction")
  } else {
    x <- c(x, "Neutral prediction")
  }
}
#add prediction to plotting dataframe
df8$Prediction = x
#label for plotting
num8plot = res8[[2]]
df.annotations <- data.frame(label = c(paste0("~italic(R^{2})==", round(num8plot$Rsqr, 2), "~italic(m)==", round(num8plot$m,4))))
vertical_adjustment = ifelse(grepl("\\bR\\b", df.annotations$label), 1.5, 3)

#plot the results
p8 = ggplot(data=df8) +
  geom_point(data=df8, aes(x=log10(p), y=freq, color=Prediction),
             alpha=.5, size=1) +
  geom_line(data=df8, aes(x=log10(p), y=freq.pred), color="black") +
  geom_line(data=df8, aes(x=log10(p), y=pred.lwr), linetype=2, color="black") +
  geom_line(data=df8, aes(x=log10(p), y=pred.upr), linetype=2, color="black") +
  geom_text(data = df.annotations, aes(x=-Inf, y=+Inf, label=label),
            hjust = 0, vjust = vertical_adjustment, size=4, parse = TRUE)+
  labs(x="log10(mean relative abundance)", y="Occurence frequency") +
  xlim(-6.3, 1)+
  theme_bw() +
  theme(axis.line = element_line(color="black"),
        legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=7),
        axis.title = element_text(size=7),
        plot.title = element_text(size=10))+
  #scale_colour_brewer(palette="Set1")+
  scale_color_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  ggtitle("Cryoconite")

print(p8)

#create inset plot
pie.df8 = data.frame(table(x))
pie.df8$value = round(pie.df8$Freq/sum(pie.df8$Freq)*100, 1)
pie.df8$group = pie.df8$x

pie.df8b <- pie.df8 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1))
pie.df8b <- pie.df8b %>%
  mutate(pos = dplyr::if_else(is.na(pos), value/2, pos))

pie8 = ggplot(pie.df8, aes(x = "" , y = value, fill = fct_inorder(group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#f8c4be", "#f18579", "#88d2d6"))+
  geom_label_repel(data = pie.df8b,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 3, nudge_x = 1.75, show.legend = FALSE) +
  guides(fill = "none") +
  theme(legend.position = "none", axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()
pie8

p8p <- p8 + annotation_custom(ggplotGrob(pie8), xmin = -1.8, xmax = 1.2, 
                              ymin = 0, ymax = 1)

print(p8p)

#add legend
leg2add <- get_legend(p1 + theme(legend.position = "bottom", 
                                 legend.title = element_blank(),
                                 legend.text = element_text(size=8))+
                        guides(color = guide_legend(override.aes = list(size = 6, alpha=1))))

pro.all.plot = egg::ggarrange(p1p,p2p,p3p,p4p)
pro.all.plot


euk.all.plot = egg::ggarrange(p5p, p6p, p7p, p8p)
euk.all.plot

title1 <- ggdraw() + draw_label("Prokaryote", hjust = 4, size=12, fontface = "bold")
title2 <- ggdraw() + draw_label("Microbial Eukaryote", hjust = 2.2, size=12, fontface = "bold")
final.p = plot_grid(title1, pro.all.plot, title2, euk.all.plot, leg2add, ncol=1, rel_heights = c(0.1, 1, 0.1, 1.1, 0.1))

pdf("../results/neutral-model-results-with-leg.pdf")
print(final.p)
dev.off()


#Proportions of ASVs above, below and within neutral predictions


snow.df = data.frame(Group = "Prokaryote", Habitat = "Snow", Above=round(pro.snow.plot[[2]][[1]]/sum(pro.snow.plot[[2]])*100,2),
                     Below=round(pro.snow.plot[[2]][[2]]/sum(pro.snow.plot[[2]])*100,2),
                     Neutral=round(pro.snow.plot[[2]][[3]]/sum(pro.snow.plot[[2]])*100,2))
sp.df = data.frame(Group = "Prokaryote", Habitat = "Spring Ice", Above=round(pro.sp.plot[[2]][[1]]/sum(pro.sp.plot[[2]])*100,2),
                   Below=round(pro.sp.plot[[2]][[2]]/sum(pro.sp.plot[[2]])*100,2),
                   Neutral=round(pro.sp.plot[[2]][[3]]/sum(pro.sp.plot[[2]])*100,2))
sum.df = data.frame(Group = "Prokaryote", Habitat = "Summer Ice", Above=round(pro.sum.plot[[2]][[1]]/sum(pro.sum.plot[[2]])*100,2),
                    Below=round(pro.sum.plot[[2]][[2]]/sum(pro.sum.plot[[2]])*100,2),
                    Neutral=round(pro.sum.plot[[2]][[3]]/sum(pro.sum.plot[[2]])*100,2))
cry.df = data.frame(Group = "Prokaryote", Habitat = "Cryoconite", Above=round(pro.cry.plot[[2]][[1]]/sum(pro.cry.plot[[2]])*100,2),
                    Below=round(pro.cry.plot[[2]][[2]]/sum(pro.cry.plot[[2]])*100,2),
                    Neutral=round(pro.cry.plot[[2]][[3]]/sum(pro.cry.plot[[2]])*100,2))

pro.perc.df = rbind(snow.df, sp.df, sum.df, cry.df)

snow.df = data.frame(Group = "Eukaryote", Habitat = "Snow", Above=round(euk.snow.plot[[2]][[1]]/sum(euk.snow.plot[[2]])*100,2),
                     Below=round(euk.snow.plot[[2]][[2]]/sum(euk.snow.plot[[2]])*100,2),
                     Neutral=round(euk.snow.plot[[2]][[3]]/sum(euk.snow.plot[[2]])*100,2))
sp.df = data.frame(Group = "Eukaryote", Habitat = "Spring Ice", Above=round(euk.sp.plot[[2]][[1]]/sum(euk.sp.plot[[2]])*100,2),
                   Below=round(euk.sp.plot[[2]][[2]]/sum(euk.sp.plot[[2]])*100,2),
                   Neutral=round(euk.sp.plot[[2]][[3]]/sum(euk.sp.plot[[2]])*100,2))
sum.df = data.frame(Group = "Eukaryote", Habitat = "Summer Ice", Above=round(euk.sum.plot[[2]][[1]]/sum(euk.sum.plot[[2]])*100,2),
                    Below=round(euk.sum.plot[[2]][[2]]/sum(euk.sum.plot[[2]])*100,2),
                    Neutral=round(euk.sum.plot[[2]][[3]]/sum(euk.sum.plot[[2]])*100,2))
cry.df = data.frame(Group = "Eukaryote", Habitat = "Cryoconite", Above=round(euk.cry.plot[[2]][[1]]/sum(euk.cry.plot[[2]])*100,2),
                    Below=round(euk.cry.plot[[2]][[2]]/sum(euk.cry.plot[[2]])*100,2),
                    Neutral=round(euk.cry.plot[[2]][[3]]/sum(euk.cry.plot[[2]])*100,2))

euk.perc.df = rbind(snow.df, sp.df, sum.df, cry.df)

final.perc.df = rbind(pro.perc.df, euk.perc.df)

#export dataframe
write.csv(final.perc.df, "../results/perc-neutral-model.csv")

####NULL MODEL PLOT BELOW#############
######COMBINE NEUTRAL AND NULL MODELS#########

#load null model results

#bnti
# pro.bnti = read.csv("../results/16S-bNTI-results-table.csv")
# euk.bnti = read.csv("../results/18S-bNTI-results-table.csv")

pro.bnti = read.csv("../results/2024-01-24-null-results/16S-bNTI-results-table.csv")
euk.bnti = read.csv("../results/2024-01-24-null-results/18S-bNTI-results-table.csv")


#bind together
full.bNTI.df = rbind(pro.bnti, euk.bnti)

#export full table
#write.csv(full.bNTI.df, "../results/full-bNTI-results-table.csv")
write.csv(full.bNTI.df, "../results/2024-01-24-null-results/full-bNTI-results-table.csv")

#rcbray results
# pro.rcbray = read.csv("../results/16S-RCbray-results-table.csv")
# euk.rcbray = read.csv("../results/18S-RCbray-results-table.csv")

pro.rcbray = read.csv("../results/2024-01-24-null-results/16S-RCbray-results-table.csv")
euk.rcbray = read.csv("../results/2024-01-24-null-results/18S-RCbray-results-table.csv")

#bind together
full.rcbray.df = rbind(pro.rcbray, euk.rcbray)

#export full table
#write.csv(full.rcbray.df, "../results/full-RCbray-results-table.csv")

write.csv(full.rcbray.df, "../results/2024-01-24-null-results/full-RCbray-results-table.csv")


#GET RESULTS READY FOR PLOTTING

#create dataframes
#make sure everything is formatted correctly for merging
RC_bray <- full.rcbray.df
colnames(RC_bray) <- c("X", "Sample_2", "Sample_1", "RCb", "Group", "Habitat")

turnover.df = inner_join(full.bNTI.df, RC_bray)

turnover.df = turnover.df %>%
  mutate(process = ifelse(abs(bNTI) < 2, #if absolute value of BNTI is less than 2
                          ifelse(abs(RCb) < 0.95, "Drift", #and absolute value of RCb is less than 0.95 = Drift
                                 ifelse(RCb >= 0.95, "Dispersal Limited", #if RCb is > 0.95 = Dispersal limited
                                        ifelse(RCb <= -0.95, "Homogenizing Dispersal", "ERROR"))), #< -0.95 = homo disp
                          ifelse(bNTI >= 2, "Variable Selection", #more than 2 = variable selection
                                 ifelse(bNTI <= -2, "Homogeneous Selection", "ERROR")))) #< -2 = homogenous selection

turnover.df$process = factor(turnover.df$process, levels = c("Drift", 
                                                             "Dispersal Limited", "Homogenizing Dispersal", 
                                                             "Variable Selection", "Homogeneous Selection"))

#for each habitat
#find the number of number of site pairs
#and the percentage of those pairs that shows that process

hab.list <- unique(turnover.df$Habitat)

#Prokaryotes
pro_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Habitat=character())
new.df <- subset(turnover.df, Group == "Prokaryote")

for (i in 1:length(hab.list)){
  x <- new.df[new.df$Habitat == hab.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Habitat<- hab.list[i]
  pro_final_data <- rbind(pro_final_data, x)
}
pro_final_data$Group = "Prokaryote"

#Eukaryotes
euk_final_data <- data.frame(process=character(), n_sites=integer(), perc=numeric(), Habitat=character())
new.df <- subset(turnover.df, Group == "Microbial Eukaryote")

for (i in 1:length(hab.list)){
  x <- new.df[new.df$Habitat == hab.list[i], ]
  x <- x %>%
    group_by(process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/nrow(x))*100) %>%
    as.data.frame
  x$Habitat<- hab.list[i]
  euk_final_data <- rbind(euk_final_data, x)
}
euk_final_data$Group = "Microbial Eukaryote"


#bind results together
final_data = rbind(pro_final_data, euk_final_data) #, mm_final_data)

#add percentages column
final_data.edit <- final_data %>% dplyr::group_by(Habitat, Group) %>%
  dplyr::mutate(frac = n_sites / sum(n_sites))
final_data.edit$frac = round(final_data.edit$frac,2)

#remove rows with zero perc
final_data.edit = final_data.edit[final_data.edit$frac >0,]

#save dataframe
write.csv(final_data.edit, "../results/null-model-results.csv")

#PLOT RESULTS

# final_data.edit = final_data.edit[final_data.edit$Group != "Micrometazoan",]

#make sure everything is in the right order
final_data.edit$Group[final_data.edit$Group == "Eukaryote"] <- "Microbial Eukaryote"
final_data.edit$Group = factor(final_data.edit$Group, 
                               levels = c("Prokaryote", "Microbial Eukaryote"))


final_data.edit$Habitat = factor(final_data.edit$Habitat,
                                 levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))


null.p = ggplot(final_data.edit, aes(x = Habitat,y = n_sites, 
                                 fill = factor(process, levels=c("Homogeneous Selection", "Variable Selection", "Homogenizing Dispersal", "Dispersal Limited", "Drift")))) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels=scales::percent)+
  geom_text(
    aes(label = percent(frac)), position = position_fill(0.5), size=4) +
  facet_grid(~ Group)+
  theme_bw()+
  ylab("Relative Contribution %")+
  theme(text=element_text(size=15), legend.position = "bottom", 
        legend.title = element_blank(), axis.title.x=element_blank(), 
        axis.text=element_text(size=15), legend.text = element_text(size=15), strip.text.x = element_text(size = 15))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_manual(values=c("#f8c4be","#f18579", "#c1e8ea", "#88d2d6","#4fbcc2"))

print(null.p)

pdf("../results/null-model.pdf", height = 10, width = 12)
print(null.p)
dev.off()


####put the plots together

#titles
title1 <- ggdraw() + draw_label("Neutral Model", hjust = 0.95, size=18, fontface = "bold")
title2 <- ggdraw() + draw_label("Null Model", hjust = 1.25, size=18, fontface = "bold")
pro.title <- ggdraw() + draw_label("Prokaryote", hjust = 1.6, size=15, fontface = "bold")
euk.title <- ggdraw() + draw_label("Microbial Eukaryote", hjust = 0.9, size=15, fontface = "bold")

#edit the legend
leg2add <- get_legend(p1 + theme(legend.position = "bottom", 
                                 legend.title = element_blank(),
                                 legend.text = element_text(size=12))+
                        guides(color = guide_legend(override.aes = list(size = 7, alpha=1))))

#edit the axis of the plots
p1p2 = p1p+theme(axis.title.y=element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=15))
p2p2 = p2p+theme(plot.title = element_text(size=15))
p3p2 = p3p+theme(axis.title.y=element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=15))
p4p2 = p4p+theme(plot.title = element_text(size=15))
p5p2 = p5p+theme(axis.title.y=element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=15))
p6p2 = p6p+theme(plot.title = element_text(size=15))
p7p2 = p7p+theme(axis.title=element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=15))
p8p2 = p8p+theme(axis.title.x=element_text(size=10), axis.text.x = element_text(size=10), plot.title = element_text(size=15))

null.p2 = null.p + theme(legend.text = element_text(size=13), axis.text=element_text(size=10), axis.title.y = element_text(size=10))

#put everything together
p_indiv <- wrap_plots(title1+plot_spacer()+plot_spacer()+plot_spacer()+pro.title+plot_spacer()+p1p2+p2p2+p3p2+p4p2+euk.title+plot_spacer()+p5p2+p6p2+p7p2+p8p2+ plot_layout(nrow = 8, heights = c(0.1,0.1,0.1,1,1,0.1,1,1))) 
combi.p <- (p_indiv/leg2add/(title2+plot_spacer())/null.p2) + plot_layout(nrow = 4, heights = c(1.2,0.1,0.05,0.5))
combi.p

pdf("../results/neutral-null-model.pdf", height = 15, width = 8)
print(combi.p)
dev.off()