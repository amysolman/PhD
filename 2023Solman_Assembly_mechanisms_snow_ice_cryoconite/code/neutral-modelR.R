## Method

# Code developed from https://github.com/mrebolleda/OrganFilters_MimulusMicrobiome/blob/master/Neutral_pol_nopol.md and Burns et al., 2016
# 
# To determine the potential contribution of stochastic processes on the assemblage of abundant, intermediate and rare taxa, Sloan’s neutral community model (Sloan et al., 2006, 2007) was applied to predict the association between ASV observation frequency (i.e. the number of samples each ASV is observed in) and their relative abundance in the metacommunity. Code developed from Gómez & Ashman (2019) (https://github.com/mrebolleda/OrganFilters_MimulusMicrobiome/blob/master/Neutral_pol_nopol.md) was used for fitting the model. In a system shaped by stochastic processes of drift and dispersal, the model predicts taxa with greater relative abundance in the metacommunity will be observed in a greater number of sites, as a large population increases the likelihood of random dispersal between sites. Conversely, rare taxa are less likely to disperse between sites and are more prone to local elimination through ecological drift. With a single free parameter m (migration rate), the model predicts the frequency with which an ASV will be observed, considering it's mean relative abundance in the metacommunity. Akaike informaition criterion (AIC) was used to compare the fit of the neutral model to a binomial distribution, serving as a null model of random sampling from the metacommunity without drift or dispersal limitation. 95% confidence intervals were calculated around model predictions. ASVs within the 95% CI were considered to be distributed according to neutral expectations.
# 
# To determine the potential contribution of stochastic processes (birth, death and migration) to the assemblage of snow, spring ice, summer ice and cryoconite communities, Sloan’s neutral community model was applied (Sloan et al., 2006). The model asserts that, under neutral assumptions, taxa that are more abundant in the metacommunity will be observed in a greater number of areas (samples). Using a single free parameter m (migration rate), Sloan’s community model predicts the association between taxa relative abundance in the metacommunity (i.e. all samples) and their observation frequency (i.e. the number of samples in which they are observed). In brief, for each ASV a metacommunity abundance value was calculated (the relative abundance of that ASV within the entire dataset), as well as a frequency value (the percentage of samples that ASV was observed in). Non-linear least squares (NLLS) was used to fit the model to the data and predict m. 
# •	The model takes a vector of observation frequencies for each ASV and uses the limit of detection, the mean number of reads per sample, and a vector of mean relative abundances of each ASV to fit the model parameter m to best predict the vector of observation frequencies. 
# •	Confint function of the stats package was used to compute confidence intervals for m
# •	The free parameter m was also fit to the model using maximum likelihood estimation
# •	AIC and BIC values were calculated for the NLLS and MLE model fits
# •	The frequency predictions of the model were calculated using the estimated value of m.
# •	R2 and root mean squares errors (RMSE) was calculated for the NLLS model fit.
# •	Frequency prediction confidence intervals were calculated using binconf function.
# •	Binomial model fit to the data and AIC and BIC calculated.
# •	R2 and RMSE calculated for the binomial model.
# •	Poisson model fit to the data and AIC and BIC calculated.
# •	R2 and RMSE calculated for the poisson model.
# 
# The fit of the neutral model was compared to that of a binomial distribution null model representative of the random distribution of taxa without the impact of ecological drift (random birth/death events) and migration. The models were compared using Akaike information criterion. 95% confidence intervals were calculated around the model frequency predictions for each ASV Code developed from Rebolleda Gómez and Ashman (2019) was used for fitting the model (https://github.com/mrebolleda/OrganFilters_MimulusMicrobiome/blob/master/Neutral_pol_nopol.md). 
# 
# 
# 
# In the neutral model, we calculate the mean relative abundance of each ASV by:
# 1. Calculating the mean number of reads per sample
# 2. Calculating the mean number of reads in each sample for each ASV (total number of reads/total number of samples)
# 3. Dividing the mean number of reads for each ASV by the mean number of reads per sample
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

# 1. Clear workspace and load packages
rm(list=ls())
graphics.off()

source("00-solman-functions.R")
library(phyloseq)
# library(RColorBrewer)
# library(tidyr) #trans data between wide and long format
 library(ggplot2)
# library(vegan) #for vegdist
# source("00-solman-functions.R")
# library(BiodiversityR) #for rank abundance curves
 library(dplyr) #for %>% function
 library(cowplot) #plot_grid
# library(reshape2) #for function melt
# library(funrar) #for make relative
# library(stringr) #mutate function, to split columns
# library(gridExtra) #for exporting as pdf
# library(scales)
# library(microbiome) #for summarize_phyloseq
# library(picante) #for faith's pd
# library(ggpmisc) #for adding lm stats to the plot
# library(broom) #for adding stats to the plot
# library(fossil) #for earth.dist function
library(Hmisc) #binconf function
library(minpack.lm) #for non-linear model fitting
library(tibble) #function rownames_to_column
library(stats4) #mle function
library(patchwork)
library(egg)
library(ggrepel)
library(tidyverse)

#prokaryotes
ps.pro <- readRDS("../results-proportional/16S-ps-controls-removed.rds") 
#eukaryotes without micrometazoans
ps.euk <- readRDS("../results-proportional/18S-ps-no-mm-controls-removed.rds") 
#micrometazoans only
ps.mm <- readRDS("../results-proportional/18S-ps-mm-only-controls-removed.rds") 

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


ps = ps.pro
hab = data.frame(sample_data(subset_samples(ps.pro, Habitat == "Snow")))$SampleID

#neutral_model <- function(ps, hab){

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

#plot our data
plot(p, freq)

#Make N an integer for binomial model fitting
N.in = as.integer(N)

#1. FIT THE MODEL TO THE DATA

#SLOAN'S BETA MODEL

#the model takes two vectors of values, that means a p parameter and q parameter is given for each ASV
#using these parameters for each ASV the model predicts the proportion of sites it will be observed in
#these parameters are varied (by changing m) until the difference between model predictions 
#and true frequencies is minimised this results in a prediction for m

#Fit model parameter m (immigration rate) using Non-linear least squares (NLS)
#increase in m = increase in both alpha and beta shapes for the beta distribution function
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001)) 

#using the estimated parameter for m get the predicted ASV frequencies under stochastic processes
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE) 

#BINOMIAL MODEL

#1.A Using the binomial distribution function, 
#the size/number of trials (N) and the probability/mean relative abundance (p)
#to estimate the proportion of trials/sites in which each ASV will be found
bino.pred <- pbinom(d, N.in, p, lower.tail=FALSE)

#1.B A binomial model is created by... 
#a) finding predicted frequencies using the binomial distribution function
#b) finding the difference between true and predicted values
#c) Use maximum likelihood estimation with those differences as quantiles
# and varying mean and sigma to reduce the difference between true and predicted values
# to fit a binomial model to the data

#this will minimise the difference between real frequency and 
#predicted frequency by varying values of mean and standard deviation.
bino.LL <- function(mu, sigma){ #the mean (mu) and standard deviation (sigma)
  
  #using the number of trials/size/population size of sample (N) and 
  #the probability of success/relative abundance
  #of each ASV (p) give the proportion of successes = frequency of observation
  bi = pbinom(d, N.in, p, lower.tail=FALSE) 
  
  #what is the difference between the true frequency and the predictions?
  R = freq - b
  
  #reduce differences by changing mean and standard deviation
  R = dnorm(R, mu, sigma) 
  -sum(log(R))
  
}

#maximum likelihood estimation model fitting by estimating values of 
#mean (mu) and standard deviation (sigma) 
bino.mle <- mle(bino.LL, start=list(mu=0.1, sigma=0.1), nobs=length(p))


#POSSION MODEL

#Get model predictions
pois.pred <- ppois(d, N*p, lower.tail=FALSE)

#Define the model function
pois.LL <- function(mu, sigma){
  R = freq - ppois(d, N*p, lower.tail=FALSE)
  R = dnorm(R, mu, sigma)
  -sum(log(R))
}

#fit using MLE
pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))

#2. QUANTIFYING MODEL ACCURACY/PERFORMANCE

#2.A plot true data against predictions
plot(freq.pred, freq)
plot(bino.pred, freq)
plot(pois.pred, freq)

#2.B Using the true frequencies and the predicted frequencies, calculate the R2 value = the proportion of variance in the
#the dependent variable (true frequency) explained by the independent variable (predicted frequency)
#R2 is calculated as the residual sum of squares / total sum of squares
#Residual sum of squares = variation between the real data and predictions
#Total sum of squares = variation within the real data (data points from the mean)

#SLOANS BETA 
RSS <- sum(residuals(m.fit)^2) #Residual sum of squares of our NLLS model
TSS <- sum((freq - mean(freq))^2) #Total sum of squares 
Rsqr <- 1 - (RSS/TSS) #R-squared value
#BINOMIAL
Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
#POISSON 
Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))


#2.C Using the true frequencies and the predicted frequencies calculate the root mean squared error (RMSE)
# = square root of the residual sum of squares divided by the number of datapoints
# = a measure of variation not captured by the model

#SLOANS BETA 
RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
#BINOMIAL
RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
#POISSON
RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

#2.D Calculate AIC/BIC as measures of how well the model fit the data, 
#considering the number of parameters in the model
#lower values are better

#Beta model
aic.fit = AIC(m.fit)
bic.fit = BIC(m.fit)

#Binomial model
aic.bino <- AIC(bino.mle, k=2)
bic.bino <- BIC(bino.mle)

#POISSON MODEL
aic.pois <- AIC(pois.mle, k=2)
bic.pois <- BIC(pois.mle)

#3. QUANTIFYING MODEL UNCERTAINTY

#BETA model
#get the confidence intervals for m
m.ci <- confint(m.fit, 'm', level=0.95)

# Get confidence interval for predictions
freq.pred.ci <- binconf(freq.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)

#Binomial model
#calculate 95% confidence intervals for the predicted values 
bino.pred.ci <- binconf(bino.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)

#POISSON MODEL
#calculate 95% confidence intervals for the predicted values 
pois.pred.ci <- binconf(pois.pred*nrow(ASV.table), nrow(ASV.table), alpha=0.05, method="wilson", return.df=TRUE)


```

```{r}
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

#       return(res.list)
# }

```

```{r, include=FALSE}
#fit the models
pro.snow <- neutral_model(ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Snow")))$SampleID)
pro.sp <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Spring Ice")))$SampleID)
pro.sum <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Summer Ice")))$SampleID)
pro.cry <- neutral_model(ps=ps.f.pro, hab = data.frame(sample_data(subset_samples(ps.f.pro, Habitat == "Cryoconite")))$SampleID)

euk.snow <- neutral_model(ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Snow")))$SampleID)
euk.sp <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Spring Ice")))$SampleID)
euk.sum <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Summer Ice")))$SampleID)
euk.cry <- neutral_model(ps=ps.f.euk, hab = data.frame(sample_data(subset_samples(ps.f.euk, Habitat == "Cryoconite")))$SampleID)

# mm.snow <- neutral_model(ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Snow")))$SampleID)
# mm.sp <- neutral_model(ps=ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Spring Ice")))$SampleID)
# mm.sum <- neutral_model(ps=ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Summer Ice")))$SampleID)
# mm.cry <- neutral_model(ps=ps.mm, hab = data.frame(sample_data(subset_samples(ps.mm, Habitat == "Cryoconite")))$SampleID)

```