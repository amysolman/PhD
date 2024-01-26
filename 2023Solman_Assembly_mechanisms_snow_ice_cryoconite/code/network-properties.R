#Script breakdown
#1. Calculate network-level proporties
#2. Random network analysis
#3. Scale-free characteristics analysis

#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(Hmisc)
library(igraph)
library(vegan) #rarecurve function
library(dplyr) #for summarising data 
library(ggpubr)
library(car) #levene test
library(plyr)
library(cowplot)
library(tibble) #rownames to column
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(microViz) #for tax_filter
source("00-solman-functions.R")
library(patchwork)

#Import networks
#positive and negative with negative values to absolute values for calculating topological features
g.sn1 = read.graph("../results/full-snow-network-for-calc-properties.gml", format="gml")
g.sp1 = read.graph("../results/full-springice-network-for-calc-properties.gml", format="gml")
g.sm1 = read.graph("../results/full-summerice-network-for-calc-properties.gml", format="gml")
g.cr1 = read.graph("../results/full-cryoconite-network-for-calc-properties.gml", format="gml")

#positive and negative networks
g.sn2 = read.graph("../results/positive-and-negative-snow-network.gml", format="gml")
g.sp2 = read.graph("../results/positive-and-negative-springice-network.gml", format="gml")
g.sm2 = read.graph("../results/positive-and-negative-summerice-network.gml", format="gml")
g.cr2 = read.graph("../results/positive-and-negative-cryoconite-network.gml", format="gml")

#positive only
g.sn3 = read.graph("../results/positive-snow-network.gml", format="gml")
g.sp3 = read.graph("../results/positive-springice-network.gml", format="gml")
g.sm3 = read.graph("../results/positive-summerice-network.gml", format="gml")
g.cr3 = read.graph("../results/positive-cryoconite-network.gml", format="gml")

#combine for analysis
g.sn <- list(g.sn1, g.sn2, g.sn3)
g.sp <- list(g.sp1, g.sp2, g.sp3)
g.sm <- list(g.sm1, g.sm2, g.sm3)
g.cr <- list(g.cr1, g.cr2, g.cr3)


#1. Calculate network-level properties

g1 = g.sn[[1]]
g2 = g.sn[[3]]
habitat = "Snow"

network_level_properties <- function(g1, g2, habitat){
  
  cor_g <- g1
  cor_p <- g2
  
  #number of nodes (a.k.a. vertices/ASVs)
  gorder(cor_g) 
  
  #number of edges
  gsize(cor_g) 
  
  #number of positive nodes
  gorder(cor_p)
  
  #number of positive edges
  gsize(cor_p)
  
  #mean node degree (number of edges), clustering coefficient (probability that the adjacent vertices of a vertex   are connected), average path length, modularity, density, network diameter 
  #create a dataframe to store out network features
  net_df<-as.data.frame(matrix(NA,ncol=11,nrow=1))
  
  net_df[1,]<-c(habitat, round(mean(degree(cor_g)),3), 
                round(transitivity(cor_g),3), 
                round(average.path.length(cor_g), 3),
                round(graph.density(cor_g),3),
                round(diameter(cor_g),3),
                round(modularity(walktrap.community(cor_g)),3), 
                gorder(cor_g), gorder(cor_g)-gorder(cor_p), gsize(cor_g), gsize(cor_g)-gsize(cor_p))
  
  colnames(net_df)<-c("Habitat", "AveDegree","ClustCoef","AvePathLen","Density","Diameter","Modularity", "Nodes", "Negative_Nodes", "Edges", "Negative_Edges")
  
  
  return(net_df)
}

sn.net <- network_level_properties(g.sn[[1]], g.sn[[3]], "Snow")
sp.net <- network_level_properties(g.sp[[1]], g.sp[[3]], "Spring Ice")
sm.net <- network_level_properties(g.sm[[1]], g.sm[[3]], "Summer Ice")
cr.net <- network_level_properties(g.cr[[1]], g.cr[[3]], "Cryoconite")

#bind together
t.net = rbind(sn.net, sp.net, sm.net, cr.net)


#2. Random network analysis

# g = g.sp[[1]]
# net_df = sp.net
# habitat = "Spring Ice"

random_network <- function(g, net_df, habitat){
  
  cor_g <- g
  
  #1) generate 1000 random networks with the same number of nodes (asvs) and edges (connections) as our real network 
  set.seed(1000)
  gs <- list()
  for (x in 1:1000) {
    gs[[x]] <- erdos.renyi.game(gorder(cor_g), gsize(cor_g), type = "gnm", directed = FALSE,loops = FALSE)
  }
  
  #2) create a dataframe and store our random network features
  net_df2<-as.data.frame(matrix(NA,ncol=6,nrow=1000))
  
  colnames(net_df2)<-c("AveDegree","ClustCoef","AvePathLen","Density","Diameter","Modularity")
  
  for (x in 1:1000){
    net_df2[x,]<-c(mean(degree(gs[[x]])),transitivity(gs[[x]]),average.path.length(gs[[x]]),graph.density(gs[[x]]),diameter(gs[[x]]),modularity(walktrap.community(gs[[x]])))
  }
  
  #3) Table the real network properties with the mean random network properties and their sd
  res = data.frame(Modularity= c(net_df$Modularity, paste0(round(mean(net_df2$Modularity), 3), " (\u00B1", round(sd(net_df2$Modularity), 3),")") ),
                   
                   ClustCoeff= c(net_df$ClustCoef, paste0(round(mean(net_df2$ClustCoef), 3), " (\u00B1", round(sd(net_df2$ClustCoef), 3),")") ),
                   
                   AvePathLen= c(net_df$AvePathLen, paste0(round(mean(net_df2$AvePathLen), 3), " (\u00B1", round(sd(net_df2$AvePathLen), 3),")") )
                   
                   
  )
  
  
  rownames(res) = c("Real Network", "Random Network Mean (SD)")
  
  #3) Is the mean score of properties for the random network significantly lower than real network properties?
  
  res.1 = wilcox.test(net_df2$Modularity, mu=as.numeric(net_df$Modularity), alternative = "less")
  res.1
  res.2 = wilcox.test(net_df2$ClustCoef, mu=as.numeric(net_df$ClustCoef), alternative = "less")
  res.2
  res.3 = wilcox.test(net_df2$AvePathLen, mu=as.numeric(net_df$AvePathLen), alternative = "less")
  res.3
  
  df = data.frame(Habitat = habitat, ModularityPVal = res.1$p.value, ClustCoefPVal = res.2$p.value, AvePathLenPVal = res.3$p.value)
  
  res.list = list(res, df)
  
  return(res.list)
}

sn.rand = random_network(g.sn[[1]], sn.net, "Snow")
sp.rand = random_network(g.sp[[1]], sp.net, "Spring Ice")
sm.rand = random_network(g.sm[[1]], sm.net, "Summer Ice")
cr.rand = random_network(g.cr[[1]], cr.net, "Cryoconite")

w.res = rbind(sn.rand[[2]], sp.rand[[2]], sm.rand[[2]], cr.rand[[2]])

#Merge data frames to make table

#Snow
r = data.frame(AveDegree = NA, ClustCoef = sn.rand[[1]][2,2], AvePathLen = sn.rand[[1]][2,3], Density = NA, Diameter = NA, Modularity = sn.rand[[1]][2,1], Nodes = sn.net$Nodes[1], Negative_Nodes = NA, Edges = sn.net$Edges[1], Negative_Edges = NA)
sn.tab = rbind(sn.net[,-c(1)], r)
rownames(sn.tab) <- c("Full", "Random Network Mean (SD)")
df1 <- tibble::rownames_to_column(sn.tab, "Network")
df1 = cbind(rep("Snow", 2), df1)
names(df1)[1] = "Habitat"

#Spring Ice
r = data.frame(AveDegree = NA, ClustCoef = sp.rand[[1]][2,2], AvePathLen = sp.rand[[1]][2,3], Density = NA, Diameter = NA, Modularity = sp.rand[[1]][2,1], Nodes = sp.net$Nodes[1], Negative_Nodes = NA, Edges = sp.net$Edges[1], Negative_Edges = NA)
sp.tab = rbind(sp.net[,-c(1)], r)
rownames(sp.tab) <- c("Full", "Random Network Mean (SD)")
df2 <- tibble::rownames_to_column(sp.tab, "Network")
df2 = cbind(rep("Spring Ice", 2), df2)
names(df2)[1] = "Habitat"

#Summer Ice
r = data.frame(AveDegree = NA, ClustCoef = sm.rand[[1]][2,2], AvePathLen = sm.rand[[1]][2,3], Density = NA, Diameter = NA, Modularity = sm.rand[[1]][2,1], Nodes = sm.net$Nodes[1], Negative_Nodes = NA, Edges = sm.net$Edges[1], Negative_Edges = NA)
sm.tab = rbind(sm.net[,-c(1)], r)
rownames(sm.tab) <- c("Full", "Random Network Mean (SD)")
df3 <- tibble::rownames_to_column(sm.tab, "Network")
df3 = cbind(rep("Summer Ice", 2), df3)
names(df3)[1] = "Habitat"

#Cryoconite
r = data.frame(AveDegree = NA, ClustCoef = cr.rand[[1]][2,2], AvePathLen = cr.rand[[1]][2,3], Density = NA, Diameter = NA, Modularity = cr.rand[[1]][2,1], Nodes = cr.net$Nodes[1], Negative_Nodes = NA, Edges = cr.net$Edges[1], Negative_Edges = NA)
cr.tab = rbind(cr.net[,-c(1)], r)
rownames(cr.tab) <- c("Full", "Random Network Mean (SD)")
df4 <- tibble::rownames_to_column(cr.tab, "Network")
df4 = cbind(rep("Cryoconite", 2), df4)
names(df4)[1] = "Habitat"


final.net.df = rbind(df1, df2, df3, df4)

#add edge/node ratio
final.net.df$EdgeNodeRatio = round(as.numeric(final.net.df$Edges)/as.numeric(final.net.df$Nodes), 3)

#export dataframe
write.csv(final.net.df, "../results/network-properties.csv")

#Are differences between the random networks and real networks significant?
rand_real_sig_dif = rbind(sn.rand[[2]], sp.rand[[2]], sm.rand[[2]], cr.rand[[2]])
write.csv(rand_real_sig_dif, "../results/wilcoxon-rand-real-networks.csv")


#3. Scale-free characteristics analysis

# Function for finding if our network exhibits scale-free characteristics
# Code taken from http://chengjun.github.io/web_data_analysis/demo2_simulate_networks/

scale_free_characteristics <- function(g, habitat){
  
  cor_g <- g
  
  # 1) power law distribution function
  fit_power_law = function(g, network) {
    
    set.seed(666)
    # calculate degree
    d = degree(g, mode = "all")
    dd = degree.distribution(g, mode = "all", cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    # delete blank values
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    
    res.list = list(probability, degree, cozf, d, alpha, R.square)
    
    return(res.list)
  }
  
  
  #2) generate scale-free network
  g.big.ba = barabasi.game(gorder(cor_g))
  
  #3) generate random network
  g.big.er = erdos.renyi.game(gorder(cor_g), 0.1)
  
  #4) Get the info for plotting
  #true network
  real.res = fit_power_law(cor_g, "real")
  power.law.fit1 = function(x) exp(real.res[[3]][[1]] + real.res[[3]][[2]] * log(x))
  #scale-free network
  scale.res = fit_power_law(g.big.ba, "scale-free")
  power.law.fit2 = function(x) exp(scale.res[[3]][[1]] + scale.res[[3]][[2]] * log(x))
  #random network
  random.res = fit_power_law(g.big.er, "random")
  power.law.fit3 = function(x) exp(random.res[[3]][[1]] + random.res[[3]][[2]] * log(x))
  
  #plot
  pdf(paste0("../results/", habitat, "power-law-model-fits.pdf"), width=10, height=6)
  par(mfrow=c(1,3), mai = c(0.6, 0.6, 0.4, 0.1))
  plot(real.res[[1]] ~ real.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "A", cex.lab=1.5, cex.axis=1.5, cex.main=3)
  curve(power.law.fit1, col = "red", add = T, n = length(real.res[[4]]))
  text(cex=1.6, x=min(real.res[[2]]), y=max(real.res[[1]]), adj = 0, labels=paste0("R2 = ", round(real.res[[6]], 2), "  alpha = ", round(real.res[[5]], 2)))
  plot(scale.res[[1]] ~ scale.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "B", cex.lab=1.5, cex.axis=1.5, cex.main=3)
  curve(power.law.fit2, col = "red", add = T, n = length(scale.res[[4]]))
  text(cex=1.6, x=min(scale.res[[2]]), y=max(scale.res[[1]]), adj = 0,labels=paste0("R2 = ", round(scale.res[[6]], 2), "  alpha = ", round(scale.res[[5]], 2)))
  plot(random.res[[1]] ~ random.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "C", cex.lab=1.5, cex.axis=1.5, cex.main=3)
  curve(power.law.fit3, col = "red", add = T, n = length(random.res[[4]]))
  text(cex=1.6, x=min(random.res[[2]]), y=max(random.res[[1]]), adj = 0,labels=paste0("R2 = ", round(random.res[[6]], 2), "  alpha = ", round(random.res[[5]], 2)))
  dev.off()

  return(res.list=list(real.res[[6]], scale.res[[6]], random.res[[6]]))
}



sn.plm = scale_free_characteristics(g.sn[[1]], "Snow")
sp.plm = scale_free_characteristics(g.sp[[1]], "SpringIce")
sm.plm = scale_free_characteristics(g.sm[[1]], "SummerIce")
cr.plm = scale_free_characteristics(g.cr[[1]], "Cryoconite")