#Script breakdown
#1. Clear workspace and load packages
#2. Load data
#3. Get abundance dataframes
#4. Merge 16S and 18S datasets
#5. Filter ASVs by percentage of samples
#6. Carry out correlation analysis
#7. Calculate node-level properties
#8. Calculate network-level proporties
#9. Random network analysis
#10. Scale-free characteristics analysis
#11. Get taxonomy of nodes + report
#12. Find connections between subcommunities
#13. Find differences between subcommunities
#14. Identify keystone taxa

# Method 

#16S and 18S rRNA datasets were merged prior to network analysis. A co-occurrence network was constructed based on Spearman’s rank correlation between ASVs using rcorr function of Hmisc package (Harrell and Dupont, 2022). The analysis was carried out using ASVs present in at least 20% of samples to reduce false high correlations (Berry and Widder, 2014). Strong positive correlations (ρ > 0.7, fdr-adjusted p value < 0.01) between ASVs were exported using igraph package (Csardi and Nepusz, 2006). Node (ASV) level properties including degree, betweenness centrality, eigenvector centrality and closeness centrality were calculated. Subgraphs were generated for abundant and rare subcommunities to calculate total and subcommunity network level properties including average degree, clustering coefficient, average path length (APL), density, diameter, modularity, number of nodes and number of edges. The number of connections between subcommunities was calculated. The stat_compare_means function from the ggpubr package (Kassambara and Kassambara, 2020) was used to compare mean node properties between abundant and rare subcommunities. 1000 Erdös–Réyni random networks were constructed with the same number of nodes and edges as the real network. The mean degree, clustering coefficient, APL, density, diameter and modularity were calculated for each random network. Structural properties of the real and random network were compared to determine the degree of randomness in the real network. Most real-world networks are scale-free. A scale-free network is a network whose degrees follow a power-law distribution, that is, a small fraction of nodes has many connections and a large number of nodes have a small number of connections. In a random network there is an equitable distribution of connections between nodes, that is, most nodes have an average number of degrees (i.e. Poisson distribution). To test if the real network exhibited scale-free characteristics the power law distribution was fit to the real network, a random network and a scale-free network. Network visualisation and modular analysis were carried out using the interactive platform Gephi (Bastian et al., 2009) and the ggtern package in R. The Louvain algorithm used to characterise modules. Keystone taxa were identified as those with standardized high degree abundance, high closeness centrality and low betweenness centrality (Berry and Widder, 2014; Banerjee et al., 2018).  


# Results

#1. Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
library(Hmisc)
library(igraph)
#library(agricolae)
#library(fdrtool)
library(vegan) #rarecurve function
library(dplyr) #for summarising data 
library(ggpubr)
library(car) #levene test
library(plyr)
library(cowplot)
library(tibble) #rownames to column
source("00-solman-functions.R")

#2. Load data
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

#3. Get abundance dataframes

#Abundance df function
abun_dataframe_func <- function(phylo.abun, phylo.int, phylo.rare){
  
  #Get data frame with classifications
  AT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.abun))), Abundance = rep("AT", length(rownames(data.frame(otu_table(phylo.abun))))))
  MT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.int))), Abundance = rep("MT", length(rownames(data.frame(otu_table(phylo.int))))))
  RT_DATA = data.frame(ID=rownames(data.frame(otu_table(phylo.rare))), Abundance = rep("RT", length(rownames(data.frame(otu_table(phylo.rare))))))
  
  
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


# phylo1 = pro.ant
# phylo2 = euk.ant
#4. Merge 16S and 18S datasets
ant.ps = merge_my_phylo_same_pole(pro.ant, euk.ant)
arc.ps = merge_my_phylo_same_pole(pro.arc, euk.arc)

#5. Filter ASVs by percentage of samples

filter_me <- function(my_df, perc){
  
  #remove ASVs present in less than X% samples
  df <- 1*(my_df>0) #presence/abundance df
  pres_vec <- vector() 
  asv_vec <- vector()
  for (k in 1:ncol(df)){  #for each ASV
    pres <- as.numeric(colSums(df)[k])/nrow(df) #find the percentage of samples the ASV is present in
    pres_vec <- c(pres_vec, pres)
    asv_vec <- c(asv_vec, colnames(df)[k])
  }
  df_2 <- data.frame(ASV=asv_vec, Present=pres_vec) #df of percentage of samples asvs are present in
  
  trim_me <- subset(df_2, df_2$Present >= perc) #only keep ASVs in >= X% of samples
  
  #use ASV IDs to subset count table
  asv_pres_trim <- my_df[, (colnames(my_df) %in% trim_me$ASV)]
  
  return(asv_pres_trim)
}  


# my_df = data.frame(t(otu_table(arc.ps)), check.names=FALSE)
# num = 6

filter_me2 <- function(my_df, num){
  
  #remove ASVs present in less than X% samples
  df <- 1*(my_df>0) #presence/abundance df
  pres_vec <- vector() 
  asv_vec <- vector()
  for (k in 1:ncol(df)){  #for each ASV
    num_pres <- as.numeric(colSums(df)[k]) #find the num of samples the ASV is present in
    pres_vec <- c(pres_vec, num_pres)
    asv_vec <- c(asv_vec, colnames(df)[k])
  }
  df_2 <- data.frame(ASV=asv_vec, Present=pres_vec) #df of num of samples asvs are present in
  
  trim_me <- subset(df_2, df_2$Present >= num) #only keep ASVs in > num of samples
  
  #use ASV IDs to subset count table
  asv_pres_trim <- my_df[, (colnames(my_df) %in% trim_me$ASV)]
  
  return(asv_pres_trim)
}  


#6. Carry out correlation analysis

get_my_network <- function(phylo, num, pole){
  
  set.seed(666) 
  
  #get our asv table
  asv_tab <- data.frame(t(otu_table(phylo)), check.names=FALSE)
  asv.tab.filter <- filter_me2(asv_tab, num)
  
  #find the spearmans correlations between asvs
  cor_analysis <-rcorr(as.matrix(asv.tab.filter),type="spearman")
  
  #assign the matrix of correlations and the p-values to r and p, respectively
  cor_r<-cor_analysis$r
  cor_p<-cor_analysis$P 
  
  #matrix diagonals - apply the value 1 to the diagonal values of the p matrix (instead of NA as they were before)
  diag(cor_p)<-1
  
  #make p-values matrix into vector
  cor_pp<-as.vector(cor_p)
  
  #how many of our p values are < 0.01
  length(cor_pp[cor_pp > 0 & cor_pp < 0.01]) 
  
  #adjust p-values using the false discovery rate method
  cor_p.adj <- p.adjust(cor_p, method="fdr")
  
  #make adjusted p-values matrix into vector
  cor_pp.adj<-as.vector(cor_p.adj)
  
  #how many of our adjusted p values are < 0.01
  length(cor_pp.adj[cor_pp.adj > 0 & cor_pp.adj < 0.01]) 
  
  #make any value of p > 0.01 a 0
  cor_p.adj[cor_p.adj>0.01]<-0
  
  #make any value <= 0.01 & greater than 0 = 1
  cor_p.adj[cor_p.adj<=0.01&cor_p.adj>0]<-1
  
  
  #let's look at the range of correlation coefficients
  hist(cor_r)
  
  #are there any with coefficients > | 0.8|?
  length(cor_r[cor_r < -0.8]) #there are  negative correlation coefficients
  length(cor_r[cor_r > 0.8]) #there  positive correlation coefficients
  
  #change the value of any correlation coefficients <  0.8  to 0
  cor_r[(cor_r)<0.8]<-0
  
  #multiple the correlation coefficients by adjusted p values so that insignificant coefficients become 0
  cor<-cor_r*cor_p.adj
  
  #are there any with coefficients > | 0.8| after removing non-significant correlations?
  length(cor[cor > 0.8]) 
  
  #create igraph graph from adjacency matrix
  cor_g <- graph.adjacency(cor, weighted=TRUE, mode="undirected")
  
  #simplify the graph so it does not contain loops or multiple edges
  cor_g <- simplify(cor_g)
  
  #delete vertices from the graph
  cor_g<-delete.vertices(cor_g,names(degree(cor_g)[degree(cor_g)==0]))
  
  #export graph
  write.graph(cor_g, paste0("../results/", pole, "-network.gml"), format="gml")
  
  return(cor_g)
}

g.ant = get_my_network(ant.ps, 6, "Antarctic")
g.arc = get_my_network(arc.ps, 6, "Arctic")


#7. Calculate node-level properties

node_level_properties <- function(graph, abun_df){
  
  cor_g <- graph
  
  rare <- subset(abun_df, (abun_df$Abundance == "RT"))
  intermediate <- subset(abun_df, (abun_df$Abundance == "MT"))
  abundant <- subset(abun_df, (abun_df$Abundance == "AT"))
  
  #create a dataframe with 5 columns and a row for each node (ASV)
  df<-as.data.frame(matrix(NA,ncol=5,nrow=length(degree(cor_g))))
  
  #make ASV IDs row names of df
  rownames(df)<-names(degree(cor_g))
  
  #name the colums the node-level topological features
  colnames(df)<-c("degree","betweenness","closeness","eigenvector","category")
  
  #categorise the ASVs as rare, abundant or intermediate
  df[intersect(names(degree(cor_g)),rare$ID),5]<-"rare"
  df[intersect(names(degree(cor_g)),abundant$ID),5]<-"abundant"
  df[intersect(names(degree(cor_g)),intermediate$ID),5]<-"intermediate" 
  
  #get betweenness
  btw<-betweenness(cor_g)
  #get closeness centrality
  cls<-closeness(cor_g)
  #get eigenvector centrality
  egv<-evcent(cor_g)
  
  #put the topological features into the dataframe
  df[,1]<-degree(cor_g)
  df[,2]<-btw
  df[,3]<-cls
  df[,4]<-egv$vector
  
  return(df)
  
}

node.prop.ant = node_level_properties(g.ant, ant.df)
node.prop.arc = node_level_properties(g.arc, arc.df)

node.prop.ant %>% 
  dplyr::group_by(category) %>% 
  dplyr::summarise(across(everything(), list(mean)))

node.prop.arc %>% 
  dplyr::group_by(category) %>% 
  dplyr::summarise(across(everything(), list(mean)))


#8. Calculate network-level properties

# graph = g.ant
# node_df = node.prop.ant

network_level_properties <- function(graph, node_df){
  
  cor_g <- graph
  df <- node_df
  
  #subset the node level dataframe by rare ASVs
  a<-subset(df,category=="rare")
  g_ra<-induced_subgraph(cor_g,rownames(a))
  #abundant
  a<-subset(df,category=="abundant")
  g_abd<-induced_subgraph(cor_g,rownames(a))
  #intermediate
  a<-subset(df,category=="intermediate")
  g_int<-induced_subgraph(cor_g,rownames(a))
  
  #number of nodes (a.k.a. vertices/ASVs)
  gorder(cor_g) 
  gorder(g_ra) 
  gorder(g_abd) 
  gorder(g_int) 
  
  #number of edges
  gsize(cor_g) 
  gsize(g_ra) 
  gsize(g_abd) 
  gsize(g_int)
  
  #mean node degree (number of edges), clustering coefficient (probability that the adjacent vertices of a vertex   are connected), average path length, modularity, density, network diameter 
  #create a dataframe to store out network features
  net_df<-as.data.frame(matrix(NA,ncol=8,nrow=4))
  
  net_df[1,]<-c(mean(degree(cor_g)),transitivity(cor_g),average.path.length(cor_g),
                graph.density(cor_g),diameter(cor_g),modularity(walktrap.community(cor_g)), gorder(cor_g), gsize(cor_g))
  
  net_df[2,]<-c(mean(degree(g_abd)),transitivity(g_abd),average.path.length(g_abd),
                graph.density(g_abd),diameter(g_abd),modularity(walktrap.community(g_abd)), gorder(g_abd), gsize(g_abd))
  
  net_df[3,]<-c(mean(degree(g_ra)),transitivity(g_ra),average.path.length(g_ra),
                graph.density(g_ra),diameter(g_ra),modularity(walktrap.community(g_ra)), gorder(g_ra), gsize(g_ra))
  
  net_df[4,]<-c(mean(degree(g_int)),transitivity(g_int),average.path.length(g_int),
                graph.density(g_int),diameter(g_int),modularity(walktrap.community(g_int)), gorder(g_int), gsize(g_int))
  
  colnames(net_df)<-c("AveDegree","ClustCoef","AvePathLen","Density","Diameter","Modularity", "Nodes", "Edges")
  rownames(net_df)<-c("Total", "Abundant","Rare","Intermediate")
  
  new.net = net_df %>% 
    mutate_if(is.numeric, round, 3)
  
  return(new.net)
}

net.df.ant = network_level_properties(g.ant, node.prop.ant)
net.df.arc = network_level_properties(g.arc, node.prop.arc)

#9. Random network analysis

random_network <- function(graph, net_df){
  
  cor_g <- graph
  net_df <- net_df
  
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
  res = data.frame(Modularity= c(round(net_df$Modularity[1], 3), paste0(round(mean(net_df2$Modularity), 3), " (\u00B1", round(sd(net_df2$Modularity), 3),")") ),
                   
                   ClustCoeff= c(round(net_df$ClustCoef[1], 3), paste0(round(mean(net_df2$ClustCoef), 3), " (\u00B1", round(sd(net_df2$ClustCoef), 4),")") ),
                   
                   AvePathLen= c(round(net_df$AvePathLen[1], 3), paste0(round(mean(net_df2$AvePathLen), 3), " (\u00B1", round(sd(net_df2$AvePathLen), 5),")") )
                   
                   
  )
  
  
  rownames(res) = c("Real Network", "Random Network Mean (SD)")
  
  #3) Is the mean score of properties for the random network significantly lower than real network properties?
  
  #First, make sure the data are normally distributed using the Shapiro-Wilk test then perform either a t-test or wilcox test
  norm.test1 = shapiro.test(net_df2$Modularity) #p > 0.05 so data are normally distributed
  if (norm.test1$p.value > 0.05){
    res.1 = t.test(net_df2$Modularity, mu=net_df$Modularity[1], alternative = "less")
  } else {
    res.1 = wilcox.test(net_df2$Modularity, mu=net_df$Modularity[1], alternative = "less")
    
  }
  
  norm.test2 = shapiro.test(net_df2$ClustCoef) #p > 0.05 so data are normally distributed
  if (norm.test2$p.value > 0.05){
    res.2 = t.test(net_df2$ClustCoef, mu=net_df$ClustCoef[1], alternative = "less")
  } else {
    res.2 = wilcox.test(net_df2$ClustCoef, mu=net_df$ClustCoef[1], alternative = "less")
    
  }
  
  norm.test3 = shapiro.test(net_df2$ClustCoef) #p > 0.05 so data are normally distributed
  if (norm.test3$p.value > 0.05){
    res.3 = t.test(net_df2$AvePathLen, mu=net_df$AvePathLen[1], alternative = "less")
  } else {
    res.3 = wilcox.test(net_df2$AvePathLen, mu=net_df$AvePathLen[1], alternative = "less")
    
  }
  
  
  df = data.frame(ModularityPVal = res.1$p.value, ClustCoefPVal = res.2$p.value, AvePathLenPVal = res.3$p.value)
  
  res.list = list(res, df)
  
  return(res.list)
}

rand.ant = random_network(g.ant, net.df.ant)
rand.arc = random_network(g.arc, net.df.arc)

#Merge data frames to make table
r = data.frame(AveDegree = NA, ClustCoef = rand.ant[[1]][2,2], AvePathLen = rand.ant[[1]][2,3], Density = NA, Diameter = NA, Modularity = rand.ant[[1]][2,1], Nodes = net.df.ant$Nodes[1], Edges = net.df.ant$Edges[1])

ant.tab = rbind(round(net.df.ant[1,], 3), r, round(net.df.ant[2,],3), round(net.df.ant[4,], 3), round(net.df.ant[3,], 3))

rownames(ant.tab) <- c("Full", "Random Network Mean (SD)", "Abundant", "Intermediate", "Rare")
df1 <- tibble::rownames_to_column(ant.tab, "Network")
df1 = cbind(rep("Antarctic", 5), df1)
names(df1)[1] = "Pole"

r = data.frame(AveDegree = NA, ClustCoef = rand.arc[[1]][2,2], AvePathLen = rand.arc[[1]][2,3], Density = NA, Diameter = NA, Modularity = rand.arc[[1]][2,1], Nodes = net.df.arc$Nodes[1], Edges = net.df.arc$Edges[1])

arc.tab = rbind(round(net.df.arc[1,], 3), r, round(net.df.arc[2,],3), round(net.df.arc[4,], 3), round(net.df.arc[3,], 3))

rownames(arc.tab) <- c("Full", "Random Network Mean (SD)", "Abundant", "Intermediate", "Rare")
df2 <- tibble::rownames_to_column(arc.tab, "Network")
df2 = cbind(rep("Arctic", 5), df2)
names(df2)[1] = "Pole"

final.net.df = rbind(df1, df2)

#add edge/node ratio
final.net.df$EdgeNodeRatio = round(final.net.df$Edges/final.net.df$Nodes, 3)

#10. Scale-free characteristics analysis

# Function for finding if our network exhibits scale-free characteristics
# Code taken from http://chengjun.github.io/web_data_analysis/demo2_simulate_networks/

scale_free_characteristics <- function(graph){
  
  cor_g <- graph
  
  # 1) power law distribution function
  fit_power_law = function(graph, network) {
    
    set.seed(666)
    # calculate degree
    d = degree(graph, mode = "all")
    dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
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
  par(mfrow=c(1,3), mai = c(0.2, 0.2, 0.2, 0.2))
  plot(real.res[[1]] ~ real.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "A")
  curve(power.law.fit1, col = "red", add = T, n = length(real.res[[4]]))
  text(x=min(real.res[[2]]), y=min(real.res[[1]])*1.3, adj = 0, labels=paste0("alpha = ", round(real.res[[5]], 2)))
  text(x=min(real.res[[2]]), y=min(real.res[[1]]), adj = 0, labels=paste0("R2 = ", round(real.res[[6]], 2)))
  plot(scale.res[[1]] ~ scale.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "B")
  curve(power.law.fit2, col = "red", add = T, n = length(scale.res[[4]]))
  text(x=min(scale.res[[2]]), y=min(scale.res[[1]])*1.3, adj = 0,labels=paste0("alpha = ", round(scale.res[[5]], 2)))
  text(x=min(scale.res[[2]]), y=min(real.res[[1]]), adj = 0,labels=paste0("R2 = ", round(scale.res[[6]], 2)))
  plot(random.res[[1]] ~ random.res[[2]], log = "xy", xlab = "Degree (log)", ylab = "Probability (log)",
       col = 1, main = "C")
  curve(power.law.fit3, col = "red", add = T, n = length(random.res[[4]]))
  text(x=min(random.res[[2]]), y=min(random.res[[1]])*1.3, adj = 0,labels=paste0("alpha = ", round(random.res[[5]], 2)))
  text(x=min(random.res[[2]]), y=min(random.res[[1]]), adj = 0,labels=paste0("R2 = ", round(random.res[[6]], 2)))
  
  return(res.list=list(real.res[[6]], scale.res[[6]], random.res[[6]]))
}


ant.scale.ran.res = scale_free_characteristics(g.ant)
arc.scale.ran.res = scale_free_characteristics(g.arc)

#functions for filling in my results
r_val_fit <- function(r_val){
  if (r_val > 0.7){
    out.note = "was a good fit"
  } else if (r_val > 0.4){
    out.note = "was an adequate fit"
  } else if (r_val <= 0.4){
    out.note = "was a poor fit"
  }
  return(out.note)
}

r_val_rand <- function(r_val){
  if (r_val > 0.4){
    out.note = "was"
  } else if (r_val <= 0.4){
    out.note = "was not"
  }
  return(out.note)
}

r_val_rand2 <- function(r_val1, r_val2){
  if (r_val1 > r_val2){
    out.note = "better"
  } else if (r_val1 < r_val2){
    out.note = "worse"
  }
  return(out.note)
}


#The degree distribution of the Antarctic network `r r_val_fit(ant.scale.ran.res[[1]])` to the power-law model (R2=`r round(ant.scale.ran.res[[1]], 3)`), indicating that the network `r r_val_rand(ant.scale.ran.res[[1]])` non-random and with scale-free characteristics. It performed `r r_val_rand2(ant.scale.ran.res[[1]], ant.scale.ran.res[[3]])` than the random network (R2=`r round(ant.scale.ran.res[[3]], 3)`), suggesting it `r r_val_rand(ant.scale.ran.res[[1]])` non-random and scale free (Fig. 1). The degree distribution of the Arctic network `r r_val_fit(arc.scale.ran.res[[1]])` to the power-law model (R2=`r round(arc.scale.ran.res[[1]], 3)`), indicating that the network `r r_val_rand(arc.scale.ran.res[[1]])` non-random and with scale-free characteristics. It performed `r r_val_rand2(arc.scale.ran.res[[1]], arc.scale.ran.res[[3]])` than the random network (R2=`r round(arc.scale.ran.res[[3]], 3)`), suggesting it `r r_val_rand(arc.scale.ran.res[[1]])` non-random and scale free (Fig. 2). Compared with the 1000 Erdös-Réyni random networks, both the empirical networks exhibited greater modularity, clustering coefficient and  average path length (Table 1). This indicates the Arctic and Antarctic networks had “small world” properties and modular structure.

#In the Antarctic network, `r round(final.net.df$Nodes[3]/final.net.df$Nodes[1]*100, 2)`% of nodes were classified as abundant, `r round(final.net.df$Nodes[4]/final.net.df$Nodes[1]*100, 2)`% of nodes were classified as intermediate and `r round(final.net.df$Nodes[5]/final.net.df$Nodes[1]*100, 2)`% as rare. In the Arctic network, `r round(final.net.df$Nodes[8]/final.net.df$Nodes[6]*100, 2)`% of nodes were classified as abundant, `r round(final.net.df$Nodes[9]/final.net.df$Nodes[6]*100, 2)`% of nodes were classified as intermediate and `r round(final.net.df$Nodes[10]/final.net.df$Nodes[6]*100, 2)`% as rare.


#11. Get taxonomy of nodes + report

network_tax <- function(phylo, node_df){
  
  df <- node_df
  #get taxa info of ASVs
  taxonomy <- data.frame(tax_table(phylo))
  network_taxa <- subset(taxonomy, rownames(taxonomy) %in% rownames(df))
  network_taxa$Subcommunity = node_df$category
  
  out <- network_taxa %>%
    dplyr::group_by(Subcommunity, Class) %>% 
    dplyr::summarise(n= n())
  
  out.a = subset(out, Subcommunity == "abundant")
  out.i = subset(out, Subcommunity == "intermediate")
  out.r = subset(out, Subcommunity == "rare")
  
  out.a$perc = round(out.a$n/sum(out.a$n)*100, 2)
  out.i$perc = round(out.i$n/sum(out.i$n)*100, 2)
  out.r$perc = round(out.r$n/sum(out.r$n)*100, 2)
  
  return(res.list = list(network_taxa, out.a, out.i, out.r))
  
}

ant.net.tax = network_tax(ant.ps, node.prop.ant)
ant.df1 = ant.net.tax[[2]]
#get classes and relative abundance of top three + NA
ant.x1 = na.omit(ant.df1)
ant.data_new1 <- ant.x1[order(ant.x1$perc, decreasing = TRUE), ]
ant.y1 = subset(ant.df1, is.na(Class))
ant.df2 = ant.net.tax[[3]]
#get classes and relative abundance of top three + NA
ant.x2 = na.omit(ant.df2)
ant.data_new2 <- ant.x2[order(ant.x2$perc, decreasing = TRUE), ]
ant.y2 = subset(ant.df2, is.na(Class))
ant.df3 = ant.net.tax[[4]]
#get classes and relative abundance of top three + NA
ant.x3 = na.omit(ant.df3)
ant.data_new3 <- ant.x3[order(ant.x3$perc, decreasing = TRUE), ]
ant.y3 = subset(ant.df3, is.na(Class))

arc.net.tax = network_tax(arc.ps, node.prop.arc)
arc.df1 = arc.net.tax[[2]]
#get classes and relative abundance of top three + NA
arc.x1 = na.omit(arc.df1)
arc.data_new1 <- arc.x1[order(arc.x1$perc, decreasing = TRUE), ]
arc.y1 = subset(arc.df1, is.na(Class))
arc.df2 = arc.net.tax[[3]]
#get classes and relative abundance of top three + NA
arc.x2 = na.omit(arc.df2)
arc.data_new2 <- arc.x2[order(arc.x2$perc, decreasing = TRUE), ]
arc.y2 = subset(arc.df2, is.na(Class))
arc.df3 = arc.net.tax[[4]]
#get classes and relative abundance of top three + NA
arc.x3 = na.omit(arc.df3)
arc.data_new3 <- arc.x3[order(arc.x3$perc, decreasing = TRUE), ]
arc.y3 = subset(arc.df3, is.na(Class))


#Within the Antarctic dataset the abundant community was dominated by `r ant.data_new1$Class[[1]]` 
#(`r ant.data_new1$perc[[1]]`%), `r ant.data_new1$Class[[2]]` 
#(`r ant.data_new1$perc[[2]]`%) and `r ant.data_new1$Class[[3]]` (`r ant.data_new1$perc[[3]]`%).
#`r ant.y1$perc`% of the abundant community could not be classified at the class-level. 
#The intermediate community was dominated by `r ant.data_new2$Class[[1]]` (`r ant.data_new2$perc[[1]]`%), 
#`r ant.data_new2$Class[[2]]` (`r ant.data_new2$perc[[2]]`%) and `r ant.data_new2$Class[[3]]` 
#(`r ant.data_new2$perc[[3]]`%). `r ant.y2$perc`% of the intermediate community could not be 
#classified at the class-level. The rare community was dominated by `r ant.data_new3$Class[[1]]` 
#(`r ant.data_new3$perc[[1]]`%), `r ant.data_new3$Class[[2]]` (`r ant.data_new3$perc[[2]]`%) and 
#`r ant.data_new3$Class[[3]]` (`r ant.data_new3$perc[[3]]`%). `r ant.y3$perc`% of the rare community 
#could not be classified at the class-level. Within the Arctic dataset the abundant community was
#dominated by `r arc.data_new1$Class[[1]]` (`r arc.data_new1$perc[[1]]`%), `r arc.data_new1$Class[[2]]` 
#(`r arc.data_new1$perc[[2]]`%) and `r arc.data_new1$Class[[3]]` (`r arc.data_new1$perc[[3]]`%). 
#`r arc.y1$perc`% of the abundant community could not be classified at the class-level. The intermediate 
#community was dominated by `r arc.data_new2$Class[[1]]` (`r arc.data_new2$perc[[1]]`%), 
#`r arc.data_new2$Class[[2]]` (`r arc.data_new2$perc[[2]]`%) and `r arc.data_new2$Class[[3]]` 
#(`r arc.data_new2$perc[[3]]`%). `r arc.y2$perc`% of the intermediate community could not be classified 
#at the class-level. The rare community was dominated by `r arc.data_new3$Class[[1]]` (`r arc.data_new3$perc[[1]]`%), 
#`r arc.data_new3$Class[[2]]` (`r arc.data_new3$perc[[2]]`%) and `r arc.data_new3$Class[[3]]` 
#(`r arc.data_new3$perc[[3]]`%). `r arc.y3$perc`% of the rare community could not be classified at the class-level. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
#12. Find connections between subcommunities
connections_between_communities<- function(graph, abun_df){
  
  # cor_g <- get_my_network(phylo1, phylo2, perc)
  # output <- network_data_prep(phylo1, phylo2, perc)
  
  cor_g = graph
  
  rare <- subset(abun_df, (abun_df$Abundance == "RT"))
  intermediate <- subset(abun_df, (abun_df$Abundance == "MT"))
  abundant <- subset(abun_df, (abun_df$Abundance == "AT"))
  
  #extract the adjacency matric from the simplified igraph
  mat <- as.data.frame(as_adjacency_matrix(cor_g, sparse = FALSE))
  
  
  #make a dataframe with all the pairs of ASVs with significant correlations
  x <- tidyr::gather(mat) #gather data 
  x$ASV_Two <- rep(colnames(mat), nrow(mat)) #add second asv
  x[x==0] <- NA #set 0 to NA
  x2<-x[complete.cases(x),]#remove rows with NA
  x2 <- subset(x2, select = -c(value)) #remove value column
  colnames(x2) <- c("ASV_One", "ASV_Two") #rename columns
  
  #this leaves us with repeat pairs so we need to remove them
  #remove rows that have the repeated pairs of ASVs
  #create a new dataframe with var 1 as ASV One and var2 as ASV Two
  dat <- data.frame(var1 = x2$ASV_One,var2 = x2$ASV_Two, cor = rep(1, nrow(x2)))
  #remove any rows that have the same two ASVs as another row
  dat1 <- dat[!duplicated(apply(dat,1,function(x) paste(sort(x),collapse=''))),]
  #remove cor column
  dat1 <- dat1[,-c(3)]
  
  #get abundance of ASV One
  for (i in 1:nrow(dat1)){
    if(dat1$var1[i] %in% rare$ID){
      dat1$var3[i] <- "rare"
    } else if(dat1$var1[i] %in% abundant$ID){
      dat1$var3[i] <- "abundant"
    } else {
      dat1$var3[i] <- "intermediate"
    }
  }
  
  #get abundance of ASV Two
  for (i in 1:nrow(dat1)){
    if(dat1$var2[i] %in% rare$ID){
      dat1$var4[i] <- "rare"
    } else if(dat1$var2[i] %in% abundant$ID){
      dat1$var4[i] <- "abundant"
    } else {
      dat1$var4[i] <- "intermediate"
    }
  }
  
  #rename the columns
  colnames(dat1) <- c("ASV_One", "ASV_Two", "Abund_One", "Abund_Two")
  
  #table the data
  community_links_df <- data.frame(table(dat1[,3:4]))
  
  return(community_links_df)
  
}

com.links.ant = connections_between_communities(g.ant, ant.df)
com.links.arc = connections_between_communities(g.arc, arc.df)

#13. Find differences between subcommunities

#node_df = node.prop.ant

signif_dif_nodes <- function(node_df){
  
  df <- node_df
  
  # I reorder the groups order : I change the order of the factor 
  df$category <- factor(df$category , levels=c("abundant", "intermediate", "rare"))
  
  #Degree
  my_comparisons <- list( c("abundant", "intermediate"), c("intermediate", "rare"), c("abundant", "rare") )
  
  p1 <- ggboxplot(df, x = "category", y = "degree",
                  fill = "category")+ 
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_brewer(palette="Set1")
  p1 = p1 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p1)$y$range$range[2]*1.1)  # Add global p-value
  
  #Between
  p2 <- ggboxplot(df, x = "category", y = "betweenness",
                  fill = "category")+
    #ylim(0, 3000)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_brewer(palette="Set1")
  p2 = p2 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p2)$y$range$range[2]*1.1)  # Add global p-value
  
  #Closeness
  # df.close = df[order(df$closeness, decreasing = TRUE), ]
  # df.close = df.close[27:nrow(df.close),]
  p3 <- ggboxplot(df, x = "category", y = "closeness",
                  fill = "category")+ 
    #ylim(0, 0.6)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_brewer(palette="Set1")
  p3 = p3 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p3)$y$range$range[2]*1.1)  # Add global p-value
  
  #Eigenvector
  p4 <- ggboxplot(df, x = "category", y = "eigenvector",
                  fill= "category")+ 
    #ylim(0, 0.01)+
    stat_compare_means(comparisons = my_comparisons)+
    scale_fill_brewer(palette="Set1")
  p4 = p4 + # Add pairwise comparisons p-value
    stat_compare_means(label.y = layer_scales(p4)$y$range$range[2]*1.1)  # Add global p-value
  
  
  
  return(res.list=list(p1, p2, p3, p4))
  
}

ant.plots = signif_dif_nodes(node.prop.ant)
arc.plots = signif_dif_nodes(node.prop.arc)

prow <- plot_grid(
  ant.plots[[1]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  ant.plots[[2]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  ant.plots[[3]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  ant.plots[[4]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  arc.plots[[1]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  arc.plots[[2]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  arc.plots[[3]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  arc.plots[[4]] + theme(legend.position="none", axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()),
  #align = 'vh',
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  hjust = -1,
  nrow = 2
)

legend_b <- get_legend(ant.plots[[1]] + 
                         guides(color = guide_legend(nrow = 1, override.aes = list(size = 3),title=NULL)) +
                         theme(legend.position = "bottom",
                               axis.text.y=element_text(size=7),
                               axis.title.y=element_text(size=7,face="bold"),
                               legend.text = element_text(size=7), 
                               legend.title = element_blank())
)

final.p = plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

final.p

#14. Identify keystone taxa
# phylo = ant.ps
# node_df = node.prop.ant
# num = 20
# pole = "Antarctic"

get_keys <- function(phylo, node_df, num, pole){
  
  #Get taxonomy of most highly connected/influential nodes
  tax.df = network_tax(phylo, node_df)
  
  #remove abundance category 
  node_edit = node_df[1:ncol(node_df)-1]
  
  #standardize the values by transforming into Z scores
  z_scores = as.data.frame(sapply(node_edit, function(df) (df-mean(df))/sd(df)))
  row.names(z_scores) = row.names(node_edit)
  
  #Get keystone taxa score - highest degree, closeness and lowest betweenness 
  keystone_scores = as.data.frame(z_scores$degree+z_scores$closeness-z_scores$betweenness)
  row.names(keystone_scores) = row.names(node_edit)
  keystone_scores$category = node_df$category
  
  #Add taxonomy 
  new.tax.df = cbind(keystone_scores, tax.df[[1]])
  
  high.z <- new.tax.df[order(-new.tax.df$`z_scores$degree + z_scores$closeness - z_scores$betweenness`),] #sort by number of degrees
  high.z.keep = head(high.z, num)
  high.z.keep2 = high.z.keep[,1:ncol(high.z.keep)-1]
  names(high.z.keep2) = c("Score", "Subcommunity", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  write.csv(high.z.keep, paste0("../results/", pole , "-keystone-taxa.csv"))
  
  return(high.z.keep2)
  
}

# phylo = ant.ps
# node_df = node.prop.ant
# pole = "Antarctic"

get_keys2 <- function(phylo, node_df, num1, num2, pole){
  
  #Get taxonomy of most highly connected/influential nodes
  tax.df = network_tax(phylo, node_df)
  
  #remove abundance category 
  node_edit = node_df[1:ncol(node_df)-1]
  new.tax.df = cbind(node_edit, tax.df[[1]])
  
  #Key taxa with more than num degrees and less than 5000 betweenness centrality
  node_edit_final = new.tax.df[new.tax.df$degree > num1 & new.tax.df$betweenness < num2,]
  
  write.csv(node_edit_final, paste0("../results/", pole , "-keystone-taxa.csv"))
  
  return(node_edit_final)
  
}

ant.keys = get_keys2(ant.ps, node.prop.ant, 30, 5000, "Antarctic")
arc.keys = get_keys2(arc.ps, node.prop.arc, 30, 5000, "Arctic")

table(arc.keys$Class)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           