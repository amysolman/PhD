###@author: Feng Ju
###@email: richieju520@gmail.com
###@cite Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. 
###@Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants.
###@Environmental Microbiology. 16(8):2421-2432

################## Correlation-based network analysis ###############################
#install.packages("vegan")
#install.packages("igraph")
#install.packages("Hmisc")
rm(list=ls())

load("results/cleaned_files/cleanedfiles_3.RData")


library(vegan)
library(igraph)
library(Hmisc)

co_occurrence_network<-function(matrix,cor.cutoff,p.cutoff){
  
  #is this supposed to be used for the rest of the script?
  # matrix1<-matrix #copy asv table to new variable matrix 1
  # matrix1[matrix1>0]<-1 #make any occurrence of an ASV in a sample more than 0 = 1 (a.k.a turn it into a presence/abudance table)
  
  #correlation analysis based on spearman's co-efficient
  matrix.dist<-rcorr(t(matrix),type="spearman") #transpose the original count table (so now ASVs are columns and samples are rows)
  ###matrix.dist<-rcorr(t(matrix),type="pearson") #comput a matrix of spearman;s rho rank correlation coefficients for all possible pairs of columns (ASVs) in a matrix
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  matrix.cor.p <- p.adjust(matrix.cor.p, method="BH")
  
  #1.Consider positive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor1<-matrix.cor
  matrix.cor1.p<-matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= cor.cutoff)]=0
  matrix.cor1[which(matrix.cor1.p>p.cutoff)]=0
  # delete those rows and columns with sum = 0
  matrix.cor1<-matrix.cor1[which(rowSums(matrix.cor1)!=1),]
  matrix.cor1<-matrix.cor1[,which(colSums(matrix.cor1)!=0)]
  
  #2.Consider netagive cooccurence at given coefficient (-cor.cutoff) and p-value cutoffs
  matrix.cor2<-matrix.cor
  matrix.cor2.p<-matrix.cor.p
  matrix.cor2[which(matrix.cor2 > (-cor.cutoff))]=0
  matrix.cor2[which(matrix.cor2.p>p.cutoff)]=0
  # delete those rows and columns with sum = 0
  matrix.cor2<-matrix.cor2[which(rowSums(matrix.cor2)!=1),]
  matrix.cor2<-matrix.cor2[,which(colSums(matrix.cor2)!=0)]
  
  #3.Consider both positive and netagive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor3<-matrix.cor
  matrix.cor3.p<-matrix.cor.p
  matrix.cor3[which(matrix.cor3>=(-cor.cutoff) & matrix.cor3 <= cor.cutoff)]=0
  matrix.cor3[which(matrix.cor3.p>p.cutoff)]=0
  
  # delete those rows and columns with sum = 0
  matrix.cor3<-matrix.cor3[which(rowSums(matrix.cor3)!=1),]
  matrix.cor3<-matrix.cor3[,which(colSums(matrix.cor3)!=0)]
  
  # generate graph using igraph
  g1<-graph.adjacency(matrix.cor1,weight=T,mode="undirected")
  g1<-simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)
  
  ###g2<-graph.adjacency(matrix.cor2,weight=T,mode="undirected")
  ###g2<-simplify(g2)
  ###V(g2)$label <- V(g2)$name
  ###V(g2)$degree <- degree(g2)
  
  g3<-graph.adjacency(matrix.cor3,weight=T,mode="undirected")
  g3<-simplify(g3)
  V(g3)$label <- V(g3)$name
  V(g3)$degree <- degree(g3)
  
  # append the output into results
  result<-list()
  result$matrix.cor<-matrix.cor
  result$matrix.cor.p<-matrix.cor.p
  
  result$matrix.cor1<-matrix.cor1
  result$graph1<-g1
  
  ###result$matrix.cor2<-matrix.cor2
  ###result$graph2<-g2
  
  result$matrix.cor3<-matrix.cor3
  result$graph3<-g3
  return(result)
}



# Co-occurrence-network-analysis
################## OTU filtering, network generation, topological analysis and export OTU table ###############################
library(igraph)
library(Hmisc)

#Abu=read.table('NW.txt',header=T) #this is where you load your data
#Abu<-as.matrix(Abu)
Abu <- as.matrix(count_table_clean)

###1. Filtering OTUs by occurrence frequency (i.e.,number of samples an OTU is Present)
table<-Abu #copy my abundance table to a new variable called table
table[table>0]<-1 #for every value in the table that is more than 0, make it a 1
table.generalist<-Abu[which(rowSums(table)>=12),] #make a new variable from Abu, only including rows (ASVs) that occure in 12 or more samples
Abu<-table.generalist #make the variable Abu only full of generalise taxa (those in 12 or more samples)

###2. Creating gml files of network (to be visulized in Gephi or Cytoscape)
pattern<-co_occurrence_network(Abu,0.6,0.01)  ## cutoffs for correlation coefficient and P-value 

write.graph(pattern$graph1,'Pos0.6-NW.gml',format='gml')    #network file for positive association
#write.graph(pattern$graph2,'Neg0.6-NW.gml',format='gml')   #network file for negative association (if any) there were no negative associations
write.graph(pattern$graph3,'PosNeg0.6-NW.gml',format='gml') #network file for all association

###3. Calculating network topological properties
g<-pattern$graph1   ###positive network
#g<-pattern$graph2   ###negative network

c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
gd  <- graph.density(g, loops=FALSE)
nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)

node.degree <- degree(g, v = V(g), mode="all")
ad  <- mean(node.degree)

e <- ecount(g)
v <- vcount(g)
global.topology <- data.frame(e,v,cc,spl,md,gd,nd,ad)
write.csv(global.topology, file="Pos0.6-NW-global.topology.csv")

# Node toplogical features
betweenness.centrality <- betweenness(g, v=V(g), 
                                      directed = FALSE, weights = NA,
                                      nobigint = TRUE, normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NA, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NA)

node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Pos0.6-NW-node.topology.csv")

# Ploting node degreee distribution in a log-log plot
degree.df <- data.frame(table(degree=factor(node.degree, levels=seq_len(max(node.degree)))))
degree.df$degree <- as.numeric(as.character(degree.df$degree))

#4. Creating an abundance table for OTUs present in the positive and negative network
my.list1 <- row.names(pattern$matrix.cor1)
###my.list2 <- row.names(pattern$matrix.cor2)

logical1 <- row.names(Abu)  %in% my.list1
###logical2 <- row.names(Abu)  %in% my.list2

tab.subset1 <- subset(Abu,logical1)
###tab.subset2 <- subset(Abu,logical2)

write.table(tab.subset1,'Pos0.6-NW.txt',sep="\t")
###write.table(tab.subset2,'Neg0.6-NW.txt',sep="\t")

