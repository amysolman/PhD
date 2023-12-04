#Control Samples

#Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)

#Import data

#prokaryotes
pro.DNA <- readRDS("../results/16S-phylo-object.rds") 
pro.RNA <- readRDS("../results/mt-pro-phylo-object.rds") 

#eukaryotes without micrometazoans
euk.DNA <- readRDS("../results/18S-phylo-object-micro-keep.rds") 
euk.RNA <- readRDS("../results/mt-euk-phylo-object.rds") 

#load plotting colours
pro.cols <- read.csv("../../../thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../../../thesis-plotting-colours-eukaryotes.csv")

#remove control samples
pro.DNA.c = subset_samples(pro.DNA, Habitat2 == "Control")
pro.DNA.c = prune_samples(sample_sums(pro.DNA.c) >0, pro.DNA.c) #remove samples with zero counts
pro.DNA.c = filter_taxa(pro.DNA.c, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
pro.RNA = subset_samples(pro.RNA, Habitat2 == "Control")
pro.RNA = prune_samples(sample_sums(pro.RNA) >0, pro.RNA) #remove samples with zero counts
pro.RNA = filter_taxa(pro.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
euk.DNA.c = subset_samples(euk.DNA, Habitat2 == "Control")
euk.DNA.c = prune_samples(sample_sums(euk.DNA.c) >0, euk.DNA.c) #remove samples with zero counts
euk.DNA.c = filter_taxa(euk.DNA.c, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
euk.RNA = subset_samples(euk.RNA, Habitat2 == "Control")
euk.RNA = prune_samples(sample_sums(euk.RNA) >0, euk.RNA) #remove samples with zero counts
euk.RNA = filter_taxa(euk.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts

#How many reads?
sum(taxa_sums(pro.DNA.c)) #158775
sum(taxa_sums(pro.RNA)) #594
sum(taxa_sums(euk.DNA.c)) #39398
sum(taxa_sums(euk.RNA)) #19

#How many ASVs?
#16S DNA.c
ntaxa(pro.DNA.c) #How many prokaryote ASVs? 103 ASVs
# arc.DNA.c <- subset_taxa(pro.DNA.c, Domain %in% c("Archaea")) #of which are archaea? NONE
# arc.DNA.c = prune_samples(sample_sums(arc.DNA.c) >0, arc.DNA.c) #remove samples with zero counts
# arc.DNA.c = filter_taxa(arc.DNA.c, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
# ntaxa(arc.DNA.c) #number of ASVs? 
# ntaxa(pro.DNA.c) - ntaxa(arc.DNA.c) #so number of bacterial ASVs only? 
# arc.DNA.c.tax = data.frame(tax_table(arc.DNA.c)) 
#16S RNA
ntaxa(pro.RNA) #How many prokaryote NTUs? 193 NTUs
# arc.RNA <- subset_taxa(pro.RNA, Domain %in% c("Archaea")) #of which are archaea? NONE
# arc.RNA = prune_samples(sample_sums(arc.RNA) >0, arc.RNA) #remove samples with zero counts
# arc.RNA = filter_taxa(arc.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
# ntaxa(arc.RNA) #number of NTUs
# ntaxa(pro.RNA) - ntaxa(arc.RNA) #so number of bacterial NTUs only?
# arc.RNA.tax = data.frame(tax_table(arc.RNA)) #archaea taxonomy
#18S DNA.c
ntaxa(euk.DNA.c) #How many eukaryote ASVs? 16 ASVs
# mf.DNA.c <- subset_taxa(euk.DNA.c, Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa")) #of which are microfauna? NONE
# mf.DNA.c = prune_samples(sample_sums(mf.DNA.c) >0, mf.DNA.c) #remove samples with zero counts
# mf.DNA.c = filter_taxa(mf.DNA.c, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
# ntaxa(mf.DNA.c) #number of ASVs?
# ntaxa(euk.DNA.c) - ntaxa(mf.DNA.c) #so number of non-microfauna eukaryotic ASVs only?
# mf.DNA.c.tax = data.frame(tax_table(mf.DNA.c)) 
#18S RNA
ntaxa(euk.RNA) #How many eukaryote NTUs? 11 NTUs
# mf.RNA <- subset_taxa(euk.RNA, Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa")) #of which are microfauna? NONE
# mf.RNA = prune_samples(sample_sums(mf.RNA) >0, mf.RNA) #remove samples with zero counts
# mf.RNA = filter_taxa(mf.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
# ntaxa(mf.RNA) #number of NTUs 26
# ntaxa(euk.RNA) - ntaxa(mf.RNA) #so number of bacterial NTUs only?
# mf.RNA.tax = data.frame(tax_table(mf.RNA)) #microfauna taxonomy



#community structure
# ps = pro.DNA.c
# lev = "Class"

#function for pulling out phylum-level relative abundance
com_structure <- function(ps, lev){
  
  glom <- tax_glom(ps, taxrank = lev)
  t.glom = data.frame(otu_table(glom))
  rownames(t.glom) = unlist(as.vector(data.frame(tax_table(glom))[c(lev)]))
  x = sort(rowSums(t.glom))/sum(t.glom)*100
  
  return(x)
  
}

#16S Relative abundances

#16S DNA.c
pro.DNA.coms = com_structure(pro.DNA.c, "Class")
pro.DNA.coms
#16S RNA
pro.RNA.coms = com_structure(pro.RNA, "Class")
pro.RNA.coms

#18S Relative abundances

#18S DNA.c
euk.DNA.coms = com_structure(euk.DNA.c, "Class")
euk.DNA.coms
#18S RNA
euk.RNA.coms = com_structure(euk.RNA, "Class")
euk.RNA.coms

#Malasseziomycetes (human skin fungi) represented 30% of sequences in our 18S DNA.c control samples.
#How abundant is this class in our environmental samples?

euk.DNA.e = subset_samples(euk.DNA, Habitat2 != "Control")
euk.DNA.e = prune_samples(sample_sums(euk.DNA.e) >0, euk.DNA.e) #remove samples with zero counts
euk.DNA.e = filter_taxa(euk.DNA.e, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts

euk.DNA.e.m = subset_taxa(euk.DNA.e, Class %in% c("Malasseziomycetes"))
euk.DNA.e.m = prune_samples(sample_sums(euk.DNA.e.m) >0, euk.DNA.e.m) #remove samples with zero counts
euk.DNA.e.m = filter_taxa(euk.DNA.e.m, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
euk.DNA.e.m
sum(taxa_sums(euk.DNA.e.m))

#Plot control samples

#function

com_plot_function <- function(ps, cols, name){
  
  #modify tax names
  df = data.frame(tax_table(ps))
  for (i in 1:nrow(df)){
    data = unlist(as.vector(df[i,]))
    x = na.locf(data)
    y = names(data[is.na(data)]) #which columns had NA?
    for (j in 1:length(y)){
      col = y[j]
      x[col] = paste0(x[col], " (",  y[j], " Unknown)")
    }
    df[i,] = x
    
  }
  
  
  tax_table(ps) = as.matrix(df)
  
  #how many phyla?
  unique(data.frame(tax_table(ps))$Phylum)
  
  #find the 15 most abundant Phyla
  top.phy <- top_taxa(ps, 
                      tax_level = "Phylum", 
                      n_taxa = 15,
                      grouping = "Habitat")
  top.phy2 = top.phy[[2]]
  length(unique(top.phy2$Phylum))
  phy.keep = unique(top.phy2$Phylum)
  
  #find the three most abundant genus within each phyla
  #get our data
  data <-
    ps %>%
    tax_glom("Class") %>%
    psmelt() %>%
    as_tibble()

  
  #get 3 most abundant genus by phylum
  out = data %>%
    filter(Phylum %in% phy.keep) %>%
    group_by(Phylum, Class) %>%
    dplyr::summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance)%>%
    top_n(n = 3)
  
  
  #sort the out dataframe and add colours
  out.col = out[with(out, order(Phylum)), ]
  out.col$Class = paste0(out.col$Phylum, ": ", out.col$Class)
  
  save.cols <- vector()
  df2keep = data.frame(Phylum=as.character(), Class=as.character(), Abundance=as.numeric())
  
  for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
    
    #get our phyla data only
    mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
    
    #which colours are assigned to that phyla?
    col1 = cols[cols$Phylum == mini.df$Phylum,]$Colour1
    col2 = cols[cols$Phylum == mini.df$Phylum,]$Colour2
    
    #get colour function for that phylum
    cols.fun = colorRampPalette(c(col1, col2))
    
    #get colours for each class
    save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
    
    df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Class=paste0(unique(out.col$Phylum)[i], ": Other")))
    df2keep = rbind(df2keep, df1)
  }
  
  #out phya and classes to plot with their assigned colours
  df2keep$col = save.cols
  
  #save data as data2plot for additional wrangling
  data2plot = data
  
  
  #remove unwanted classes and replace with other
  for (i in 1:nrow(data2plot)){
    if(!data2plot$Class[i] %in% out$Class){
      data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": Other")
    } else if (data2plot$Class[i] %in% out$Class){
      data2plot$Class[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Class[i])
    }
  }
  
  
  #remove unwanted phyla
  for (i in 1:nrow(data2plot)){
    if(!data2plot$Phylum[i] %in% unique(out$Phylum)){
      data2plot$Class[i] <- "Other"
    }
  }

  
  #add those colours to the corresponding row in our main dataframe
  
  #add "Other" to out.col
  out.col2 = rbind(df2keep, data.frame(Phylum="Other", Class="Other", col="#000000"))
  length(unique(data2plot$Class))
  length(unique(out.col2$Class))
  data2plot2 = merge(data2plot, out.col2[, c("Class", "col")], by="Class")
  
  #plot data
  df = data2plot2

  #sort dataframe according to phylum
  df.plot = df[with(df, order(Phylum)), ]
  #set phylum as factor 
  levs = unique(df.plot$Phylum)
  levs = levs[levs != "Other"]
  df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))
  
  #set Class as factors to keep their order when plotting
  #sort dataframe according to phylum
  gen.df = out.col2[with(out.col2, order(Class)), ]
  levs = gen.df$Class
  levs = levs[levs != "Other"]
  df.plot$Class <- factor(df.plot$Class, levels = out.col2$Class)
  
  #sort out colours 
  col <- as.character(df.plot$col)
  names(col) <- as.character(df.plot$Class)
  
  df.plot2 = 
    df.plot %>%
    group_by(Class, Sample) %>%
    dplyr::summarise(across(c(Abundance), sum))
  
  
  p = ggplot(df.plot2, aes(fill=Class, y=Abundance, x=Sample)) + 
    geom_bar(position="fill", stat="identity", alpha=1)+
    ylab("Relative Abundance")+ 
    theme_bw()+
    scale_fill_manual(values = col)+
    guides(fill=guide_legend(ncol=3, byrow=FALSE))+
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          legend.title=element_blank())
  
  p
  
  pdf(paste0("../results/", name, "-control-sample-plot-rel.pdf"), width=15, height=10)
  print(p)
  dev.off()
  
  return(p)
  
}


#16S DNA
pro.DNA.p = com_plot_function(pro.DNA.c, pro.cols, "16S-DNA")
pro.DNA.p

#16S RNA
pro.RNA.p = com_plot_function(pro.RNA, pro.cols, "16S-RNA")
pro.RNA.p

#18S DNA
euk.DNA.p = com_plot_function(euk.DNA.c, euk.cols, "18S-DNA")
euk.DNA.p

#18S RNA
euk.RNA.p = com_plot_function(euk.RNA, euk.cols, "18S-RNA")
euk.RNA.p

#combine plots

p.final = (pro.DNA.p + pro.RNA.p) /
          (euk.DNA.p + euk.RNA.p) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 25))
p.final


pdf("../results/all-community-control-plots.pdf", width=22, height=20)
print(p.final)
dev.off()