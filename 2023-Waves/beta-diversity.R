rm(list=ls())
graphics.off()

library(ggplot2)
library(phyloseq)
# library(patchwork)
# library(ggpubr)
# library(glue)
# library(tibble) # for add_column
# library(reshape2)
# library(cowplot)
# library(betapart)
# library(dplyr)

# 2. Import data
ps.pro <- readRDS("../data/16S-phylo-object-rel.rds") 
ps.euk <- readRDS("../data/18S-phylo-object-rel.rds") 
#remove blanks
ps.euk = subset_samples(ps.euk, SampleType != "Control")
ps.euk = filter_taxa(ps.euk, function(x) sum(x) > 0, TRUE)

#remove the 2 most abundant taxa from ps.euk
#which ASVs have the most counts?
tax = data.frame(tax_table(ps.euk))
x = data.frame(sort(taxa_sums(ps.euk), decreasing = FALSE))
names(x) = "value"
y = x %>% 
  top_n(2, value)
z = tax[rownames(tax) %in% rownames(y),]
allTaxa = taxa_names(ps.euk)
allTaxa <- allTaxa[!(allTaxa %in% rownames(z))]
ps.euk.rm = prune_taxa(allTaxa, ps.euk)

# ps = ps.pro
# d = "wunifrac"
# lab="Weighted Unifrac Distances"

beta_function <- function(ps, d, lab){
  
  #get distance matrix and melt
  dis <- phyloseq::distance(ps, method = d) 
  dis.melt = melt(as.matrix(dis))
  names(dis.melt) = c("Sample1", "Sample2", "Dis")
  #remove samples being compared to themseves
  df = dis.melt[dis.melt['Sample1'] != dis.melt['Sample2'],]
  
  #get sample data
  sd = data.frame(sample_data(ps))
  sd = sd %>%
    select("SampleID", "SampleType") %>%
    mutate_if(is.factor,as.character)
  
  colnames(sd) = c("Sample1", "Type1")
  df = left_join(df, sd, by = "Sample1")
  
  colnames(sd) = c("Sample2", "Type2")
  df = left_join(df, sd, by = "Sample2")
  
  #only keep rows where samples from the same type are compared
  final.df = df %>% 
    filter(str_sub(Type1) == (str_sub(Type2)))
  
  #plot
  
  my_comparisons <- list( c("Interface", "Meltwater"), c("Interface", "Surface"), c("Surface", "Meltwater") )
  
  p = ggplot(final.df, aes(x = Type1, y = Dis, fill=Type1)) +
    theme_bw() +
    geom_boxplot() +
    scale_color_identity() +
    stat_compare_means(comparisons = my_comparisons)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20)) +
    ylab(lab) +
    xlab("Sample Type")
  
  p
  
  return(p)
  
}


euk.p = beta_function(ps.euk, "wunifrac", "Weighted Unifrac Distances")
euk.p

euk.p.rm = beta_function(ps.euk.rm, "wunifrac", "Weighted Unifrac Distances")
euk.p.rm

#save plot

pdf("../results/wunifrac-boxplots.pdf", height=5, width=10)
print(euk.p)
dev.off()

