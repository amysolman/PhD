#Control Samples

#Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)
# devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic") #for top_taxa
library(patchwork) #for formatting my plots
library(cowplot) #for ggdraw

#Import data

#prokaryotes
pro.DNA <- readRDS("../results/mt-dna-16S-phylo-object.rds") 
pro.RNA <- readRDS("../results/mt-rna-16S-phylo-object.rds") 

#eukaryotes without micrometazoans
euk.DNA <- readRDS("../results/mt-dna-18S-phylo-object-micro-remove.rds") 
euk.RNA <- readRDS("../results/mt-rna-18S-phylo-object.rds") 

#load plotting colours
pro.cols <- read.csv("../data/thesis-plotting-colours-prokaryotes.csv")
euk.cols <- read.csv("../data/thesis-plotting-colours-eukaryotes.csv")

#Histograms of environmental samples and control samples to see how many reads are in our controls

#function
sample_hist_func <- function(ps, data_type, label_letter){
  
  df = data.frame(otu_table(ps))
  count_tab_df = as.data.frame(colSums(df)) #make into df
  count_tab_df$Sample = rownames(count_tab_df) #make sample names a column
  names(count_tab_df) = c("Reads", "Sample")
  p<-ggplot(data=count_tab_df, aes(x=Sample, y=Reads)) +
    geom_bar(stat="identity")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.title = element_text(size=20))+
    ggtitle(paste0(label_letter, ": ", data_type, " Control Read Depth"))
  
  p
  
  return(p)
}

#run function
pro.dna.p = sample_hist_func(pro.DNA, "DNA", label_letter = "A1")
pro.dna.p
pro.rna.p = sample_hist_func(pro.RNA, "RNA", label_letter = "A3")
pro.rna.p
euk.dna.p = sample_hist_func(euk.DNA, "DNA", label_letter = "B1")
euk.dna.p
euk.rna.p = sample_hist_func(euk.RNA, "RNA", label_letter = "B3")
euk.rna.p

#relative abundance bar plots of control samples
# ps = pro.RNA
# col.df = pro.cols
# data_type = "RNA"

control_com_plot <- function(ps, col.df, data_type, leg_col, label_letter){
  
#subset to control samples only
ps.c = subset_samples(ps, Habitat == "Control")
#remove samples/taxa with zero counts
ps.c = prune_samples(sample_sums(ps.c) >0, ps.c)
ps.c = filter_taxa(ps.c, function(x) sum(x) > 0, TRUE)

#STEP ONE: find the 10 most abundant phyla in each habitat
#setting to top 100 taxa because we want to include everything
top.phy <- top_taxa(ps.c, 
                    tax_level = "Phylum", 
                    n_taxa = 100,
                    grouping = c("Habitat"))

top.phy2 = top.phy[[2]]
length(unique(top.phy2$Phylum))
phy.keep = unique(top.phy2$Phylum)
phy.keep


#STEP TWO: modify tax names so it looks better when we plot them
ps.m = ps.c
x = data.frame(tax_table(ps.m)) %>% 
  mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
  mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
  mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
  mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
  mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
  mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus))

#replace anything that says NA with Genus Unknown
y = x
y$Phylum <- gsub("P.NA", "(Genus Unknown)", y$Phylum)
y$Class <- gsub("P.NA C.NA", "(Genus Unknown)", y$Class)
y$Class <- gsub("C.NA", "(Genus Unknown)", y$Class)
y$Order <- gsub("P.NA C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("C.NA O.NA", "(Genus Unknown)", y$Order)
y$Order <- gsub("O.NA", "(Genus Unknown)", y$Order)
y$Family <- gsub("P.NA C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("C.NA O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("O.NA F.NA", "(Genus Unknown)", y$Family)
y$Family <- gsub("F.NA", "(Genus Unknown)", y$Family)
y$Genus <- gsub("P.NA C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("C.NA O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("O.NA F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("F.NA G.NA", "(Genus Unknown)", y$Genus)
y$Genus <- gsub("G.NA", "(Genus Unknown)", y$Genus)

tax_table(ps.m) = as.matrix(y)

#find the three most abundant genus within each phyla
#get our data
data <-
  ps.m %>%
  tax_glom("Genus") %>%
  psmelt() %>%
  as_tibble()

#get 3 most abundant genus by phylum
out = data %>%
  filter(Phylum %in% phy.keep) %>%
  group_by(Phylum, Genus) %>%
  dplyr::summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)%>%
  top_n(n = 1)

#does our new data frame have only the top phyla in it? Should all be TRUE
sort(unique(out$Phylum)) == sort(phy.keep)

#sort the out dataframe and add colours
out.col = out[with(out, order(Phylum)), ]
out.col$Genus = paste0(out.col$Phylum, ": ", out.col$Genus)
save.cols <- vector()
df2keep = data.frame(Phylum=as.character(), Genus=as.character(), Abundance=as.numeric())

for (i in 1:length(unique(out.col$Phylum))){ #for each unique phyla
  
  #get our phyla data only
  mini.df = out.col[out.col$Phylum == unique(out.col$Phylum)[i],]
  
  #which colours are assigned to that phyla?
  col1 = col.df[col.df$Phylum == unique(mini.df$Phylum),]$Colour1
  col2 = col.df[col.df$Phylum == unique(mini.df$Phylum),]$Colour2
  
  #get colour function for that phylum
  cols.fun = colorRampPalette(c(col1, col2))
  
  #get colours for each class
  save.cols = c(save.cols, cols.fun(nrow(mini.df)+1))
  
  df1 = rbind(mini.df, data.frame(Phylum=unique(out.col$Phylum)[i], Genus=paste0(unique(out.col$Phylum)[i], ": Other")))
  df2keep = rbind(df2keep, df1)
}

#add phya and genera to plot with their assigned colours
df2keep$col = save.cols

#save data as data2plot for additional wrangling
data2plot = data

#sanity check
#number of all phyla in our dataframe
length(unique(data2plot$Phylum))

#number of all genera in our dataframe
length(unique(data2plot$Genus))

#in the main dataframe with all our plotting data, create an additional colour with genera and phyla data combines
data2plot$Genus2 = paste0(data2plot$Phylum, ": ", data2plot$Genus)

#if they are not in our keep pile turn them into "Other"
for (i in 1:nrow(data2plot)){ #for each row of our main dataframe
  
  if(!data2plot$Genus2[i] %in% out.col$Genus){
    #if the genus' in the main dataframe are not in the top 3 of their phyla they are classed as "other"
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": Other")
    
  } else if (data2plot$Genus2[i] %in% out.col$Genus){ 
    #else they are just given the phyla name with the genus name
    data2plot$Genus[i] <- paste0(data2plot$Phylum[i], ": ", data2plot$Genus[i])
    
  }
}

#remove unwanted phyla
for (i in 1:nrow(data2plot)){
  if(!data2plot$Phylum[i] %in% unique(out$Phylum)){ #if phyla aren't in our keep selection then replace their "genus" with Other
    data2plot$Genus[i] <- "Other"
  }
}

#add those colours to the corresponding row in our main dataframe
#add "Other" to out.col
out.col2 = rbind(df2keep, data.frame(Phylum="Other", Genus="Other", col="#000000"))
#check we have the same values
length(intersect(out.col2$Genus, data2plot$Genus)) == length(unique(data2plot$Genus)) #should = TRUE
data2plot2 = merge(data2plot, out.col2[, c("Genus", "col")], by="Genus")

#plot data
df = data2plot2
df$Habitat <- factor(df$Habitat, levels = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite"))
df$Data = data_type

#sort dataframe according to phylum
df.plot = df[with(df, order(Phylum)), ]
#set phylum as factor 
levs = unique(df.plot$Phylum)
levs = levs[levs != "Other"]
df.plot$Phylum <- factor(df.plot$Phylum, levels = c(levs, "Other"))

#set Genus as factors to keep their order when plotting
#sort dataframe according to phylum
gen.df = out.col2[with(out.col2, order(Genus)), ]
levs = gen.df$Genus
levs = levs[levs != "Other"]
df.plot$Genus <- factor(df.plot$Genus, levels = out.col2$Genus)

#sort out colours 
col <- as.character(df.plot$col)
names(col) <- as.character(df.plot$Genus)

df.plot2 = 
  df.plot %>%
  group_by(Genus, Habitat, Data, Sample) %>%
  dplyr::summarise(across(c(Abundance), sum))

#make axis title
# my_y_title <- expression(paste(italic("tDNA"), " Relative Abundance"))
#control_labels = c("Control" = "Control")
p = ggplot(df.plot2, aes(fill=Genus, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity", alpha=1)+
  # facet_grid(cols = vars(Habitat), scales = "free", space = "free",
  #            labeller = labeller(Habitat = control_labels))+
  theme_bw()+
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(ncol= leg_col, byrow=FALSE))+
  theme(legend.position="bottom",
        #axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.title=element_blank(), 
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        plot.title = element_text(size=20))+
  ggtitle(paste0(label_letter, ": ", data_type, " Control Community Profiles"))

p

return(p)

}

#run function
pro.dna.com = control_com_plot(pro.DNA, pro.cols, "DNA", leg_col = 3, label_letter = "A2")
pro.dna.com
pro.rna.com = control_com_plot(pro.RNA, pro.cols, "RNA", leg_col = 3, label_letter = "A4")
pro.rna.com
euk.dna.com = control_com_plot(euk.DNA, euk.cols, "DNA", leg_col = 4, label_letter = "B2")
euk.dna.com
euk.rna.com = control_com_plot(euk.RNA, euk.cols, "RNA", leg_col = 4, label_letter = "B4")
euk.rna.com

#titles
title1 <- ggdraw() + draw_label("A: Prokaryote", size=40, y=0, x=0, vjust=0, hjust=0)
title2 <- ggdraw() + draw_label("B: Microbial Eukaryote", size=40, y=0, x=0, vjust=0, hjust=0)
all.plots = (title1 | title2) / (plot_spacer() | plot_spacer()) / (pro.dna.p | euk.dna.p) / (pro.dna.com | euk.dna.com) / (pro.rna.p | euk.rna.p) / (pro.rna.com | euk.rna.com) + plot_layout(heights = c(0.1, 0.1, 1, 1, 1, 1))
all.plots
pdf("../results/mt-all-control-plots.pdf", height=25, width=18)
print(all.plots)
dev.off()

#Remove control samples and Malassezia from our phyloseq objects

pro.DNA.e = subset_samples(pro.DNA, Habitat != "Control")
pro.DNA.e = prune_samples(sample_sums(pro.DNA.e) >0, pro.DNA.e)
pro.DNA.e = filter_taxa(pro.DNA.e, function(x) sum(x) > 0, TRUE)

pro.RNA.e = subset_samples(pro.RNA, Habitat != "Control")
pro.RNA.e = prune_samples(sample_sums(pro.RNA.e) >0, pro.RNA.e)
pro.RNA.e = filter_taxa(pro.RNA.e, function(x) sum(x) > 0, TRUE)

euk.DNA.e = subset_samples(euk.DNA, Habitat != "Control")
euk.DNA.e = prune_samples(sample_sums(euk.DNA.e) >0, euk.DNA.e)
euk.DNA.e = filter_taxa(euk.DNA.e, function(x) sum(x) > 0, TRUE)
euk.DNA.e = subset_taxa(euk.DNA.e, !Genus %in% c("Malassezia"))

euk.RNA.e = subset_samples(euk.RNA, Habitat != "Control")
euk.RNA.e = prune_samples(sample_sums(euk.RNA.e) >0, euk.RNA.e)
euk.RNA.e = filter_taxa(euk.RNA.e, function(x) sum(x) > 0, TRUE)
euk.RNA.e = subset_taxa(euk.RNA.e, !Genus %in% c("Malassezia"))

#save our new phyloseq objects
saveRDS(pro.DNA.e, "../results/mt-dna-16S-phylo-object-no-controls.rds")
saveRDS(pro.RNA.e, "../results/mt-rna-16S-phylo-object-no-controls.rds")
saveRDS(euk.DNA.e, "../results/mt-dna-18S-phylo-object-no-controls.rds")
saveRDS(euk.RNA.e, "../results/mt-rna-18S-phylo-object-no-controls.rds")
