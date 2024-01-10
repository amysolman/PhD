#Community stats

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

#remove control samples
pro.DNA = subset_samples(pro.DNA, Habitat2 != "Control")
pro.DNA = prune_samples(sample_sums(pro.DNA) >0, pro.DNA) #remove samples with zero counts
pro.DNA = filter_taxa(pro.DNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
pro.RNA = subset_samples(pro.RNA, Habitat2 != "Control")
pro.RNA = prune_samples(sample_sums(pro.RNA) >0, pro.RNA) #remove samples with zero counts
pro.RNA = filter_taxa(pro.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
euk.DNA = subset_samples(euk.DNA, Habitat2 != "Control")
euk.DNA = prune_samples(sample_sums(euk.DNA) >0, euk.DNA) #remove samples with zero counts
euk.DNA = filter_taxa(euk.DNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
euk.RNA = subset_samples(euk.RNA, Habitat2 != "Control")
euk.RNA = prune_samples(sample_sums(euk.RNA) >0, euk.RNA) #remove samples with zero counts
euk.RNA = filter_taxa(euk.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts

#How many reads?
sum(taxa_sums(pro.DNA)) #905021
sum(taxa_sums(pro.RNA)) #2015700
sum(taxa_sums(euk.DNA)) #840816
sum(taxa_sums(euk.RNA)) #697354

#How many ASVs?
#16S DNA
ntaxa(pro.DNA) #How many prokaryote ASVs? 2671 ASVs
arc.DNA <- subset_taxa(pro.DNA, Domain %in% c("Archaea")) #of which are archaea?
arc.DNA = prune_samples(sample_sums(arc.DNA) >0, arc.DNA) #remove samples with zero counts
arc.DNA = filter_taxa(arc.DNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
ntaxa(arc.DNA) #number of ASVs? 8
ntaxa(pro.DNA) - ntaxa(arc.DNA) #so number of bacterial ASVs only? 2663
arc.DNA.tax = data.frame(tax_table(arc.DNA)) 
#16S RNA
ntaxa(pro.RNA) #How many prokaryote NTUs? 4174 NTUs
arc.RNA <- subset_taxa(pro.RNA, Domain %in% c("Archaea")) #of which are archaea?
arc.RNA = prune_samples(sample_sums(arc.RNA) >0, arc.RNA) #remove samples with zero counts
arc.RNA = filter_taxa(arc.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
ntaxa(arc.RNA) #number of NTUs 26
ntaxa(pro.RNA) - ntaxa(arc.RNA) #so number of bacterial NTUs only? 4148
arc.RNA.tax = data.frame(tax_table(arc.RNA)) #archaea taxonomy
#18S DNA
ntaxa(euk.DNA) #How many eukaryote ASVs? 825 ASVs
mf.DNA <- subset_taxa(euk.DNA, Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa")) #of which are microfauna?
mf.DNA = prune_samples(sample_sums(mf.DNA) >0, mf.DNA) #remove samples with zero counts
mf.DNA = filter_taxa(mf.DNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
ntaxa(mf.DNA) #number of ASVs? 19
ntaxa(euk.DNA) - ntaxa(mf.DNA) #so number of non-microfauna eukaryotic ASVs only? 814
mf.DNA.tax = data.frame(tax_table(mf.DNA)) 
#18S RNA
ntaxa(euk.RNA) #How many eukaryote NTUs? 1180 NTUs
# mf.RNA <- subset_taxa(euk.RNA, Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa")) #of which are microfauna? NONE
# mf.RNA = prune_samples(sample_sums(mf.RNA) >0, mf.RNA) #remove samples with zero counts
# mf.RNA = filter_taxa(mf.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
# ntaxa(mf.RNA) #number of NTUs 26
# ntaxa(euk.RNA) - ntaxa(mf.RNA) #so number of bacterial NTUs only? 4151
# mf.RNA.tax = data.frame(tax_table(mf.RNA)) #microfauna taxonomy



#community structure

#function for pulling out phylum-level relative abundance
com_structure <- function(ps, Habitat){
  
  sub <- prune_samples(Habitat, ps)
  glom <- tax_glom(sub, taxrank = 'Phylum')
  t.glom = data.frame(otu_table(glom))
  rownames(t.glom) = data.frame(tax_table(glom))$Phylum
  x = sort(rowSums(t.glom))/sum(t.glom)*100
  
  return(x)
  
}

#16S Relative abundances

#16S Snow DNA
pro.snow.dna = com_structure(pro.DNA, data.frame(sample_data(subset_samples(pro.DNA, Habitat2 == "Snow")))$SampleID)
pro.snow.dna
#16S Snow RNA
pro.snow.rna = com_structure(pro.RNA, data.frame(sample_data(subset_samples(pro.RNA, Habitat2 == "Snow")))$SampleID)
pro.snow.rna
#16S Spring Ice DNA
pro.sp.dna = com_structure(pro.DNA, data.frame(sample_data(subset_samples(pro.DNA, Habitat2 == "Spring Ice")))$SampleID)
pro.sp.dna
#16S Spring Ice RNA
pro.sp.rna = com_structure(pro.RNA, data.frame(sample_data(subset_samples(pro.RNA, Habitat2 == "Spring Ice")))$SampleID)
pro.sp.rna
#16S Summer Ice DNA
pro.sm.dna = com_structure(pro.DNA, data.frame(sample_data(subset_samples(pro.DNA, Habitat2 == "Summer Ice")))$SampleID)
pro.sm.dna
#16S Summer Ice RNA
pro.sm.rna = com_structure(pro.RNA, data.frame(sample_data(subset_samples(pro.RNA, Habitat2 == "Summer Ice")))$SampleID)
pro.sm.rna
#16S Cryoconite DNA
pro.cry.dna = com_structure(pro.DNA, data.frame(sample_data(subset_samples(pro.DNA, Habitat2 == "Cryoconite")))$SampleID)
pro.cry.dna
#16S Cryoconite RNA
pro.cry.rna = com_structure(pro.RNA, data.frame(sample_data(subset_samples(pro.RNA, Habitat2 == "Cryoconite")))$SampleID)
pro.cry.rna


#18S Relative abundances

#remove microfauna from calculcating relative abundance
euk.DNA.no.mf = subset_taxa(euk.DNA, !Phylum %in% c("Tardigrada", "Rotifera", "Nematozoa"))

#18S Snow DNA
euk.snow.dna = com_structure(euk.DNA.no.mf, data.frame(sample_data(subset_samples(euk.DNA.no.mf, Habitat2 == "Snow")))$SampleID)
euk.snow.dna
#18S Snow RNA
euk.snow.rna = com_structure(euk.RNA, data.frame(sample_data(subset_samples(euk.RNA, Habitat2 == "Snow")))$SampleID)
euk.snow.rna
#18S Spring Ice DNA
euk.sp.dna = com_structure(euk.DNA.no.mf, data.frame(sample_data(subset_samples(euk.DNA.no.mf, Habitat2 == "Spring Ice")))$SampleID)
euk.sp.dna
#18S Spring Ice RNA
euk.sp.rna = com_structure(euk.RNA, data.frame(sample_data(subset_samples(euk.RNA, Habitat2 == "Spring Ice")))$SampleID)
euk.sp.rna
#18S Summer Ice DNA
euk.sm.dna = com_structure(euk.DNA.no.mf, data.frame(sample_data(subset_samples(euk.DNA.no.mf, Habitat2 == "Summer Ice")))$SampleID)
euk.sm.dna
#18S Summer Ice RNA
euk.sm.rna = com_structure(euk.RNA, data.frame(sample_data(subset_samples(euk.RNA, Habitat2 == "Summer Ice")))$SampleID)
euk.sm.rna
#18S Cryoconite DNA
euk.cry.dna = com_structure(euk.DNA.no.mf, data.frame(sample_data(subset_samples(euk.DNA.no.mf, Habitat2 == "Cryoconite")))$SampleID)
euk.cry.dna
#18S Cryoconite RNA
euk.cry.rna = com_structure(euk.RNA, data.frame(sample_data(subset_samples(euk.RNA, Habitat2 == "Cryoconite")))$SampleID)
euk.cry.rna
