#Analysis of microbial communities obtained from total RNA using phyloflash analysis

#clean workspace and load package
rm(list=ls())
graphics.off()

# library(data.table) #for list.files
library(dplyr) #for the infix operator %>%
library(purrr) #for map_df and reduce functions
library(tidyr) #for separate function
library(stringr) #for string_replace_all 
library(phyloseq) #for phyloseq functions
library(ggplot2)
library(tibble) #rows_to_columns

#Load metadata
metadata <- read.csv(file="../data/mt-metadata.csv", sep=",") #Metadata
rownames(metadata) = metadata$SampleID

#remove amplicon field blanks from metadata
metadata = metadata[!(metadata$SampleID %in% c("Blank4", "Blank5", "Blank6", "Blank7", "Blank8")),]

#import all files
list_csv_files =   list.files(path = "../data/mt/", pattern = "*.csv") 
#ignore samples S21.28 and S21.88 because we did not recover any reads from them.
list_csv_files = list_csv_files[!list_csv_files %in% c("S21-28.phyloFlash.NTUfull_abundance.csv", "S21-88.phyloFlash.NTUfull_abundance.csv")]
#put all our csv files into a list
df = lapply(list_csv_files, function(x) read.csv(paste0("../data/mt/", x), header=FALSE))

#start by merging our dfs
df.join <- df %>%
  reduce(full_join, by="V1")

#make sure they have the correct sample names
sam.nam = unlist(lapply(list_csv_files, function(x) paste0("S21.", substr(x, 5, 6))))
names(df.join) = c("Taxon", sam.nam)


NTUtable = df.join

#create tax-table
path_split <- strsplit(NTUtable[,1], ";")

#tax_slv_138.1.txt downloaded from https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/
silva <- read.table("../data/tax_slv_ssu_138.1.txt", h = F, sep = "\t", stringsAsFactors = F) # import taxa map from version used to annotate with phyloFLash

#silva_map code and SILVAtaxopath function taken from https://github.com/ChrisTrivedi/Trivedi_et_al_DNA_RNA_Preservation/blob/main/TotalRNA_phyloFlash_NTUabundance_to_phyloseq.R
silva_map <- data.frame(  # prepare taxa map in right format for parsing function
  path = gsub(";$", "", silva$V1),
  node = sapply(strsplit(silva$V1, ";"), function(x) x[length(x)]),
  rank = silva$V3,
  stringsAsFactors = T
)

# Create new function for SILVA taxonomy
SILVAtaxopath <- function(tax, SILVA){    # parsing function provided by Christiane Hassenrück @chassenr on github
  output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank)))
  colnames(output) <- levels(SILVA$rank)
  for (i in 1:length(tax)) {
    for (j in 1:length(levels(SILVA$rank))) {
      if (paste(tax[[i]][1:j], collapse = ";") %in% SILVA$path) {
        output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
      }
    }
  }
  return(output)
}

# Use function to organize taxonomic calls
TAXmat <- SILVAtaxopath(path_split,silva_map)
prefix<- "NTU"  #create rownames corresponding to NTUs
suffix<- seq(1:nrow(NTUtable))
NTU.names<- paste(prefix,suffix, sep = "")
row.names(TAXmat)<-NTU.names
TAXmat <- TAXmat[,c("domain","major_clade","kingdom","phylum","class","order","family","genus")]
#change column names to first letter upper case
colnames(TAXmat) = str_to_title(names(as.data.frame(TAXmat)))
#remove Major clade and kingdom columns
TAXmat<- TAXmat[,colnames(TAXmat) != "Major_clade"]
TAXmat<- TAXmat[,colnames(TAXmat) != "Kingdom"]

#We have a lot of NA in the phylum column, even though we can identify the organism at lower tax levels
#I'm going at manually add in the phylums for plotting
tax.df = as.data.frame(TAXmat)
tax.check = tax.df[is.na(tax.df$Phylum),]
sort(table(tax.check$Class)) 

#import the taxonomy table from SILVA
tax <- read_tsv(file="../data/taxonomy.tsv") #Taxonomy table

#format
tax <- tax %>%
  #mutate(tax=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(tax=str_replace_all(string=Taxon, pattern=".__", replacement="")) %>%
  mutate(tax=str_replace_all(string=tax, pattern=";$", replacement="")) %>%
  mutate(tax=str_replace_all(string=tax, pattern=" ", replacement="")) %>%
  separate(tax, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";") %>%
  dplyr::select(-Taxon) %>%
  column_to_rownames(var='Feature ID')

#only keep unique phylum and class pairings
check.list = unique(tax[,c('Phylum', 'Class')])

#based on these classifications, fill in the missing phyla for our phyloflash data
for (i in 1:nrow(tax.df)){
  if (is.na(tax.df[i,]$Phylum)){ #if the phylum is na
    if (!is.na(tax.df[i,]$Class)){ #is the class is not na
      if (tax.df[i,]$Class %in% check.list$Class){ #if the class is in the amplicon list
        if (tax.df[i,]$Class != "uncultured"){ #if the class is not "uncultured"
        tax.df[i,]$Phylum = check.list[which(check.list$Class == tax.df[i,]$Class),]$Phylum #replace the phylum according to the taxonomy of our classes in amplicon data
        }
      }
      
    }
  }
}

#check how many ASVs are still NA at the phylum level
tax.check2 = tax.df[is.na(tax.df$Phylum),]
sort(table(tax.check2$Class)) 

#remove taxonomy from df.join and add NTU labels
df.join = df.join[,-c(1)]
rownames(df.join) = rownames(tax.df)

#change names of blanks
names(df.join)[(ncol(df.join)-2):(ncol(df.join))] = c("Blank1", "Blank2", "Blank3")

#turn count table NA into zero
df.join[is.na(df.join)] = 0

#tax.df to matrix
tax.df = as.matrix(tax.df)

#make into phyloseq objects
ps.pro = phyloseq(otu_table(as.matrix(df.join), taxa_are_rows = TRUE), tax_table(tax.df), sample_data(metadata))
ps.euk = phyloseq(otu_table(as.matrix(df.join), taxa_are_rows = TRUE), tax_table(tax.df), sample_data(metadata))
ps.pro
ps.euk

pro.edit <- subset_taxa(ps.pro, !is.na(Domain) & !Domain %in% c("Unassigned", "Eukaryota") & !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria"))
pro.edit <- prune_taxa(taxa_sums(pro.edit) >= 1, pro.edit)
pro.edit

#eukaryote dataset with microfauna removed.
euk.edit <- subset_taxa(ps.euk, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !is.na(Phylum) & !Phylum %in% c("Vertebrata", "Tardigrada", "Rotifera", "Nematozoa", "Arthropoda") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
euk.edit <- prune_taxa(taxa_sums(euk.edit) >= 1, euk.edit) 
euk.edit

#total eukaryote dataset this is the same because we had not microfauna
euk.edit.tot <- subset_taxa(ps.euk, !is.na(Domain) & !Domain %in% c("Unassigned", "Bacteria", "Archaea") & !is.na(Phylum) & !Phylum %in% c("Vertebrata", "Arthropoda") & !Class %in% c("Embryophyta") & !Family %in% c("Mitochondria"))
euk.edit.tot <- prune_taxa(taxa_sums(euk.edit.tot) >= 1, euk.edit.tot) 
euk.edit.tot

#export phyloseq objects
saveRDS(pro.edit, "../results/mt-rna-16S-phylo-object.rds")
saveRDS(euk.edit, "../results/mt-rna-18S-phylo-object.rds")
