#Community stats

#How many reads in each dataset?
#How many ASVs/NTUs in each dataset?
#Relative abundance of phyla?
#Relative abundance of class?
#Relative abundance of orders?
#Relative abundance of families?
#Relative abundance of genera?

#Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)

#Import data

#prokaryotes
pro.DNA <- readRDS("../results/mt-dna-16S-phylo-object-no-controls.rds") 
pro.RNA <- readRDS("../results/mt-rna-16S-phylo-object-no-controls.rds") 

#eukaryotes without micrometazoans
euk.DNA <- readRDS("../results/mt-dna-18S-phylo-object-no-controls.rds") 
euk.RNA <- readRDS("../results/mt-rna-18S-phylo-object-no-controls.rds") 

#modify taxonomic names

#ps = pro.DNA
name_mod <- function(ps){
  
  ps.m = ps
  x = data.frame(tax_table(ps.m)) %>% 
    mutate(Phylum = ifelse(is.na(Phylum), paste0(Domain, " P.",Phylum), Phylum)) %>%
    mutate(Class = ifelse(is.na(Class), paste0(Phylum, " C.", Class), Class)) %>%
    mutate(Class = ifelse(Class == "Cyanobacteriia", "Cyanophyceae", Class)) %>%
    mutate(Class = ifelse(Class == "uncultured", paste0(Phylum, " uncult."), Class)) %>%
    mutate(Class = ifelse(Class == "Incertae_Sedis", paste0(Phylum, " (incertae sedis)"), Class)) %>%
    mutate(Class = ifelse(Class == "Incertae Sedis", paste0(Phylum, " (incertae sedis)"), Class)) %>%
    mutate(Order = ifelse(is.na(Order), paste0(Class, " O.", Order), Order)) %>%
    mutate(Order = ifelse(Order == "Leptolyngbyales", "Synechococcales", Order)) %>%
    mutate(Order = ifelse(Order == "uncultured", paste0(Class, " uncult."), Order)) %>%
    mutate(Order = ifelse(Order == "Incertae_Sedis", paste0(Class, " (incertae sedis)"), Order)) %>%
    mutate(Order = ifelse(Order == "Incertae Sedis", paste0(Class, " (incertae sedis)"), Order)) %>%
    mutate(Family = ifelse(is.na(Family), paste0(Order, " F.", Family), Family)) %>%
    mutate(Family = ifelse(Family == "uncultured", paste0(Order, " uncult."), Family)) %>%
    mutate(Family = ifelse(Family == "Unknown_Family", paste0(Order, " (family unknown)"), Family)) %>%
    mutate(Family = ifelse(Family == "Unknown Family", paste0(Order, " (family unknown)"), Family)) %>%
    mutate(Family = ifelse(Family == "Incertae_Sedis", paste0(Order, " (incertae sedis)"), Family)) %>%
    mutate(Family = ifelse(Family == "Incertae Sedis", paste0(Order, " (incertae sedis)"), Family)) %>%
    mutate(Family = ifelse(Family == "Pleosporales", "Pleosporales (Genus Unknown)", Family)) %>%
    mutate(Family = ifelse(Family == "Rhizophydiales", "Rhizophydiales (Genus Unknown)", Family)) %>%
    mutate(Family = ifelse(Family == "Thecofilosea", "Thecofilosea (Genus Unknown)", Family)) %>%
    mutate(Genus = ifelse(is.na(Genus), paste0(Family, " G.", Genus), Genus)) %>%
    mutate(Genus = ifelse(Genus == "uncultured", paste0(Family, " uncult."), Genus)) %>%
    mutate(Genus = ifelse(Genus == "Unknown_Family", paste0(Order, " (family and genus unknown)"), Genus)) %>%
    mutate(Phylum = ifelse(Class == "Labyrinthulomycetes", "Bigyra", Phylum)) %>%
    mutate(Phylum = ifelse(Phylum == "Labyrinthulomycetes", "Bigyra", Phylum)) %>%
    mutate(Class = ifelse(Class == "Labyrinthulomycetes (Genus Unknown)", "Labyrinthulomycetes", Class)) %>%
    mutate(Class = ifelse(Order == "Burkholderiales", "Betaproteobacteria", Class)) %>%
    mutate(Phylum = ifelse(Phylum == "Actinobacteria", "Actinobacteriota", Phylum)) %>%
    mutate(Class = ifelse(Class == "Actinobacteria", "Actinomycetia", Class)) %>%
    mutate(Class = ifelse(Genus == "Ferruginibacter", "Chitinophagia", Class)) %>%
    mutate(Class = ifelse(Genus == "Hymenobacter", "Cytophagia", Class))%>%
    mutate(Class = ifelse(Genus == "Arcicella", "Cytophagia", Class))%>%
    mutate(Class = ifelse(Order == "Cytophagales", "Cytophagia", Class)) %>%
    mutate(Class = ifelse(Order == "Flavobacteriales", "Flavobacteriia", Class)) %>%
    mutate(Class = ifelse(Order == "Sphingobacteriales", "Sphingobacteriia", Class)) %>%
    mutate(Class = ifelse(Order == "Chitinophagales", "Chitinophagia", Class))%>%
    mutate(Order = ifelse(Genus == "Nostoc_PCC-73102", "Nostocales", Order)) %>%
    mutate(Order = ifelse(Family == "Phormidiaceae", "Oscillatoriales", Order)) %>%
    mutate(Order = ifelse(Family == "Nostocaceae", "Nostocales", Order))
    
    

  
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
  y$Class <- gsub("Labyrinthulomycetes (Genus Unknown)", "Labyrinthulomycetes", y$Class)
  y$Class <- ifelse(y$Class %in% c("Labyrinthulomycetes (Genus Unknown)"), "Labyrinthulomycetes", y$Class)
  y$Class <- ifelse(y$Class %in% c("Labyrinthulomycetes uncult."), "Labyrinthulomycetes", y$Class)
  y$Order <- ifelse(y$Order %in% c("Thecofilosea (Genus Unknown)"), "Thecofilosea", y$Order)
  
  
  tax_table(ps.m) = as.matrix(y)
  
  return(ps.m)
}

#get community structure

com_structure <- function(ps1, ps2, rank){
  
  habs = c("Snow", "Spring Ice", "Summer Ice", "Cryoconite")
  
  sample_list = list(sam1a = data.frame(sample_data(subset_samples(ps1, Habitat2 == "Snow")))$SampleID,
  sam2a = data.frame(sample_data(subset_samples(ps1, Habitat2 == "Spring Ice")))$SampleID,
  sam3a = data.frame(sample_data(subset_samples(ps1, Habitat2 == "Summer Ice")))$SampleID,
  sam4a = data.frame(sample_data(subset_samples(ps1, Habitat2 == "Cryoconite")))$SampleID,
  sam1b = data.frame(sample_data(subset_samples(ps2, Habitat2 == "Snow")))$SampleID,
  sam2b = data.frame(sample_data(subset_samples(ps2, Habitat2 == "Spring Ice")))$SampleID,
  sam3b = data.frame(sample_data(subset_samples(ps2, Habitat2 == "Summer Ice")))$SampleID,
  sam4b = data.frame(sample_data(subset_samples(ps2, Habitat2 == "Cryoconite")))$SampleID)
  
  df2save = data.frame(Group = as.character())
  
  for (i in 1:length(habs)){
    
    sub <- prune_samples(sample_list[[i]], ps1)
    glom <- tax_glom(sub, taxrank = rank)
    t.glom = data.frame(otu_table(glom))
    df = data.frame(tax_table(glom))
    row.nam = unlist(as.vector(df[names(df) %in% rank]))
    row.nam[duplicated(row.nam)]
    rownames(t.glom) = row.nam
    x = data.frame(sort(rowSums(t.glom), decreasing = TRUE)/sum(t.glom)*100)
    x = cbind(rownames(x), x)
    names(x) = c("Group", "Relative_Abundance")
    x$Relative_Abundance = round(x$Relative_Abundance, 2)
    x = x[x$Relative_Abundance > 0.01,]
    rownames(x) <- NULL
    

    #sub2 <- subset_samples(ps2, Habitat2 == habs[[i]])
    sub2 <- prune_samples(sample_list[[i+4]], ps2)
    glom2 <- tax_glom(sub2, taxrank = rank)
    t.glom2 = data.frame(otu_table(glom2))
    df2 = data.frame(tax_table(glom2))
    row.nam2 = unlist(as.vector(df2[names(df2) %in% rank]))
    row.nam2[duplicated(row.nam2)]
    rownames(t.glom2) = row.nam2
    y = data.frame(sort(rowSums(t.glom2), decreasing = TRUE)/sum(t.glom2)*100)
    y = cbind(rownames(y), y)
    names(y) = c("Group", "Relative_Abundance")
    y$Relative_Abundance = round(y$Relative_Abundance, 2)
    y = y[y$Relative_Abundance > 0.01,]
    rownames(y) <- NULL
    
    #join the datafames
    z = full_join(x, y, by="Group")
    names(z) = c("Group", paste0(habs[[i]], " DNA RA"), paste0(habs[[i]], " RNA RA"))
    
    df2save = full_join(df2save, z, by="Group")
    
  }

  
  return(df2save)
  
}

#run function
pro.dom = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Domain")
pro.phy = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Phylum")
pro.cla = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Class")
pro.ord = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Order")
# pro.fam = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Family")
# pro.gen = com_structure(ps1 = name_mod(pro.DNA), ps2 = name_mod(pro.RNA), rank = "Genus")

euk.dom = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Domain")
euk.phy = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Phylum")
euk.cla = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Class")
euk.ord = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Order")
# euk.fam = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Family")
# euk.gen = com_structure(ps1 = name_mod(euk.DNA), ps2 = name_mod(euk.RNA), rank = "Genus")

#most abundant genera per phyla

ps = name_mod(euk.RNA)
rank = "Class"
phyla = c("Ascomycota", "Cercozoa", "Chlorophyta")

ps = name_mod(sub)
rank = "Class"
phyla = c("Phragmoplastophyta", "Amoebozoa")

phy_structure <- function(ps, rank, phyla){
  
  #phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria")
  
  res_list = list()
  
  for (i in 1:length(phyla)){
    
    check.data = data.frame(tax_table(ps))
    ps.phy = prune_taxa(data.frame(tax_table(ps))$Phylum %in% phyla[[i]], ps)
    check.data = data.frame(tax_table(ps.phy))
    glom <- tax_glom(ps.phy, taxrank = rank)
    t.glom = data.frame(otu_table(glom))
    df = data.frame(tax_table(glom))
    row.nam = unlist(as.vector(df[names(df) %in% rank]))
    row.nam[duplicated(row.nam)]
    rownames(t.glom) = row.nam
    x = data.frame(sort(rowSums(t.glom), decreasing = TRUE)/sum(t.glom)*100)
    x = cbind(rownames(x), x)
    names(x) = c("Group", "Relative_Abundance")
    x$Relative_Abundance = round(x$Relative_Abundance, 2)
    x = x[x$Relative_Abundance > 1,]
    rownames(x) <- NULL
    
    res_list[[i]] = x
  }
  
  return(res_list)
}

#Community structure numbers

####################################################################################################################

#output results to text document
sink("../results/mt-community-profiles.txt", type="output")
writeLines("===============================================================
COMMUNITY PROFILES
===============================================================")

writeLines("Total DNA prokaryote ASVs:")
ntaxa(pro.DNA)
writeLines("Total RNA prokaryote NTUs:")
ntaxa(pro.RNA)
writeLines("Total DNA microbial eukaryote ASVs:")
ntaxa(euk.DNA)
writeLines("Total RNA microbial eukaryote NTUs:")
ntaxa(euk.RNA)

writeLines("Total DNA prokaryote reads:")
sum(taxa_sums(pro.DNA))
writeLines("Total RNA prokaryote reads:")
sum(taxa_sums(pro.RNA))
writeLines("Total DNA microbial eukaryote reads:")
sum(taxa_sums(euk.DNA))
writeLines("Total RNA microbial eukaryote reads:")
sum(taxa_sums(euk.RNA))

writeLines("Number of DNA Archaea ASVs:")
arc.DNA <- subset_taxa(pro.DNA, Domain %in% c("Archaea")) #of which are archaea?
arc.DNA = prune_samples(sample_sums(arc.DNA) >0, arc.DNA) #remove samples with zero counts
arc.DNA = filter_taxa(arc.DNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
ntaxa(arc.DNA)
writeLines("Taxonomy of DNA Archaea ASVs:")
data.frame(tax_table(arc.DNA))[1:5] 

writeLines("Number of RNA Archaea NTUs:")
arc.RNA <- subset_taxa(pro.RNA, Domain %in% c("Archaea")) #of which are archaea?
arc.RNA = prune_samples(sample_sums(arc.RNA) >0, arc.RNA) #remove samples with zero counts
arc.RNA = filter_taxa(arc.RNA, function(x) sum(x) > 0, TRUE) #remove taxa with zero counts
ntaxa(arc.RNA)
writeLines("Taxonomy of RNA Archaea NTUs:")
data.frame(tax_table(arc.RNA))[1:5] 

writeLines("Prokaryote Domain:")
pro.dom

writeLines("Prokaryote Phyla:")
pro.phy
writeLines("Eukaryote Phyla:")
euk.phy

writeLines("Prokaryote Class:")
pro.cla
writeLines("Eukaryote Class:")
euk.cla

writeLines("Prokaryote Order:")
pro.ord
writeLines("Eukaryote Order:")
euk.ord


pro_groups_cla <- phy_structure(name_mod(pro.DNA), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
writeLines("Proteobacteria Classes:")
pro_groups_cla[[1]]
writeLines("Actinobacteria Classes:")
pro_groups_cla[[2]]
writeLines("Bacteroidota Classes:")
pro_groups_cla[[3]]
writeLines("Cyanobacteria Classes:")
pro_groups_cla[[4]]

pro_groups_ord <- phy_structure(name_mod(pro.DNA), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
writeLines("Proteobacteria Orders:")
pro_groups_ord[[1]]
writeLines("Actinobacteria Orders:")
pro_groups_ord[[2]]
writeLines("Bacteroidota Orders:")
pro_groups_ord[[3]]
writeLines("Cyanobacteria Orders:")
pro_groups_ord[[4]]

pro_groups_fam <- phy_structure(name_mod(pro.DNA), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
writeLines("Proteobacteria Families:")
pro_groups_fam[[1]]
writeLines("Actinobacteria Families:")
pro_groups_fam[[2]]
writeLines("Bacteroidota Families:")
pro_groups_fam[[3]]
writeLines("Cyanobacteria Families:")
pro_groups_fam[[4]]

pro_groups_gen <- phy_structure(name_mod(pro.DNA), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
writeLines("Proteobacteria Genera:")
pro_groups_gen[[1]]
writeLines("Actinobacteria Genera:")
pro_groups_gen[[2]]
writeLines("Bacteroidota Genera:")
pro_groups_gen[[3]]
writeLines("Cyanobacteria Genera:")
pro_groups_gen[[4]]

euk_groups_cla <- phy_structure(name_mod(euk.DNA), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
writeLines("Ascomycota Classes:")
euk_groups_cla[[1]]
writeLines("Cercozoa Classes:")
euk_groups_cla[[2]]
writeLines("Chlorophyta Classes:")
euk_groups_cla[[3]]

euk_groups_ord <- phy_structure(name_mod(euk.DNA), "Order", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
writeLines("Ascomycota Orders:")
euk_groups_ord[[1]]
writeLines("Cercozoa Orders:")
euk_groups_ord[[2]]
euk_groups_ord[[3]]

euk_groups_fam <- phy_structure(name_mod(euk.DNA), "Family", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
writeLines("Ascomycota Families:")
euk_groups_fam[[1]]
writeLines("Cercozoa Families:")
euk_groups_fam[[2]]
writeLines("Chlorophyta Families:")
euk_groups_fam[[3]]

euk_groups_gen <- phy_structure(name_mod(euk.DNA), "Genus", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
writeLines("Ascomycota Genera:")
euk_groups_gen[[1]]
writeLines("Cercozoa Genera:")
euk_groups_gen[[2]]
writeLines("Chlorophyta Genera:")
euk_groups_gen[[3]]

sink()

#####################################################################################
#BREAKDOWN OF TAXONOMY BY HABITAT
#PROKARYOTES DNA CLASSES
sub = subset_samples(pro.DNA, Habitat2 == "Snow")
pro_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Spring Ice")
pro_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Summer Ice")
pro_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Cryoconite")
pro_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Classes
x = full_join(pro_groups_cla_sn[[1]], pro_groups_cla_sp[[1]], by="Group")
x = full_join(x, pro_groups_cla_sm[[1]], by="Group")
x = full_join(x, pro_groups_cla_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Classes
x2 = full_join(pro_groups_cla_sn[[2]], pro_groups_cla_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_cla_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_cla_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Classes
x3 = full_join(pro_groups_cla_sn[[3]], pro_groups_cla_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_cla_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_cla_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Classes
x4 = full_join(pro_groups_cla_sn[[4]], pro_groups_cla_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_cla_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_cla_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df = rbind(x, x2, x3, x4)
write_csv(df, "../results/mt-dna-pro-key-phyla-class-habitat-profiles.csv")

#PROKARYOTES DNA ORDERS
sub = subset_samples(pro.DNA, Habitat2 == "Snow")
pro_groups_ord_sn <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Spring Ice")
pro_groups_ord_sp <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Summer Ice")
pro_groups_ord_sm <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Cryoconite")
pro_groups_ord_cr <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Order
x = full_join(pro_groups_ord_sn[[1]], pro_groups_ord_sp[[1]], by="Group")
x = full_join(x, pro_groups_ord_sm[[1]], by="Group")
x = full_join(x, pro_groups_ord_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Order
x2 = full_join(pro_groups_ord_sn[[2]], pro_groups_ord_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_ord_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_ord_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Order
x3 = full_join(pro_groups_ord_sn[[3]], pro_groups_ord_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_ord_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_ord_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Order
x4 = full_join(pro_groups_ord_sn[[4]], pro_groups_ord_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_ord_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_ord_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df2 = rbind(x, x2, x3, x4)
write_csv(df2, "../results/mt-dna-pro-key-phyla-order-habitat-profiles.csv")

#PROKARYOTES DNA FAMILY
sub = subset_samples(pro.DNA, Habitat2 == "Snow")
pro_groups_fam_sn <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Spring Ice")
pro_groups_fam_sp <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Summer Ice")
pro_groups_fam_sm <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Cryoconite")
pro_groups_fam_cr <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Family
x = full_join(pro_groups_fam_sn[[1]], pro_groups_fam_sp[[1]], by="Group")
x = full_join(x, pro_groups_fam_sm[[1]], by="Group")
x = full_join(x, pro_groups_fam_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Family
x2 = full_join(pro_groups_fam_sn[[2]], pro_groups_fam_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_fam_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_fam_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Family
x3 = full_join(pro_groups_fam_sn[[3]], pro_groups_fam_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_fam_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_fam_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Family
x4 = full_join(pro_groups_fam_sn[[4]], pro_groups_fam_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_fam_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_fam_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df3 = rbind(x, x2, x3, x4)
write_csv(df3, "../results/mt-dna-pro-key-phyla-family-habitat-profiles.csv")

#PROKARYOTES DNA GENUS
sub = subset_samples(pro.DNA, Habitat2 == "Snow")
pro_groups_gen_sn <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Spring Ice")
pro_groups_gen_sp <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Summer Ice")
pro_groups_gen_sm <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.DNA, Habitat2 == "Cryoconite")
pro_groups_gen_cr <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Genus
x = full_join(pro_groups_gen_sn[[1]], pro_groups_gen_sp[[1]], by="Group")
x = full_join(x, pro_groups_gen_sm[[1]], by="Group")
x = full_join(x, pro_groups_gen_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Genus
x2 = full_join(pro_groups_gen_sn[[2]], pro_groups_gen_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_gen_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_gen_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Genus
x3 = full_join(pro_groups_gen_sn[[3]], pro_groups_gen_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_gen_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_gen_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Genus
x4 = full_join(pro_groups_gen_sn[[4]], pro_groups_gen_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_gen_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_gen_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df4 = rbind(x, x2, x3, x4)
write_csv(df4, "../results/mt-dna-pro-key-phyla-genera-habitat-profiles.csv")


#PROKARYOTES RNA CLASSES
sub = subset_samples(pro.RNA, Habitat2 == "Snow")
pro_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Spring Ice")
pro_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Summer Ice")
pro_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Cryoconite")
pro_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Classes
x = full_join(pro_groups_cla_sn[[1]], pro_groups_cla_sp[[1]], by="Group")
x = full_join(x, pro_groups_cla_sm[[1]], by="Group")
x = full_join(x, pro_groups_cla_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Classes
x2 = full_join(pro_groups_cla_sn[[2]], pro_groups_cla_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_cla_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_cla_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Classes
x3 = full_join(pro_groups_cla_sn[[3]], pro_groups_cla_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_cla_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_cla_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Classes
x4 = full_join(pro_groups_cla_sn[[4]], pro_groups_cla_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_cla_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_cla_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df5 = rbind(x, x2, x3, x4)
write_csv(df5, "../results/mt-key-pro-rna-phyla-class-habitat-profiles.csv")

#PROKARYOTES RNA ORDERS
sub = subset_samples(pro.RNA, Habitat2 == "Snow")
pro_groups_ord_sn <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Spring Ice")
pro_groups_ord_sp <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Summer Ice")
pro_groups_ord_sm <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Cryoconite")
pro_groups_ord_cr <- phy_structure(name_mod(sub), "Order", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Order
x = full_join(pro_groups_ord_sn[[1]], pro_groups_ord_sp[[1]], by="Group")
x = full_join(x, pro_groups_ord_sm[[1]], by="Group")
x = full_join(x, pro_groups_ord_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Order
x2 = full_join(pro_groups_ord_sn[[2]], pro_groups_ord_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_ord_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_ord_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Order
x3 = full_join(pro_groups_ord_sn[[3]], pro_groups_ord_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_ord_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_ord_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Order
x4 = full_join(pro_groups_ord_sn[[4]], pro_groups_ord_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_ord_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_ord_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df6 = rbind(x, x2, x3, x4)
write_csv(df6, "../results/mt-rna-pro-key-phyla-order-habitat-profiles.csv")

#PROKARYOTES RNA FAMILY
sub = subset_samples(pro.RNA, Habitat2 == "Snow")
pro_groups_fam_sn <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Spring Ice")
pro_groups_fam_sp <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Summer Ice")
pro_groups_fam_sm <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Cryoconite")
pro_groups_fam_cr <- phy_structure(name_mod(sub), "Family", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Family
x = full_join(pro_groups_fam_sn[[1]], pro_groups_fam_sp[[1]], by="Group")
x = full_join(x, pro_groups_fam_sm[[1]], by="Group")
x = full_join(x, pro_groups_fam_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Family
x2 = full_join(pro_groups_fam_sn[[2]], pro_groups_fam_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_fam_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_fam_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Family
x3 = full_join(pro_groups_fam_sn[[3]], pro_groups_fam_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_fam_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_fam_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Family
x4 = full_join(pro_groups_fam_sn[[4]], pro_groups_fam_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_fam_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_fam_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df7 = rbind(x, x2, x3, x4)
write_csv(df7, "../results/mt-rna-pro-key-phyla-family-habitat-profiles.csv")

#PROKARYOTES RNA GENUS
sub = subset_samples(pro.RNA, Habitat2 == "Snow")
pro_groups_gen_sn <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Spring Ice")
pro_groups_gen_sp <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Summer Ice")
pro_groups_gen_sm <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))
sub = subset_samples(pro.RNA, Habitat2 == "Cryoconite")
pro_groups_gen_cr <- phy_structure(name_mod(sub), "Genus", phyla = c("Proteobacteria", "Actinobacteriota", "Bacteroidota", "Cyanobacteria"))

#Proteobacteria Genus
x = full_join(pro_groups_gen_sn[[1]], pro_groups_gen_sp[[1]], by="Group")
x = full_join(x, pro_groups_gen_sm[[1]], by="Group")
x = full_join(x, pro_groups_gen_cr[[1]], by="Group")
x = cbind(Taxa = "Proteobacteria", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Actinobacteria Genus
x2 = full_join(pro_groups_gen_sn[[2]], pro_groups_gen_sp[[2]], by="Group")
x2 = full_join(x2, pro_groups_gen_sm[[2]], by="Group")
x2 = full_join(x2, pro_groups_gen_cr[[2]], by="Group")
x2 = cbind(Taxa = "Actinobacteria", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Bacteroidota Genus
x3 = full_join(pro_groups_gen_sn[[3]], pro_groups_gen_sp[[3]], by="Group")
x3 = full_join(x3, pro_groups_gen_sm[[3]], by="Group")
x3 = full_join(x3, pro_groups_gen_cr[[3]], by="Group")
x3 = cbind(Taxa = "Bacteroidota", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cyanobacteria Genus
x4 = full_join(pro_groups_gen_sn[[4]], pro_groups_gen_sp[[4]], by="Group")
x4 = full_join(x4, pro_groups_gen_sm[[4]], by="Group")
x4 = full_join(x4, pro_groups_gen_cr[[4]], by="Group")
x4 = cbind(Taxa = "Cyanobacteria", x4)
names(x4) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df8 = rbind(x, x2, x3, x4)
write_csv(df8, "../results/mt-rna-pro-key-phyla-genera-habitat-profiles.csv")


#######EUKARYOTES############

#EUKARYOTES DNA CLASSES
sub = subset_samples(euk.DNA, Habitat2 == "Snow")
euk_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Spring Ice")
euk_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Summer Ice")
euk_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Cryoconite")
euk_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))

#Ascomycota Classes
x = full_join(euk_groups_cla_sn[[1]], euk_groups_cla_sp[[1]], by="Group")
x = full_join(x, euk_groups_cla_sm[[1]], by="Group")
x = full_join(x, euk_groups_cla_cr[[1]], by="Group")
x = cbind(Taxa = "Ascomycota", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Classes
x2 = full_join(euk_groups_cla_sn[[2]], euk_groups_cla_sp[[2]], by="Group")
x2 = full_join(x2, euk_groups_cla_sm[[2]], by="Group")
x2 = full_join(x2, euk_groups_cla_cr[[2]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Classes
x3 = full_join(euk_groups_cla_sn[[3]], euk_groups_cla_sp[[3]], by="Group")
x3 = full_join(x3, euk_groups_cla_sm[[3]], by="Group")
x3 = full_join(x3, euk_groups_cla_cr[[3]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df9 = rbind(x, x2, x3)
write_csv(df9, "../results/mt-dna-euk-phyla-class-habitat-profiles.csv")


#EUKARYOTES DNA ORDER
sub = subset_samples(euk.DNA, Habitat2 == "Snow")
euk_groups_ord_sn <- phy_structure(name_mod(sub), "Order", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Spring Ice")
euk_groups_ord_sp <- phy_structure(name_mod(sub), "Order", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Summer Ice")
euk_groups_ord_sm <- phy_structure(name_mod(sub), "Order", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Cryoconite")
euk_groups_ord_cr <- phy_structure(name_mod(sub), "Order", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))

#Ascomycota Order
x = full_join(euk_groups_ord_sn[[1]], euk_groups_ord_sp[[1]], by="Group")
x = full_join(x, euk_groups_ord_sm[[1]], by="Group")
x = full_join(x, euk_groups_ord_cr[[1]], by="Group")
x = cbind(Taxa = "Ascomycota", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Order
x2 = full_join(euk_groups_ord_sn[[2]], euk_groups_ord_sp[[2]], by="Group")
x2 = full_join(x2, euk_groups_ord_sm[[2]], by="Group")
x2 = full_join(x2, euk_groups_ord_cr[[2]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Order
x3 = full_join(euk_groups_ord_sn[[3]], euk_groups_ord_sp[[3]], by="Group")
x3 = full_join(x3, euk_groups_ord_sm[[3]], by="Group")
x3 = full_join(x3, euk_groups_ord_cr[[3]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df10 = rbind(x, x2, x3)
write_csv(df10, "../results/mt-dna-euk-phyla-order-habitat-profiles.csv")


#EUKARYOTES DNA FAMILY
sub = subset_samples(euk.DNA, Habitat2 == "Snow")
euk_groups_fam_sn <- phy_structure(name_mod(sub), "Family", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Spring Ice")
euk_groups_fam_sp <- phy_structure(name_mod(sub), "Family", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Summer Ice")
euk_groups_fam_sm <- phy_structure(name_mod(sub), "Family", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Cryoconite")
euk_groups_fam_cr <- phy_structure(name_mod(sub), "Family", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))

#Ascomycota Family
x = full_join(euk_groups_fam_sn[[1]], euk_groups_fam_sp[[1]], by="Group")
x = full_join(x, euk_groups_fam_sm[[1]], by="Group")
x = full_join(x, euk_groups_fam_cr[[1]], by="Group")
x = cbind(Taxa = "Ascomycota", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Family
x2 = full_join(euk_groups_fam_sn[[2]], euk_groups_fam_sp[[2]], by="Group")
x2 = full_join(x2, euk_groups_fam_sm[[2]], by="Group")
x2 = full_join(x2, euk_groups_fam_cr[[2]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Family
x3 = full_join(euk_groups_fam_sn[[3]], euk_groups_fam_sp[[3]], by="Group")
x3 = full_join(x3, euk_groups_fam_sm[[3]], by="Group")
x3 = full_join(x3, euk_groups_fam_cr[[3]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df11 = rbind(x, x2, x3)
write_csv(df11, "../results/mt-dna-euk-phyla-family-habitat-profiles.csv")

#EUKARYOTES DNA GENUS
sub = subset_samples(euk.DNA, Habitat2 == "Snow")
euk_groups_gen_sn <- phy_structure(name_mod(sub), "Genus", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Spring Ice")
euk_groups_gen_sp <- phy_structure(name_mod(sub), "Genus", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Summer Ice")
euk_groups_gen_sm <- phy_structure(name_mod(sub), "Genus", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.DNA, Habitat2 == "Cryoconite")
euk_groups_gen_cr <- phy_structure(name_mod(sub), "Genus", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))

#Ascomycota Genus
x = full_join(euk_groups_gen_sn[[1]], euk_groups_gen_sp[[1]], by="Group")
x = full_join(x, euk_groups_gen_sm[[1]], by="Group")
x = full_join(x, euk_groups_gen_cr[[1]], by="Group")
x = cbind(Taxa = "Ascomycota", x)
names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Genus
x2 = full_join(euk_groups_gen_sn[[2]], euk_groups_gen_sp[[2]], by="Group")
x2 = full_join(x2, euk_groups_gen_sm[[2]], by="Group")
x2 = full_join(x2, euk_groups_gen_cr[[2]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Genus
x3 = full_join(euk_groups_gen_sn[[3]], euk_groups_gen_sp[[3]], by="Group")
x3 = full_join(x3, euk_groups_gen_sm[[3]], by="Group")
x3 = full_join(x3, euk_groups_gen_cr[[3]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df12 = rbind(x, x2, x3)
write_csv(df12, "../results/mt-dna-euk-phyla-genus-habitat-profiles.csv")

#EUKARYOTES RNA CLASSES
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
#euk_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
euk_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
#euk_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
euk_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
#euk_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
euk_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
#euk_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Ascomycota", "Cercozoa", "Chlorophyta"))
euk_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Cercozoa", "Chlorophyta"))

#Ascomycota Classes

#######No ascomycota
# x = full_join(euk_groups_cla_sp[[1]], euk_groups_cla_sm[[1]], by="Group")
# #x = full_join(x, euk_groups_cla_sm[[1]], by="Group")
# x = full_join(x, euk_groups_cla_cr[[1]], by="Group")
# x = cbind(Taxa = "Ascomycota",  x[,1], Snow_RA = NA, x[,2:4])
# names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Classes
x2 = full_join(euk_groups_cla_sn[[1]], euk_groups_cla_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_cla_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_cla_cr[[1]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Classes
x3 = full_join(euk_groups_cla_sn[[2]], euk_groups_cla_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_cla_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_cla_cr[[2]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df13 = rbind(x2, x3)
write_csv(df13, "../results/mt-rna-euk-phyla-class-habitat-profiles.csv")


#EUKARYOTES RNA ORDER
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
euk_groups_ord_sn <- phy_structure(name_mod(sub), "Order", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
euk_groups_ord_sp <- phy_structure(name_mod(sub), "Order", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
euk_groups_ord_sm <- phy_structure(name_mod(sub), "Order", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
euk_groups_ord_cr <- phy_structure(name_mod(sub), "Order", phyla = c("Cercozoa", "Chlorophyta"))

#Ascomycota Order
# x = full_join(euk_groups_ord_sn[[1]], euk_groups_ord_sp[[1]], by="Group")
# x = full_join(x, euk_groups_ord_sm[[1]], by="Group")
# x = full_join(x, euk_groups_ord_cr[[1]], by="Group")
# x = cbind(Taxa = "Ascomycota", x)
# names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Order
x2 = full_join(euk_groups_ord_sn[[1]], euk_groups_ord_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_ord_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_ord_cr[[1]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Order
x3 = full_join(euk_groups_ord_sn[[2]], euk_groups_ord_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_ord_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_ord_cr[[2]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df14 = rbind(x2, x3)
write_csv(df14, "../results/mt-rna-euk-phyla-order-habitat-profiles.csv")


#EUKARYOTES RNA FAMILY
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
euk_groups_fam_sn <- phy_structure(name_mod(sub), "Family", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
euk_groups_fam_sp <- phy_structure(name_mod(sub), "Family", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
euk_groups_fam_sm <- phy_structure(name_mod(sub), "Family", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
euk_groups_fam_cr <- phy_structure(name_mod(sub), "Family", phyla = c("Cercozoa", "Chlorophyta"))

#Ascomycota Family
# x = full_join(euk_groups_fam_sn[[1]], euk_groups_fam_sp[[1]], by="Group")
# x = full_join(x, euk_groups_fam_sm[[1]], by="Group")
# x = full_join(x, euk_groups_fam_cr[[1]], by="Group")
# x = cbind(Taxa = "Ascomycota", x)
# names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Family
x2 = full_join(euk_groups_fam_sn[[1]], euk_groups_fam_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_fam_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_fam_cr[[1]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Family
x3 = full_join(euk_groups_fam_sn[[2]], euk_groups_fam_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_fam_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_fam_cr[[2]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df15 = rbind(x2, x3)
write_csv(df15, "../results/mt-rna-euk-phyla-family-habitat-profiles.csv")

#EUKARYOTES RNA GENUS
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
euk_groups_gen_sn <- phy_structure(name_mod(sub), "Genus", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
euk_groups_gen_sp <- phy_structure(name_mod(sub), "Genus", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
euk_groups_gen_sm <- phy_structure(name_mod(sub), "Genus", phyla = c("Cercozoa", "Chlorophyta"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
euk_groups_gen_cr <- phy_structure(name_mod(sub), "Genus", phyla = c("Cercozoa", "Chlorophyta"))

#Ascomycota Genus
# x = full_join(euk_groups_gen_sn[[1]], euk_groups_gen_sp[[1]], by="Group")
# x = full_join(x, euk_groups_gen_sm[[1]], by="Group")
# x = full_join(x, euk_groups_gen_cr[[1]], by="Group")
# x = cbind(Taxa = "Ascomycota", x)
# names(x) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Cercozoa Genus
x2 = full_join(euk_groups_gen_sn[[1]], euk_groups_gen_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_gen_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_gen_cr[[1]], by="Group")
x2 = cbind(Taxa = "Cercozoa", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Chlorophyta Genus
x3 = full_join(euk_groups_gen_sn[[2]], euk_groups_gen_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_gen_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_gen_cr[[2]], by="Group")
x3 = cbind(Taxa = "Chlorophyta", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df16 = rbind(x2, x3)
write_csv(df16, "../results/mt-rna-euk-phyla-genus-habitat-profiles.csv")

#additional investigation
#Which classes make up Phragmoplastophyt and Amoebozoa?
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
euk_groups_cla_sn <- phy_structure(name_mod(sub), "Class", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
euk_groups_cla_sp <- phy_structure(name_mod(sub), "Class", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
euk_groups_cla_sm <- phy_structure(name_mod(sub), "Class", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
euk_groups_cla_cr <- phy_structure(name_mod(sub), "Class", phyla = c("Phragmoplastophyta", "Amoebozoa"))

#Phragmoplastophyta Classes
x2 = full_join(euk_groups_cla_sn[[1]], euk_groups_cla_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_cla_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_cla_cr[[1]], by="Group")
x2 = cbind(Taxa = "Phragmoplastophyta", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Amoebozoa Classes
x3 = full_join(euk_groups_cla_sn[[2]], euk_groups_cla_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_cla_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_cla_cr[[2]], by="Group")
x3 = cbind(Taxa = "Amoebozoa", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df_1 = rbind(x2, x3)
write_csv(df_1, "../results/mt-rna-euk-phyla-class-habitat-profiles-phrag-amoe.csv")

##ORDERS
sub = subset_samples(euk.RNA, Habitat2 == "Snow")
euk_groups_ord_sn <- phy_structure(name_mod(sub), "Order", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Spring Ice")
euk_groups_ord_sp <- phy_structure(name_mod(sub), "Order", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Summer Ice")
euk_groups_ord_sm <- phy_structure(name_mod(sub), "Order", phyla = c("Phragmoplastophyta", "Amoebozoa"))
sub = subset_samples(euk.RNA, Habitat2 == "Cryoconite")
euk_groups_ord_cr <- phy_structure(name_mod(sub), "Order", phyla = c("Phragmoplastophyta", "Amoebozoa"))

#Phragmoplastophyta Order
x2 = full_join(euk_groups_ord_sn[[1]], euk_groups_ord_sp[[1]], by="Group")
x2 = full_join(x2, euk_groups_ord_sm[[1]], by="Group")
x2 = full_join(x2, euk_groups_ord_cr[[1]], by="Group")
x2 = cbind(Taxa = "Phragmoplastophyta", x2)
names(x2) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#Amoebozoa Order
x3 = full_join(euk_groups_ord_sn[[2]], euk_groups_ord_sp[[2]], by="Group")
x3 = full_join(x3, euk_groups_ord_sm[[2]], by="Group")
x3 = full_join(x3, euk_groups_ord_cr[[2]], by="Group")
x3 = cbind(Taxa = "Amoebozoa", x3)
names(x3) = c("Taxa", "Group", "Snow RA", "Spring Ice RA", "Summer Ice RA", "Cryoconite RA")

#combine and export as dataframe
df_1 = rbind(x2, x3)
write_csv(df_1, "../results/mt-rna-euk-phyla-order-habitat-profiles-phrag-amoe.csv")
