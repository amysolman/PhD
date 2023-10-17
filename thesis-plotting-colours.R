#assigning plotting colours for thesis

#Clear workspace and load packages
rm(list=ls())
graphics.off()

library(phyloseq)

#Import data

#chapter one amplicons
pro.c1 <- readRDS("01-chapter-one/results/16S-phylo-object.rds") 
euk.c1 <- readRDS("01-chapter-one/results/18S-phylo-object.rds")

#chapter two amplicons
#prokaryotes
pro.c2 <- readRDS("02-chapter-two/results-proportional/16S-ps-norm.rds") 
#eukaryotes
euk.c2 <- readRDS("02-chapter-two/results-proportional/18S-no-mm-ps-norm.rds") 
#microfauna
mm <- readRDS("02-chapter-two/results-proportional/18S-mm-only-ps-norm.rds") 

#metatranscriptomics
#prokaryotes
pro.DNA.c3 <- readRDS("03-chapter-three/metatranscriptomics-analysis/results/16S-phylo-object.rds") 
pro.RNA.c3 <- readRDS("03-chapter-three/metatranscriptomics-analysis/results/mt-pro-phylo-object.rds") 
#eukaryotes
euk.DNA.c3 <- readRDS("03-chapter-three/metatranscriptomics-analysis/results/18S-phylo-object-micro-keep.rds") 
euk.RNA.c3 <- readRDS("03-chapter-three/metatranscriptomics-analysis/results/mt-euk-phylo-object.rds") 

#PMA
#prokaryotes
pro.PMA.c3 <- readRDS("03-chapter-three/pma/results/16S-phylo-object-norm-rel.rds") 
#eukaryotes without micrometazoans
euk.PMA.c3 <- readRDS("03-chapter-three/pma/results/18S-phylo-object-norm-rel.rds") 

#extract all taxonomy tables, combine and keep only unique phyla

#prokaryotes

pro.col.list = c("antiquewhite", "antiquewhite4",
                 "aquamarine", "aquamarine4",
             "orange", "darkorange4",
             "blue", "darkblue",
             "lightcyan1", "lightcyan4",
             "darkgoldenrod1", "darkgoldenrod4",
             "darkolivegreen1", "darkolivegreen",
             "tan", "tan4", 
             "cyan", "cyan4",
             "lavenderblush1", "lavenderblush4",
             "chocolate", "chocolate4",
             "snow1", "snow4",
             "peachpuff1", "peachpuff4",
             "turquoise", "turquoise4",
             "red", "darkred",
             "mistyrose1", "mistyrose4",
             "deepskyblue", "deepskyblue4",
             "brown1", "brown4",
             "green", "darkgreen",
             "coral", "coral4", 
             "white", "aliceblue",
             "antiquewhite1", "antiquewhite4",
             "azure", "azure4",
             "bisque", "bisque4",
             "burlywood1", "burlywood4",
             "cadetblue1", "cadetblue4",
             "chartreuse", "chartreuse4",
             "cornsilk", "cornsilk4",
             "darkseagreen1", "darkseagreen4",
             "darkslategray1", "darkslategray4",
             "firebrick1", "firebrick4",
             "gold1", "gold4",
             "lightsteelblue1", "lightsteelblue4",
             "violetred1", "violetred4", 
             "lemonchiffon1", "lemonchiffon4",
             "khaki1", "khaki4",
             "ivory1", "ivory4",
             "indianred1", "indianred4",
             "hotpink1", "hotpink4",
             "honeydew1", "honeydew4",
             "maroon1", "maroon4",
             "magenta1", "magenta4",
             "lightyellow1", "lightyellow4",
             "gray100", "gray40",
             "lightskyblue1", "lightskyblue4",
             "lightsalmon1", "lightsalmon4",
             "lightpink1", "lightpink4",
             "yellow", "yellow4", 
             "palevioletred1", "palevioletred4",
             "paleturquoise", "paleturquoise4",
             "palegreen1", "palegreen4",
             "orchid1", "orchid4",
             "orangered1", "orangered4",
             "olivedrab1", "olivedrab4",
             "navajowhite1", "navajowhite4",
             "pink", "deeppink4", 
             "mediumpurple1", "mediumpurple4", #57
             "purple","darkorchid4", 
             "plum1", "plum4",
             "rosybrown1", "rosybrown4",
             "royalblue1", "royalblue4",
             "salmon1", "salmon4",
             "gray80", "gray100",
             "seagreen1", "seagreen4",
             "seashell1", "seashell4",
             "sienna1", "sienna4",
             "skyblue1", "skyblue4",
             "slateblue1", "slateblue4",
             "slategray1", "slategray4", #68
             "darkgoldenrod1", "darkgoldenrod4",
             "springgreen1", "springgreen4",
             "steelblue1", "steelblue4", 
             "thistle1", "thistle4",
             "tomato1", "tomato4",
             "wheat1", "wheat4",
             "gray30", "gray60")

pro.tax = rbind.fill(data.frame(tax_table(pro.c1)), data.frame(tax_table(pro.c2)), data.frame(tax_table(pro.DNA.c3)), 
      data.frame(tax_table(pro.RNA.c3)), data.frame(tax_table(pro.PMA.c3)))
#remove eukaryotes and keep only phyla data
pro.tax2 = pro.tax[pro.tax$Domain != "Eukaryota" | pro.tax$Kingdom != "Eukaryota",]
pro.phyla = data.frame(Phylum = sort(unique(na.omit(pro.tax2$Phylum))), 
                       Colour1=pro.col.list[c(TRUE, FALSE)], 
                       Colour2=pro.col.list[c(FALSE, TRUE)])

#save plotting colours
write.csv(pro.phyla, "thesis-plotting-colours-prokaryotes.csv")

#eukaryotes

euk.col.list = c("aquamarine", "aquamarine4",
                 "antiquewhite", "antiquewhite4",
                 "lavenderblush1", "lavenderblush4",
                 "snow1", "snow4",
                 "lightcyan1", "lightcyan4",
                 "darkgoldenrod1", "darkgoldenrod4",
                 "antiquewhite1", "antiquewhite4",
                 "plum1", "plum4",
                 "cyan", "cyan4",
                 "orange", "darkorange4",
                 "chocolate", "chocolate4",
                 "blue", "darkblue",
                 "peachpuff1", "peachpuff4",
                 "cornsilk", "cornsilk4",
                 "white", "aliceblue",
                 "mistyrose1", "mistyrose4",
                 "springgreen1", "springgreen4",
                 "violetred1", "violetred4", 
                 "coral", "coral4", 
                 "green", "darkgreen",
                 "red", "darkred",
                 "darkolivegreen1", "darkolivegreen",
                 "azure", "azure4",
                 "bisque", "bisque4",
                 "burlywood1", "burlywood4",
                 "purple","darkorchid4", 
                 "chartreuse", "chartreuse4",
                 "turquoise", "turquoise4",
                 "darkseagreen1", "darkseagreen4",
                 "darkslategray1", "darkslategray4",
                 "firebrick1", "firebrick4",
                 "gold1", "gold4",
                 "lightsteelblue1", "lightsteelblue4",
                 "brown1", "brown4",
                 "lemonchiffon1", "lemonchiffon4",
                 "khaki1", "khaki4",
                 "ivory1", "ivory4",
                 "indianred1", "indianred4",
                 "hotpink1", "hotpink4",
                 "honeydew1", "honeydew4",
                 "maroon1", "maroon4",
                 "magenta1", "magenta4",
                 "lightyellow1", "lightyellow4",
                 "gray100", "gray40",
                 "lightskyblue1", "lightskyblue4",
                 "lightsalmon1", "lightsalmon4",
                 "lightpink1", "lightpink4",
                 "lavenderblush1", "lavenderblush4",
                 "palevioletred1", "palevioletred4",
                 "paleturquoise", "paleturquoise4",
                 "palegreen1", "palegreen4",
                 "orchid1", "orchid4",
                 "orangered1", "orangered4",
                 "olivedrab1", "olivedrab4",
                 "navajowhite1", "navajowhite4",
                 "seashell1", "seashell4",
                 "mediumpurple1", "mediumpurple4", #57
                 "cadetblue1", "cadetblue4",
                 "tan", "tan4", 
                 "rosybrown1", "rosybrown4",
                 "royalblue1", "royalblue4",
                 "salmon1", "salmon4",
                 "gray80", "gray100",
                 "seagreen1", "seagreen4",
                 "pink", "deeppink4", 
                 "sienna1", "sienna4",
                 "skyblue1", "skyblue4",
                 "slateblue1", "slateblue4",
                 "slategray1", "slategray4", #68
                 "darkgoldenrod1", "darkgoldenrod4",
                 "deepskyblue", "deepskyblue4",
                 "steelblue1", "steelblue4", 
                 "thistle1", "thistle4",
                 "tomato1", "tomato4",
                 "wheat1", "wheat4",
                 "gray30", "gray60",
                 "antiquewhite", "antiquewhite4",
                 "aquamarine", "aquamarine4",
                 "orange", "darkorange4",
                 "blue", "darkblue",
                 "lightcyan1", "lightcyan4",
                 "darkgoldenrod1", "darkgoldenrod4",
                 "darkolivegreen1", "darkolivegreen",
                 "tan", "tan4", 
                 "cyan", "cyan4",
                 "yellow", "yellow4", 
                 "chocolate", "chocolate4",
                 "snow1", "snow4",
                 "turquoise", "turquoise4")

euk.tax = rbind.fill(data.frame(tax_table(euk.c1)), data.frame(tax_table(euk.c2)), data.frame(tax_table(euk.DNA.c3)), 
                     data.frame(tax_table(euk.RNA.c3)), data.frame(tax_table(euk.PMA.c3)))
#remove eukaryotes and keep only phyla data
euk.tax2 = euk.tax[euk.tax$Domain != "Bacteria" | euk.tax$Kingdom != "Bacteria",]
euk.phyla = data.frame(Phylum = sort(unique(na.omit(euk.tax2$Phylum))), 
                       Colour1=euk.col.list[c(TRUE, FALSE)], 
                       Colour2=euk.col.list[c(FALSE, TRUE)])

#save plotting colours
write.csv(euk.phyla, "thesis-plotting-colours-eukaryotes.csv")



