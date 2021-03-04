### RScript to go along with the tutorial found here: https://astrobiomike.github.io/amplicon/dada2_workflow_ex


library(dada2)
packageVersion("dada2") # 1.11.5 when this was put together

setwd("~/dada2_amplicon_ex_workflow")

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
  # one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples", what="character")

  # one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_sub_R1_trimmed.fq.gz")
  # and one with the reverse
reverse_reads <- paste0(samples, "_sub_R2_trimmed.fq.gz")

  # and variables holding file names for the forward and reverse
  # filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")



plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
  # and just plotting the last 4 samples of the reverse reads
plotQualityProfile(reverse_reads[17:20])




filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                rm.phix=TRUE, minLen=175, truncLen=c(250,200))


class(filtered_out) # matrix
dim(filtered_out) # 20 2

filtered_out


plotQualityProfile(filtered_forward_reads)
plotQualityProfile(filtered_reverse_reads)
plotQualityProfile(filtered_reverse_reads[17:20])



err_forward_reads <- learnErrors(filtered_forward_reads)
# err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE) # problem running this way if on Binder
err_reverse_reads <- learnErrors(filtered_reverse_reads)
# err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE) # problem running this way if on Binder


plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)




derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples



dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo")
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder


merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                    derep_reverse, trimOverhang=TRUE, minOverlap=170)

 # this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"



seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521


seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 17 bimeras out of 2521 input sequences.

  # though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9931372 # good, we barely lost any in terms of abundance




  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab




taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", tryRC=T)
# taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=T) # problem running this way if on Binder



  # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

  # tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)




library(decontam)
packageVersion("decontam") # 1.1.2 when this was put together



colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants



  # getting IDs of those identified as likely contaminants
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

  # making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

  # making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

  # making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

  ## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)


### Analysis in R ###

  # if you're following/working with the tutorial data, make sure you're in the correct place still
setwd("~/amplicon_example_workflow/")
list.files()

## here are installing packages, though if you're working in the Binder environment, these are already taken care of
# install.packages("phyloseq")
# install.packages("vegan")
# install.packages("DESeq2")
# install.packages("ggplot2")
# install.packages("dendextend")
# install.packages("tidyr")
# install.packages("viridis")
# install.packages("reshape")



# source("https://bioconductor.org/biocLite.R")
# BiocManager::install("phyloseq")
# BiocManager::install("DESeq2")

	# and here is an example using devtools:
# install.packages("devtools")
# library("devtools")
# devtools::install_github("tidyverse/ggplot2")


library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")


packageVersion("phyloseq") # 1.26.1
packageVersion("vegan") # 2.5.4
packageVersion("DESeq2") # 1.22.2
packageVersion("ggplot2") # 3.1.1
packageVersion("dendextend") # 1.10.0
packageVersion("tidyr") # 0.8.3
packageVersion("viridis") # 0.5.1
packageVersion("reshape") # 0.8.8




## NOTE ##
# if you loaded the saved R data above with `load(".RData")` and didn't run all the
# steps, you may want to write out these files here before clearing the environment
# it won't hurt if you did it already, so better to just run these here :)

write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

# ok, now moving on

rm(list=ls())
  
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
             check.names=F, sep="\t")[ , -c(1:4)]

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
           row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                   check.names=F, sep="\t")
  
  # and setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

sample_info_tab # to take a peek





  # first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~type) 
    # we have to include the "colData" and "design" arguments because they are 
    # required, as they are needed for further downstream processing by DESeq2, 
    # but for our purposes of simply transforming the data right now, they don't 
    # matter
  
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
    # NOTE: If you get this error here with your dataset: "Error in
    # estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
    # gene contains at least one zero, cannot compute log geometric means", that
    # can be because the count table is sparse with many zeroes, which is common
    # with marker-gene surveys. In that case you'd need to use a specific
    # function first that is equipped to deal with that. You could run:
      # deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
    # now followed by the transformation function:
      # deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

  # and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

  # and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))




euc_clust <- hclust(euc_dist, method="ward.D2")

  # hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
    # but i like to change them to dendrograms for two reasons:
      # 1) it's easier to color the dendrogram plot by groups
      # 2) if wanted you can rotate clusters with the rotate() 
      #    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")




  # making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

  # generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
  theme(legend.position="none")



rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)

  # and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))




  # first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

  # and now we can call the plot_richness() function on our phyloseq object
plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
  theme(legend.title = element_blank())



plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
  theme(legend.title = element_blank())



  # using phyloseq to make a count table that has summed all ASVs
    # that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum")) 

  # making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Phylum"))[,2]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

  # we also have to account for sequences that weren't assigned any
    # taxonomy even at the phylum level 
    # these came into R as 'NAs' in the taxonomy table, but their counts are
    # still in the count table
    # so we can get that value for each sample by substracting the column sums
    # of this new table (that has everything that had a phylum assigned to it)
    # from the column sums of the starting count table (that has all
    # representative sequences)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
  # and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

  # now we'll remove the Proteobacteria, so we can next add them back in
    # broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

  # making count table broken down by class (contains classes beyond the
    # Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Class")) 

  # making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="Class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)

  # making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"])

  # changing the row names like above so that they correspond to the taxonomy,
    # rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$Class) 

  # making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

  # there are also possibly some some sequences that were resolved to the level
    # of Proteobacteria, but not any further, and therefore would be missing from
    # our class table
    # we can find the sum of them by subtracting the proteo class count table
    # from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

  # now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

  # and to check we didn't miss any other sequences, we can compare the column
    # sums to see if they are the same
    # if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

  # now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

  # if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
  # we see there are currently 44 rows, which might be a little busy for a
    # summary figure
    # many of these taxa make up a very small percentage, so we're going to
    # filter some out
    # this is a completely arbitrary decision solely to ease visualization and
    # intepretation, entirely up to your data and you
    # here, we'll only keep rows (taxa) that make up greater than 5% in any
    # individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
  # checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
    # now we have 13, much more manageable for an overview figure

  # though each of the filtered taxa made up less than 5% alone, together they
    # may add up and should still be included in the overall summary
    # so we're going to add a row called "Other" that keeps track of how much we
    # filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)




  # first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

  # and add a column of the taxa names so that it is within the table, rather
  # than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

  # now we'll transform the table into narrow, or long, format (also makes
  # plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)

  # take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
    # manipulating tables like this is something you may need to do frequently in R

  # now we want a table with "color" and "characteristics" of each sample to
    # merge into our plotting table so we can use that more easily in our plotting
    # function
    # here we're making a new table by pulling what we want from the sample
    # information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab), "char"=sample_info_tab$char, "color"=sample_info_tab$color, stringsAsFactors=F)

  # and here we are merging this table with the plotting table we just made
    # (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

  # and now we're ready to make some summary figures with our wonderfully
    # constructed table

  ## a good color scheme can be hard to find, i included the viridis package
    ## here because it's color-blind friendly and sometimes it's been really
    ## helpful for me, though this is not demonstrated in all of the following :/ 

  # one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")



ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
  geom_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")




  # let's set some helpful variables first:
bw_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "water"]
rock_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "rock"]

  # first we need to subset our plotting table to include just the rock samples to plot
filt_major_taxa_proportions_rocks_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% rock_sample_IDs, ]
  # and then just the water samples
filt_major_taxa_proportions_water_samples_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% bw_sample_IDs, ]

  # and now we can use the same code as above just with whatever minor alterations we want
  # rock samples
ggplot(filt_major_taxa_proportions_rocks_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Rock samples only")

  # water samples
ggplot(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
  scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
  geom_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +
  scale_color_manual(values=unique(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$char)])) +
  geom_boxplot(fill=NA, outlier.color=NA) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
  labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Bottom water samples only")




  # notice we're leaving off the "char" and "color" columns, in the code, and be sure to peak at the tables after making them
rock_sample_major_taxa_proportion_tab <- cast(filt_major_taxa_proportions_rocks_only_tab_for_plot.g[, c(1:3)], Major_Taxa ~ Sample)
water_sample_major_taxa_proportion_tab <- cast(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g[, c(1:3)], Major_Taxa ~ Sample)
  # summing each taxa across all samples for both groups 
rock_sample_summed_major_taxa_proportions_vec <- rowSums(rock_sample_major_taxa_proportion_tab)
water_sample_summed_major_taxa_proportions_vec <- rowSums(water_sample_major_taxa_proportion_tab)

rock_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(rock_sample_summed_major_taxa_proportions_vec), "Proportion"=rock_sample_summed_major_taxa_proportions_vec, row.names=NULL)
water_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(water_sample_summed_major_taxa_proportions_vec), "Proportion"=water_sample_summed_major_taxa_proportions_vec, row.names=NULL)

  # plotting just rocks
ggplot(data.frame(rock_sample_major_taxa_summary_tab), aes(x="Rock samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Rock samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

  # and plotting just water samples
ggplot(data.frame(water_sample_major_taxa_summary_tab), aes(x="Bottom water samples", y=Proportion, fill=Major_Taxa)) + 
  geom_bar(width=1, stat="identity") +
  coord_polar("y") +
  scale_fill_viridis(discrete=TRUE) +
  ggtitle("Water samples only") +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())



anova(betadisper(euc_dist, sample_info_tab$type)) # 0.002



  # first we'll need to go back to our transformed table, and generate a
    # distance matrix only incorporating the basalt samples
      # and to help with that I'm making a variable that holds all basalt rock
      # names (just removing the single calcium carbonate sample, R7)
basalt_sample_IDs <- rock_sample_IDs[!rock_sample_IDs %in% "R7"]

  # new distance matrix of only basalts
basalt_euc_dist <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% basalt_sample_IDs]))

  # and now making a sample info table with just the basalts
basalt_sample_info_tab <- sample_info_tab[row.names(sample_info_tab) %in% basalt_sample_IDs, ]

  # running betadisper on just these based on level of alteration as shown in the images above:
anova(betadisper(basalt_euc_dist, basalt_sample_info_tab$char)) # 0.7




adonis(basalt_euc_dist~basalt_sample_info_tab$char) # 0.003



  # making our phyloseq object with transformed table
basalt_vst_count_phy <- otu_table(vst_trans_count_tab[, colnames(vst_trans_count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_sample_info_tab_phy <- sample_data(basalt_sample_info_tab)
basalt_vst_physeq <- phyloseq(basalt_vst_count_phy, basalt_sample_info_tab_phy)

  # generating and visualizing the PCoA with phyloseq
basalt_vst_pcoa <- ordinate(basalt_vst_physeq, method="MDS", distance="euclidean")
basalt_eigen_vals <- basalt_vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

  # and making our new ordination of just basalts with our adonis statistic
plot_ordination(basalt_vst_physeq, basalt_vst_pcoa, color="char") + 
  labs(col="type") + geom_point(size=1) + 
  geom_text(aes(label=rownames(basalt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  annotate("text", x=25, y=73, label="Highly altered vs glassy") +
  annotate("text", x=25, y=67, label="Permutational ANOVA = 0.005") + 
  coord_fixed(sqrt(basalt_eigen_vals[2]/basalt_eigen_vals[1])) + ggtitle("PCoA - basalts only") + 
  scale_color_manual(values=unique(basalt_sample_info_tab$color[order(basalt_sample_info_tab$char)])) + 
  theme(legend.position="none")




  # first making a basalt-only phyloseq object of non-transformed values (as that is what DESeq2 operates on
basalt_count_phy <- otu_table(count_tab[, colnames(count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_count_physeq <- phyloseq(basalt_count_phy, basalt_sample_info_tab_phy)
  
  # now converting our phyloseq object to a deseq object
basalt_deseq <- phyloseq_to_deseq2(basalt_count_physeq, ~char)

  # and running deseq standard analysis:
basalt_deseq <- DESeq(basalt_deseq)





  # pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_altered_vs_glassy <- results(basalt_deseq, alpha=0.01, contrast=c("char", "altered", "glassy"))

  # we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_altered_vs_glassy) 
    # this tells us out of ~1,800 ASVs, with adj-p < 0.01, there are 7 increased when comparing altered basalts to glassy basalts, and about 6 decreased
    # "decreased" in this case means at a lower count abundance in the altered basalts than in the glassy basalts, and "increased" means greater proportion in altered than in glassy
    # remember, this is done with a drastically reduced dataset, which is hindering the capabilities here quite a bit i'm sure

  # let's subset this table to only include these that pass our specified significance level
sigtab_res_deseq_altered_vs_glassy <- deseq_res_altered_vs_glassy[which(deseq_res_altered_vs_glassy$padj < 0.01), ]

  # now we can see this table only contains those we consider significantly differentially abundant
summary(sigtab_res_deseq_altered_vs_glassy) 

  # next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_deseq_altered_vs_glassy_with_tax <- cbind(as(sigtab_res_deseq_altered_vs_glassy, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_deseq_altered_vs_glassy), ], "matrix"))

  # and now let's sort that table by the baseMean column
sigtab_deseq_altered_vs_glassy_with_tax[order(sigtab_deseq_altered_vs_glassy_with_tax$baseMean, decreasing=T), ]
