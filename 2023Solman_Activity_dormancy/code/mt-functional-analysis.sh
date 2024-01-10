#!/bin/bash

#Workflow for functional analysis of total RNA for paper
#'Active and dormant snow, ice and cryoconite communities of a High Arctic glacier'
#Solman A. B. T.; Jungblut A. D.; Larose, C.; Mourot, R.; Fox-Powell, M.; Benning, L. G.; C. Sanchez-Cid Torres; Bradley, J. A.
#Total RNA is extracted from snow, ice and cryoconite samples from a high Arctic Glacier.

#This is based on the ASaiM workflow (Batut et al. 2018)
#Step One: Quality control with FastQC, MultiQC and Trimmomatic
#Step Two: Filter rRNA and mRNA reads with SortMeRNA
#Step Three: Interlace forward and reverse reads with Seqfu
#Step Four: Use MetaPhlAn to assemble contigs and assign taxonomy (to contigs and single reads) based on rRNA and mRNA reads.
#Step Five: Visualise the communities using hclust2 and Graphlan.
#Step Six: Extract functional information from interlaced mRNA reads using HUMAnN.
#Step Seven: Normalise the gene families and metabolic pathways by sampling depth (i.e. relative abundance per sample).
#Step Eight: Identify which gene families are involved with which pathways using HUMAnN.
#Step Nine: Group gene families into Gene Ontology (GO) terms.
#Step Ten: Combine taxonomic and functional information.

#######################################################################################################

#Step One: Quality control with FastQC, MultiQC and Trimmomatic

#navigate to directory - this directory needs to have lots of space because we'll be downloading large databases
cd mbl_2022-mt/final/Fastq-files

#if fastq are in separate folders enter each folder and copy the contents to current directory.
#cp */* . 
#rm -r */ #Delete old directories. 

#make a folder for our fastqc output
mkdir FastQC

#generate fastqc files
fastqc Fastq-files/* -o FastQC

#generate multiQC report
multiqc FastQC/. -o FastQC

#open and view report
open multiqc_report.html

#check conda is up to date
conda update conda

#install trimmomatic if not already installed
conda install -c bioconda trimmomatic

#make a directory for our trimmed reads
mkdir Fastq-trimmed

#run trimmomatic on all files using a loop

#example trimmomatic function
# trimmomatic \ call trimmomatic
# PE \ PE stands for paired end
# -trimlog Fastq-trimmed/2023-01-18-trim-log \ #create a trim log file
# Fastq-3-sets-test/S21-25_S31_L001_R1_001.fastq.gz \ #forward read
# Fastq-3-sets-test/S21-25_S31_L001_R2_001.fastq.gz \ #reverse read
# Fastq-trimmed/S21-25_trim_forward_paired.fastq.gz \ #forward read paired output
# Fastq-trimmed/S21-25_trim_forward_unpaired.fastq.gz \ #forward read unpaired output
# Fastq-trimmed/S21-25_trim_reverse_paired.fastq.gz \ #reverse read paired output
# Fastq-trimmed/S21-25_trim_reverse_unpaired.fastq.gz \ #reverse read unpaired output
# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \ #remove adapters
# LEADING:20 \ #Remove leading low quality or N bases (below quality 20)
# TRAILING:20 \ #Remove trailing low quality or N bases (below quality 20)
# SLIDINGWINDOW:4:20 \ #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20
# MINLEN:50 #Drop reads below the 50 bases long

#also note you need the TruSeq3-PE.fa file in your current folder for this
#you can get this from the Trimmomatic github

for file in Fastq-files/*R1_001.fastq.gz; 
do 
R1="${file:12:15}_R1_001.fastq.gz";
R2="${file:12:15}_R2_001.fastq.gz";
R3="${file:12:6}";
R4="${file:12:6}";
echo "$R1"
echo "$R2"
echo "$R3"
echo "$R4"
trimmomatic \
PE \
-trimlog Fastq-trimmed/"$R1"-2023-09-15-trim-log \
Fastq-files/"$R1" \
Fastq-files/"$R2" \
Fastq-trimmed/"${R3}_forward_paired.fastq.gz" \
Fastq-trimmed/"${R3}_forward_unpaired.fastq.gz" \
Fastq-trimmed/"${R4}_reverse_paired.fastq.gz" \
Fastq-trimmed/"${R4}_reverse_unpaired.fastq.gz" \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
LEADING:20 \
TRAILING:20 \
SLIDINGWINDOW:4:20 \
MINLEN:50
done

#make a directory for storing QC info on trimmed files
mkdir FastQC-files-trimmed

#run fastqc on trimmed files
fastqc Fastq-trimmed/*_paired.fastq.gz -o FastQC-files-trimmed

#generate multiQC report
multiqc FastQC-files-trimmed/. -o FastQC-files-trimmed

#view the new quality report of the trimmed files
open FastQC-files-trimmed/multiqc_report.html

#looking at our new multiqc file we can see the quality has improved 
#and adapters have been removed

#######################################################################################################

#Step Two: Filter rRNA and mRNA reads with SortMeRNA

#create environment
conda create --name sortmerna_env

#activate environment
conda activate sortmerna_env

#install sortmerna
conda install sortmerna

#test the installation
sortmerna --version #4.3.6

#download the sortmerna github repository as this contains the fasta databases
git clone https://github.com/biocore/sortmerna.git

#make a folder for the sortmerna output
mkdir sortmerna-output

#if necessary look at the help output for sortmerna
#put the help info into a text file
sortmerna -h > sortmerna-help.txt
#view the text file
more sortmerna-help.txt #q to exit

#remember to run this command between each run as files are generated in this folder for each sample and need to be deleted
#rm /home/amys1/sortmerna/run/kvdb/*

#Example SortMeRNA function
# sortmerna \
# --reads Fastq-trimmed/S21-25_forward_paired.fastq.gz \ #forward reads to analyse
# --reads Fastq-trimmed/S21-25_reverse_paired.fastq.gz \ #reverse reads to analyse
# --ref sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta \ #database to look for rRNA reads
# --ref sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \ #database to look for rRNA reads
# --fastx \ #output aligned reads to fastq file
# --aligned Fastq-sorted/S21-25-rRNA \ #name of aligned reads file
# --other Fastq-sorted/S21-25-mRNA \ #name of unaligned reads file
# --out2 \ #output paired reads into separate files. Must be used with fastx.
# --paired_in \ #flags paired reads as aligned if just one of them is aligned. Must be used with fastx.
# --paired_out \ # both paired-end reads go in --other fasta/q file if just one of them isn't alligned. Mutally exclusive with paired_in
# --otu_map \ #output OTU map (input to QIIME's make_otu_table.py)
# rm /home/amys1/sortmerna/run/kvdb/* #purge this key value storage directory prior to each new run
# done
#run sortmerna 

for file in Fastq-trimmed/*forward_paired.fastq.gz; 
do 
R1="${file:14:6}";
echo "$R1"
sortmerna \
--reads Fastq-trimmed/"${R1}"_forward_paired.fastq.gz \
--reads Fastq-trimmed/"${R1}"_reverse_paired.fastq.gz \
--ref sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
--ref sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta \
--ref sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
--ref sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta \
--ref sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
--ref sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta \
--ref sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
--ref sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta \
--fastx \
--aligned Fastq-sorted/"${R1}"-rRNA \
--other Fastq-sorted/"${R1}"-mRNA \
--out2 \
--paired_out \
--otu_map
rm /home/amys1/sortmerna/run/kvdb/* 
done

#deactivate conda evironment
conda deactivate

# #Step Three: Interlace forward and reverse reads with Seqfu

# #create a new environment
# conda create --name seqfu_env

# #activate the environment
# conda activate seqfu_env

# #install seqfu
# conda install -c conda-forge -c bioconda "seqfu>1.10"

# #make a directory for interleaves/interlaced reads
# mkdir mRNA-interleaved

# #Example seqfu function
# # seqfu ilv \ #call intervieave function with ilv
# # -1 Fastq-sorted/S21-26-mRNA_fwd.fq.gz \ #give forward reads
# # -2 Fastq-sorted/S21-26-mRNA_rev.fq.gz > mRNA-interleaved/S21-26-mRNA_int.fq.gz #give reverse reads and name of file for interleaved reads

# #interleave all files in a loop
# for file in Fastq-sorted/*-mRNA_fwd.fq.gz; 
# do 
# R1="${file:14:5}";
# echo "$R1"
# seqfu ilv -1 Fastq-sorted/"${R1}"-mRNA_fwd.fq.gz -2 Fastq-sorted/"${R1}"-mRNA_rev.fq.gz > mRNA-interleaved/"${R1}"-mRNA_int.fq.gz
# done

# #interleave all files in a loop
# for file in Fastq-sorted/*-mRNA_fwd.fq.gz; 
# do 
# R1="${file:14:5}";
# echo "$R1"
# seqfu ilv -1 Fastq-sorted/"${R1}"-mRNA_fwd.fq.gz -2 Fastq-sorted/"${R1}"-mRNA_rev.fq.gz > mRNA-interleaved/"${R1}"-mRNA_int.fastq
# done

# #deactivate conda evironment
# conda deactivate

#######################################################################################################

#Step Four: Use MetaPhlAn to assemble contigs and assign taxonomy (to contigs and single reads) based on rRNA and mRNA reads.
#Metaphlan takes metagenomes/metatranscriptomes and identifies clades
#from phyla to species, and their realtive abundance
#it uses bowtie2 which is a fast, efficient short read aligner (a.k.a. if finds where our reads match a database or reference genome)

#create a new environment
conda create --name metaphlan_env

#activate the environment
conda activate metaphlan_env

#install metaphlan
conda install -c bioconda metaphlan

#test the installation
metaphlan --version #MetaPhlAn version 4.0.6 (1 Mar 2023)

#make an output file
mkdir mpa-out

#if necessary look at the help output for metaphlan
#put the help info into a text file
metaphlan -h > metaphlan-help.txt
#view the text file
more metaphlan-help.txt #q to exit

#note: when using metaphlan to analyse environmental samples it can be useful to lower stat_q and min_mapq_val.
# min_mapq_val = Minimum mapping quality value (MAPQ) [default 5] - this indicates the quality value below which an alignment will be discarded.
# Here I'm going to use 4.
# stat_q =  Quantile value for the robust average [default 0.2] - essentially, only markers that are most abundant,
# for example above the 20th percentile when using 0.2, will be used to calculate the relative abundance 
# of a taxa. Here I'm going to use 0.01 because I want to include nearly all marker genes.

#Example execution
# metaphlan \ #call metaphlan
# Fastq-trimmed/S21-26_forward_paired.fastq.gz,Fastq-trimmed/S21-26_reverse_paired.fastq.gz \ #give the locations of QC-ed forward and reverse reads (rRNA + mRNA)
# --bowtie2out mpa-out/S21-26.bowtie2.bz2 \ #save the intermediate bowtie output for re-running metaphlan quickly
# --nproc 5 \
# --stat_q 0.01 \
# --min_mapq_val 4\
# --unclassified_estimation \ #estimate unclassified fraction of the metagenomes/metatranscriptome
# --input_type fastq \ #what is the input file type - ours is fastq files!
# --add_viruses \ include viruses in the analysis
# -o mpa-out/S21-26.txt #this is the file in which the profiles metagenome/transcriptome will be stored

# If this is your first time running MetaPhlAn, the first step in this process involves MetaPhlAn 
# downloading, decompressing, and indexing its marker gene database.

#test on one file
metaphlan \
Fastq-trimmed/S21-26_forward_paired.fastq.gz,Fastq-trimmed/S21-26_reverse_paired.fastq.gz \
--bowtie2out mpa-out/S21-26.bowtie2.bz2 \
--nproc 5 \
--unclassified_estimation \
--stat_q 0.01 \
--min_mapq_val 4 \
--input_type fastq \
--add_viruses \
-o mpa-out/S21-26.txt

#metaphlan on all files
for file in Fastq-trimmed/*forward_paired.fastq.gz; 
do 
R1="${file:14:6}";
echo "$R1"
metaphlan \
Fastq-trimmed/"${R1}"_forward_paired.fastq.gz,Fastq-trimmed/"${R1}"_reverse_paired.fastq.gz \
--bowtie2out mpa-out/"${R1}".bowtie2.bz2 \
--nproc 5 \
--unclassified_estimation \
--add_viruses \
--stat_q 0.01 \
--min_mapq_val 4 \
--input_type fastq \
-o mpa-out/"${R1}".txt
done

#Merge output tables into one table
merge_metaphlan_tables.py mpa-out/*.txt > mpa-out/merged_abundance_table.txt

#deactivate metaphlan environment
conda deactivate 

#######################################################################################################

#Step Five: Visualise the MetaPhlan output using Graphlan.

#create a new environment
conda create --name graphlan_env

#activate the environment
conda activate graphlan_env

#install graphlan
conda install -c bioconda graphlan

#also install export2graphlan for creating annotation and tree files
conda install -c bioconda export2graphlan

#make an output directory for graphlan
mkdir graphlan-out

#export2graphlan takes the merged_abundance_table from metaphlan and 
#produces to files: 
# - merged_abundance_tree, a txt file listing all the identified clades in a order that denotes their connectedness alphabetically,
#for example, bacteria is top, then bacteria, acidobacteria, then all the classes, orders, families within that,
#then actinobacteria and so on.
# - merged_abundance.annot, a txt file that lists each unique class, it's marker colour and size. 
#it then goes on to list each clade at each phylogenetic level, along with the size and the colour of the marker if it is plotted.
# Edit this file to change labels and plotting colours

#to see more export2graphlan.py features print the help info to text file
#there are lots of ways to manipulate this data for plotting
export2graphlan.py --help > export2graphlan-help.txt
more export2graphlan-help.txt

#example export2graphlan command
# export2graphlan.py \
# --skip_rows 1,2 \  #ignore the title row and sample name row
# -i mpa-out/merged_abundance_table.txt \ #input file
# --tree graphlan-out/merged_abundance.tree.txt \ #tree output file name 
# --annotation graphlan-out/merged_abundance.annot.txt \ #annotation output file name
# --most_abundant 100 \ #show the 100 most abundant clades
# --abundance_threshold 0 \ #minimum abundance threshold for a clade to be plotted
# --least_biomarkers 0 \ #minimum number of biomarkers to extract.
# --annotations 5,6 \ #list which levels should be annotated in the tree. 5,6 = family + genus
# --external_annotations 5 \ #which levels should use external annotations - here 5 means the legend is for family level identifications
# --min_clade_size 1 #minimum number of clades that are biomarkers

#create graphlan input files
export2graphlan.py \
--skip_rows 1,2 \
-i mpa-out/merged_abundance_table.txt \
--tree graphlan-out/merged_abundance.tree.txt \
--annotation graphlan-out/merged_abundance.annot.txt \
--most_abundant 200 \
--abundance_threshold 0 \
--least_biomarkers 5 \
--annotations 5,6 \
--external_annotations 5,6 \
--max_font_size 10 \
--min_font_size 9 \
--min_clade_size 1

#At this point I manipulate the merged_abundance.annot file to change plotting labels and colours
#it is now call merged_abundance_edit.annot.txt

#create a cladogram
graphlan_annotate.py \
--annot graphlan-out/merged_abundance_edit.annot.txt graphlan-out/merged_abundance.tree.txt graphlan-out/merged_abundance.xml

graphlan.py \
--dpi 300 graphlan-out/merged_abundance.xml graphlan-out/merged_abundance.png \
--size 15 

#deactivate the environment
conda deactivate

#######################################################################################################

#Step Six: Extract functional information from interlaced mRNA reads using HUMAnN.

#HUMAnN (HMP Unified Metabolic Analysis Network) is a pipeline for profiling
#the presence/absence and abundance of microbial pathways from metagenomic/metatranscriptomic data.
#This is called functional profiling.
#When using metatranscriptome data we can answer the question, 
#what are microbes in these communities doing?

#Uniref database provides gene family definitions.
#MetaCyc database provides pathway definitions by gene family.
#Uses Bowtie2 and Diamond (translated alignment) to map reads to reference databases.

#create a new environment
conda create --name humann_env

#activate the environment
conda activate humann_env

#install humann
conda install -c biobakery humann

#which version?
humann --version #v3.6

#test installation
humann_test

#metaphlan is needed for running humann
metaphlan --version #v4.0.6

#make an output directory for humann
mkdir humann-out

#make a directory for databases
mkdir humann_dbs

#upgrade annotations database
humann_databases --download utility_mapping full humann_dbs --update-config yes

#upgrade pangenome database
#the chocoplan database is a database built by clusting coding sequences from
#genomes within the NCBI. These are systematically organized and annotated 
#microbial genomes and gene family clusters.
humann_databases --download chocophlan full humann_dbs --update-config yes

#upgrade protein database
#This is a translated search database.
#A translated search is where we look at which each ORF (open reading frames, a part of a sequence that may code for a protein)
# in a nucleotide sequence,
#and translate it into amino acid sequences which can then be searched against a protein database.
#here we have downloaded the full UniRef90 database
#Uniref = UniProt Reference Clusters
humann_databases --download uniref uniref90_diamond humann_dbs --update-config yes

#before using humann we need to concatenate forward and reverse reads into a single file

#copy our files and decompress them before concatenating
mkdir tRNA
cp Fastq-trimmed/* tRNA
gunzip tRNA/*

mkdir mRNA
cp Fastq-sorted/*mRNA* mRNA
gunzip mRNA/*

mkdir rRNA
cp Fastq-sorted/*rRNA* rRNA
gunzip rRNA/*

#remove any log files or unpaired reads from the above directories after copying and unzipping

#concatenate mRNA, rRNA and total RNA files
for file in Fastq-sorted/*mRNA_fwd.fq.gz; 
do 
R1="${file:13:6}";
echo "$R1"
cat mRNA/"$R1"-mRNA_fwd.fq mRNA/"$R1"-mRNA_rev.fq > mRNA/"$R1"-mRNA.fastq
cat rRNA/"$R1"-rRNA_fwd.fq rRNA/"$R1"-rRNA_rev.fq > rRNA/"$R1"-rRNA.fastq
cat tRNA/"$R1"_forward_paired.fastq tRNA/"$R1"_reverse_paired.fastq > tRNA/"$R1"-tRNA.fastq
done

#run HUMAnN on one total RNA concatenated file
humann \
--input tRNA/S21-31-tRNA.fastq \
--input-format fastq \
--output humann-out

#humann output files:
# OUTPUT_DIR/$SAMPLENAME_0.log
# OUTPUT_DIR/$SAMPLENAME_1_metaphlan_profile.tsv: Taxonomic output file from Metaphlan, how does this compare to the files produced when we ran metaphlan independently?
# I could try running HUMAnN with my previous Metaphlan files.
# OUTPUT_DIR/$SAMPLENAME_2_genefamilies.tsv: This file details the abundance of each gene family in the community. 
# Gene families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
# HUMAnN uses ChocoPhlAn (nucleotide alignment) and UniRef90 (translated search) databases to identify genes.
# Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length.
# RPK units reflect relative gene (or transcript) copy number in the community. 
# RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples.
# OUTPUT_DIR/$SAMPLENAME_3_reactions.tsv: Abundance of reactions in the community, after regrouping gene families to reactions.
# OUTPUT_DIR/$SAMPLENAME_4_pathabundance.tsv: Abundance of metabolic pathways in the community, as a function of the abundance
# of the pathways component reactions. 

#So essentially we're looking at
# gene families > reactions > metabolic pathways

#run HUMAnN on all total RNA reads
for file in tRNA/*-tRNA.fastq; 
do 
R1="${file:13:6}";
echo "$R1"
humann \
--input tRNA/"$R1"-tRNA.fastq \
--input-format fastq \
--output humann-out
done

#######################################################################################################

#Step Seven: Normalise the gene families and metabolic pathways by sampling depth (i.e. relative abundance per sample).

#for one file normalise the abundance output files to the size of each sample (relative abundance)
humann_renorm_table \
--input S21-49_2_genefamilies.tsv \
--output S21-49_2_genefamilies_relab.tsv \
--units relab

humann_renorm_table \
--input S21-49_4_pathabundance.tsv \
--output S21-49_4_pathabundance_relab.tsv \
--units relab

#normalise all gene and pathway abundance output files
for file in tRNA/*-tRNA.fastq; 
do 
R1="${file:13:6}";
echo "$R1"
humann_renorm_table \
--input "$R1"_2_genefamilies.tsv \
--output "$R1"_2_genefamilies_relab.tsv \
--units relab

humann_renorm_table \
--input "$R1"_4_pathabundance.tsv \
--output "$R1"_4_pathabundance_relab.tsv \
--units relab
done

#join the gene family and abundance files from all samples into one
humann_join_tables \
--input humman-out \
--output humann_2_genefamilies.tsv \
--file_name genefamilies_relab

humann_join_tables \
--input $OUTPUT_DIR \
--output humann_4_pathabundance.tsv \
--file_name pathabundance_relab

#barplot of humann output
humann_barplot \
--input humann_2_genefamilies.tsv \
#--focal-feature $FEATURE \
--outfile genefamilies_barplot

#######################################################################################################

#Step Eight: Identify which gene families are involved with which pathways using HUMAnN.
#Step Nine: Group gene families into Gene Ontology (GO) terms.
#Step Ten: Combine taxonomic and functional information.