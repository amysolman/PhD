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

#Step Three: Interlace forward and reverse reads with Seqfu

#create a new environment
conda create --name seqfu_env

#activate the environment
conda activate seqfu_env

#install seqfu
conda install -c conda-forge -c bioconda "seqfu>1.10"

#make a directory for interleaves/interlaced reads
mkdir mRNA-interleaved

#Example seqfu function
# seqfu ilv \ #call intervieave function with ilv
# -1 Fastq-sorted/S21-26-mRNA_fwd.fq.gz \ #give forward reads
# -2 Fastq-sorted/S21-26-mRNA_rev.fq.gz > mRNA-interleaved/S21-26-mRNA_int.fq.gz #give reverse reads and name of file for interleaved reads

#interleave all files in a loop
for file in Fastq-sorted/*-mRNA_fwd.fq.gz; 
do 
R1="${file:14:5}";
echo "$R1"
seqfu ilv -1 Fastq-sorted/"${R1}"-mRNA_fwd.fq.gz -2 Fastq-sorted/"${R1}"-mRNA_rev.fq.gz > mRNA-interleaved/"${R1}"-mRNA_int.fq.gz
done

#deactivate conda evironment
conda deactivate

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
mkdir metaphlan-output

#if necessary look at the help output for metaphlan
#put the help info into a text file
metaphlan -h > metaphlan-help.txt
#view the text file
more metaphlan-help.txt #q to exit

#Example execution
# metaphlan \ #call metaphlan
# Fastq-trimmed/S21-26_forward_paired.fastq.gz,Fastq-trimmed/S21-26_reverse_paired.fastq.gz \ #give the locations of QC-ed forward and reverse reads (rRNA + mRNA)
# --bowtie2db metaphlan_dbs \#location of the databases
# --bowtie2out metaphlan-output/S21-26.bowtie2.bz2 \ #save the intermediate bowtie output for re-running metaphlan quickly
# --nproc 5 \
# --unclassified_estimation \ #estimate unclassified fraction of the metagenomes/metatranscriptome
# --input_type fastq \ #what is the input file type - ours is fastq files!
# -o metaphlan-output/S21-26.txt #this is the file in which the profiles metagenome/transcriptome will be stored

# If this is your first time running MetaPhlAn, the first step in this process involves MetaPhlAn 
# downloading, decompressing, and indexing its marker gene database.

# If this download stalls, try deleting the databases and starting again:
rm -r /home/amys1/conda/envs/metaphlan_env/lib/python3.10/site-packages/metaphlan/metaphlan_databases

#download the databases 
metaphlan --install --bowtie2db metaphlan_dbs
#if the above doesn't work then try this...
wget --continue --limit-rate=20M http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar

#test on one file
metaphlan \
Fastq-trimmed/S21-26_forward_paired.fastq.gz,Fastq-trimmed/S21-26_reverse_paired.fastq.gz \
--bowtie2db metaphlan_dbs \
--bowtie2out metaphlan-output/S21-26.bowtie2.bz2 \
--nproc 5 \
--unclassified_estimation \
--input_type fastq \
-o metaphlan-output/S21-26.txt

#I think the database installed correctly!

#metaphlan on all files
for file in Fastq-trimmed/*forward_paired.fastq.gz; 
do 
R1="${file:14:6}";
echo "$R1"
metaphlan \
Fastq-trimmed/"${R1}"_forward_paired.fastq.gz,Fastq-trimmed/"${R1}"_reverse_paired.fastq.gz \
--bowtie2db metaphlan_dbs \
--bowtie2out metaphlan-output/"${R1}".bowtie2.bz2 \
--nproc 5 \
--unclassified_estimation \
--input_type fastq \
-o metaphlan-output/"${R1}".txt
done

#Merge output tables into one table
merge_metaphlan_tables.py metaphlan-output/*.txt > metaphlan-output/merged_abundance_table.txt

#deactivate metaphlan environment
conda deactivate 

#Step Five: Visualise the MetaPhlan output using heatmaps and Graphlan (publication).

#install and run hclust2

#create a new environment for hclust2 using mamba and install hclust
mamba create --name hclust2_env hclust2

#activate environment
mamba activate hclust2_env

#install hclust2
#conda install -c biobakery hclust2

#make a folder for hclust2 output
mkdir hclust2-out

#generate phylum level abundance table
grep -E "(s__)|(^ID)" metaphlan-output/merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > hclust2-out/merged_abundance_table_species.txt

#Create a heatmap of the communities
hclust2.py \
-i hclust2-out/merged_abundance_table_species.txt \
-o hclust2-out/HMP.sqrt_scale.png \
--skip_rows 1 \
--ftop 50 \
--f_dist_f correlation \
--s_dist_f braycurtis \
--cell_aspect_ratio 9 \
-s --fperc 99 \
--flabel_size 4 \
--metadata_rows 2,3,4 \
--legend_file hclust2-out/HMP.sqrt_scale.legend.png \
--max_flabel_len 100 \
--metadata_height 0.075 \
--minv 0.01 \
--no_slabels \
--dpi 300 \
--slinkage complete

#deactivate hclust2 environment
conda deactivate

#install and run graphlan

#create a new environment
conda create --name graphlan_env

#activate the environment
conda activate graphlan_env

#install graphlan
conda install -c biobakery graphlan

#make an output directory for graphlan
mkdir graphlan-out

#create graphlan input files
export2graphlan.py \
--skip_rows 1,2 \
-i metaphlan-output/merged_abundance_table.txt \
--tree graphlan-out/merged_abundance.tree.txt \
--annotation graphlan-out/merged_abundance.annot.txt \
--most_abundant 100 \
--abundance_threshold 1 \
--least_biomarkers 10 \
--annotations 5,6 \
--external_annotations 7 \
--min_clade_size 1

#create a cladogram
graphlan_annotate.py \
--annot graphlan-out/merged_abundance.annot.txt graphlan-out/merged_abundance.tree.txt graphlan-out/merged_abundance.xml
graphlan.py \
--dpi 300 graphlan-out/merged_abundance.xml graphlan-out/merged_abundance.png \
--external_legends

#Step Six: Extract functional information from interlaced mRNA reads using HUMAnN.

#create a new environment
conda create --name humann_env

#activate the environment
conda activate humann_env

#install humann
conda install -c biobakery humann

#test installation
humann --version

#metaphlan is needed for running humann
metaphlan --version

#make an output directory for humann
mkdir humann-out