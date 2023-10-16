#!/bin/bash

# Workflow for quality control, ASV identification and taxonomic classification

#navigate to folder
cd seq-data/ 

#make directory for qiime data
mkdir qiime-work-raw 

#make seperate directories for 16S and 18S
mkdir qiime-work-raw/16S 
mkdir qiime-work-raw/18S
mkdir qiime-work-raw/16S/raw-fastq-files 
mkdir qiime-work-raw/18S/raw-fastq-files

#copy over 16S and 18S files
cp FASTQ_raw/*_16* qiime-work-raw/16S/raw-fastq-files 
cp FASTQ_raw/*_18* qiime-work-raw/18S/raw-fastq-files

#make directory for fastqc reports
mkdir qiime-work-raw/16S/fastqc 
mkdir qiime-work-raw/18S/fastqc

#generate fastqc reports
fastqc qiime-work-raw/16S/raw-fastq-files/*fastq.gz -o qiime-work-raw/16S/fastqc 
fastqc qiime-work-raw/18S/raw-fastq-files/*fastq.gz -o qiime-work-raw/18S/fastqc

#generate multiqc reports
multiqc qiime-work-raw/16S/fastqc/ -o qiime-work-raw/16S 
multiqc qiime-work-raw/18S/fastqc/ -o qiime-work-raw/18S

#view multiqc reports
open qiime-work-raw/16S/multiqc_report.html 
open qiime-work-raw/18S/multiqc_report.html

#rename files to match Casava formating

#forward reads
for file in qiime-work-raw/1*S/raw-fastq-files/*R1_001.fastq.gz; 
do
R1="${file:0:44}_L001_R1_001.fastq.gz"
echo "$file"
echo "$R1"
mv -- "$file" "$R1"
rm "$file";
done

#reverse reads
for file in qiime-work-raw/1*S/raw-fastq-files/*R2_001.fastq.gz; 
do
R1="${file:0:44}_L001_R2_001.fastq.gz"
echo "$file"
echo "$R1"
mv -- "$file" "$R1"
rm "$file";
done

#activate qiime environment
conda activate qiime2-2023.2 

#convert sequences to qiime artefact
qiime tools import\
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path qiime-work-raw/16S/raw-fastq-files \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path qiime-work-raw/16S/demux-paired-end.qza

qiime tools import\
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path qiime-work-raw/18S/raw-fastq-files \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path qiime-work-raw/18S/demux-paired-end.qza

#visualise the raw sequences
qiime demux summarize
--i-data qiime-work-raw/16S/demux-paired-end.qza
--o-visualization qiime-work-raw/16S/demux-paired-end.qzv

qiime demux summarize
--i-data qiime-work-raw/18S/demux-paired-end.qza
--o-visualization qiime-work-raw/18S/demux-paired-end.qzv

#view the sequences
qiime tools view qiime-work-raw/16S/demux-paired-end.qzv

qiime tools view qiime-work-raw/18S/demux-paired-end.qzv

#use Cutadapt to remove primers
#16S
qiime cutadapt trim-paired \
--i-demultiplexed-sequences qiime-work-raw/16S/demux-paired-end.qza \
--p-cores 16 \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-adapter-f ATTAGAWACCCBNGTAGTCC \
--p-adapter-r TTACCGCGGCKGCTGRCAC \ 
--p-minimum-length 190 \
--p-discard-untrimmed \
--o-trimmed-sequences qiime-work-raw/16S/demux-paired-end-trimmed.qza \
--verbose

#18S
qiime cutadapt trim-paired \
--i-demultiplexed-sequences qiime-work-raw/18S/demux-paired-end.qza \
--p-cores 16 \
--p-front-f GTACACACCGCCCGTC \
--p-front-r TGATCCTTCTGCAGGTTCACCTAC \
--p-adapter-f GTAGGTGAACCTGCAGAAGGATCA \
--p-adapter-r GACGGGCGGTGTGTAC \
--p-discard-untrimmed \
--o-trimmed-sequences qiime-work-raw/18S/demux-paired-end-trimmed.qza \
--verbose

#run DADA2 to filter/trim, dereplicate, merge paired reads and make sequence table
#16S
time nohup qiime dada2 denoise-paired \
--i-demultiplexed-seqs qiime-work-raw/16S/demux-paired-end-trimmed.qza
--p-trim-left-f 35 \
--p-trim-left-r 35 \
--p-trunc-len-f 210 \
--p-trunc-len-r 190 \
--o-table qiime-work-raw/16S/DADA2_denoising_output/table.qza \
--o-representative-sequences qiime-work-raw/16S/DADA2_denoising_output/rep-seqs.qza \
--o-denoising-stats qiime-work-raw/16S/DADA2_denoising_output/denoising-stats.qza \
--verbose

#18S
time nohup qiime dada2 denoise-paired \
--i-demultiplexed-seqs qiime-work-raw/18S/demux-paired-end-trimmed.qza
--p-trunc-len-f 120 \
--p-trunc-len-r 110 \
--o-table qiime-work-raw/18S/DADA2_denoising_output/table.qza \
--o-representative-sequences qiime-work-raw/18S/DADA2_denoising_output/rep-seqs.qza \
--o-denoising-stats qiime-work-raw/18S/DADA2_denoising_output/denoising-stats.qza \
--verbose

#visualise denoising stats
#16S
qiime metadata tabulate \
--m-input-file qiime-work-raw/16S/DADA2_denoising_output/denoising_stats.qza \
--o-visualization qiime-work-raw/16S/DADA2_denoising_output/denoising_stats.qzv

#18S
qiime metadata tabulate \
--m-input-file qiime-work-raw/18S/DADA2_denoising_output/denoising_stats.qza \
--o-visualization qiime-work-raw/18S/DADA2_denoising_output/denoising_stats.qzv

#view denoising stats
#16S
qiime tools view qiime-work-raw/16S/DADA2_denoising_output/denoising_stats.qzv

#18S
qiime tools view qiime-work-raw/18S/DADA2_denoising_output/denoising_stats.qzv

#visualise representative sequences
#16S
qiime feature-table tabulate-seqs \
--i-data qiime-work-raw/16S/DADA2_denoising_output/rep_seq.qza \
--o-visualization qiime-work-raw/16S/DADA2_denoising_output/rep-seqs.qzv

#18S
qiime feature-table tabulate-seqs \
--i-data qiime-work-raw/18S/DADA2_denoising_output/rep_seq.qza \
--o-visualization qiime-work-raw/18S/DADA2_denoising_output/rep-seqs.qzv

#view representative sequences
#16S
qiime tools view qiime-work-raw/16S/DADA2_denoising_output/rep-seqs.qzv

#18S
qiime tools view qiime-work-raw/18S/DADA2_denoising_output/rep-seqs.qzv

#visualise ASV table
#16S
time qiime feature-table summarize \
--i-table qiime-work-raw/16S/DADA2_denoising_output/table.qza \
--m-sample-metadata-file qiime-work-raw/16S/metadata.tsv \
--o-visualization qiime-work-raw/16S/DADA2_denoising_output/table.qzv

#18S
time qiime feature-table summarize \
--i-table qiime-work-raw/18S/DADA2_denoising_output/table.qza \
--m-sample-metadata-file qiime-work-raw/18S/metadata.tsv \
--o-visualization qiime-work-raw/18S/DADA2_denoising_output/table.qzv

#view ASV table
#16S
qiime tools view qiime-work-raw/16S/DADA2_denoising_output/table.qzv

#18S
qiime tools view qiime-work-raw/18S/DADA2_denoising_output/table.qzv

#extract reference reads from the SILVA 138 database
#16S
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--o-reads qiime-work-raw/16S/ref-seqs

#18S
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads qiime-work-raw/18S/ref-seqs

#use the extracted ref reads/tax files to create trained classifier 
#16S
nohup qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads qiime-work-raw/16S/ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier qiime-work-raw/16S/classifier.qza

#18S
nohup qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads qiime-work-raw/18S/ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier qiime-work-raw/18S/classifier.qza

#use the classifier to taxonomically classify the ASVs
#16S
qiime feature-classifier classify-sklearn \
--i-classifier qiime-work-raw/16S/classifier.qza \
--i-reads qiime-work-raw/16S/rep-seqs.qza \
--o-classification qiime-work-raw/16S/taxonomy.qza

#18S
qiime feature-classifier classify-sklearn \
--i-classifier qiime-work-raw/18S/classifier.qza \
--i-reads qiime-work-raw/18S/rep-seqs.qza \
--o-classification qiime-work-raw/18S/taxonomy.qza

#visualise the taxonomy
#16S
qiime metadata tabulate \
--m-input-file qiime-work-raw/16S/taxonomy.qza \
--o-visualization qiime-work-raw/16S/taxonomy.qzv

#18S
qiime metadata tabulate \
--m-input-file qiime-work-raw/18S/taxonomy.qza \
--o-visualization qiime-work-raw/18S/taxonomy.qzv

#make folders to keep our exported files
mkdir qiime-work-raw/16S/export
mkdir qiime-work-raw/18S/export

#MAFFT alignment
#16S
time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences qiime-work-raw/16S/DADA2_denoising_output/rep_seq.qza \
--o-alignment qiime-work-raw/16S/export/aligned-rep-seqs.qza \
--o-masked-alignment qiime-work-raw/16S/export/masked-aligned-rep-seqs.qza \
--o-tree qiime-work-raw/16S/export/unrooted-tree.qza \
--o-rooted-tree qiime-work-raw/16S/export/rooted-tree.qza

#18S
time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences qiime-work-raw/18S/DADA2_denoising_output/rep_seq.qza \
--o-alignment qiime-work-raw/18S/export/aligned-rep-seqs.qza \
--o-masked-alignment qiime-work-raw/18S/export/masked-aligned-rep-seqs.qza \
--o-tree qiime-work-raw/18S/export/unrooted-tree.qza \
--o-rooted-tree qiime-work-raw/18S/export/rooted-tree.qza

#export the ASV table into biom format
#16S
qiime tools export \
--input-path qiime-work-raw/16S/DADA2_denoising_output/table.qza \
--output-path qiime-work-raw/16S/export

#18S
qiime tools export \
--input-path qiime-work-raw/18S/DADA2_denoising_output/table.qza \
--output-path qiime-work-raw/18S/export

#convert the biom format into tsv
#16S
biom convert \
-i qiime-work-raw/16S/export/feature-table.biom \
-o qiime-work-raw/16S/export/table.tsv --to-tsv

#18S
biom convert \
-i qiime-work-raw/18S/export/feature-table.biom \
-o qiime-work-raw/18S/export/table.tsv --to-tsv

#export representative sequences
#16S
qiime tools export \
--input-path qiime-work-raw/16S/DADA2_denoising_output/rep_seq.qza \
--output-path qiime-work-raw/16S/export

#18S
qiime tools export \
--input-path qiime-work-raw/18S/DADA2_denoising_output/rep_seq.qza \
--output-path qiime-work-raw/18S/export

#export taxonomy file
#16S
qiime tools export \
--input-path qiime-work-raw/16S/taxonomy.qza \
--output-path qiime-work-raw/16S/export

#18S
qiime tools export \
--input-path qiime-work-raw/18S/taxonomy.qza \
--output-path qiime-work-raw/18S/export

#export the phylogenetic tree
#16S
qiime tools export \
--input-path qiime-work-raw/16S/export/rooted-tree.qza \
--output-path qiime-work-raw/16S/export

#18S
qiime tools export \
--input-path qiime-work-raw/18S/export/rooted-tree.qza \
--output-path qiime-work-raw/18S/export

#deactivate conda environment
conda deactivate