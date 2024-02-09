#!/bin/bash

# Workflow for quality control, ASV identification and taxonomic classification

#make working folders
mkdir 16S
mkdir 18S

#Process 16S rRNA V4 amplified reads 

#move to your working folder
cd 16S

#activate qiime
conda activate qiime2-2023.2

#run DADA2 to filter/trim, dereplicate, merge paired reads and make sequence table
time nohup qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux-paired-end.qza \
--p-trim-left-f 20 \
--p-trim-left-r 20 \
--p-trunc-len-f 240 \
--p-trunc-len-r 134 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza \
--verbose

#visualise ASV table (include simple metadata file)
time qiime feature-table summarize \
--i-table table.qza \
--m-sample-metadata-file metadata.tsv \
--o-visualization table.qzv

#view ASV table
qiime tools view table.qzv

#extract reference reads from the SILVA 138 database using the older EMP 16S primers
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--o-reads ref-seqs

#use the extracted ref reads/tax files to create trained classifier 
nohup qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier classifier.qza

#use the classifier to taxonomically classify the ASVs
qiime feature-classifier classify-sklearn \
--i-classifier classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

#visualise the taxonomy
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

#make a folder to keep our exported files
mkdir export

#export the ASV table into biom format
qiime tools export \
--input-path table.qza \
--output-path export

#convert the biom format into tsv
biom convert \
-i export/feature-table.biom \
-o export/table.tsv --to-tsv

#export representative sequences
qiime tools export \
--input-path rep-seqs.qza \
--output-path export

#export taxonomy file
qiime tools export \
--input-path taxonomy.qza \
--output-path export

#export the phylogenetic tree
qiime tools export \
--input-path rooted-tree.qza \
--output-path export


#Process 18S rRNA V9 amplified reads 

#move to your working folder
cd 18S

#run DADA2 to filter/trim, dereplicate, merge paired reads and make sequence table
time nohup qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux-paired-end.qza \
--p-trim-left-f 20 \
--p-trim-left-r 20 \
--p-trunc-len-f 130 \
--p-trunc-len-r 76 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza \
--verbose

#visualise ASV table (include simple metadata file)
time qiime feature-table summarize \
--i-table table.qza \
--m-sample-metadata-file metadata.tsv \
--o-visualization table.qzv

#view ASV table
qiime tools view table.qzv

#extract reference reads from the SILVA 138 database using the older EMP 16S primers
qiime feature-classifier extract-reads \
--i-sequences silva-138-99-seqs.qza \
--p-f-primer GTACACACCGCCCGTC \
--p-r-primer TGATCCTTCTGCAGGTTCACCTAC \
--o-reads ref-seqs

#use the extracted ref reads/tax files to create trained classifier 
nohup qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier classifier.qza

#use the classifier to taxonomically classify the ASVs
qiime feature-classifier classify-sklearn \
--i-classifier classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

#visualise the taxonomy
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

#make a folder to keep our exported files
mkdir export

#MAFFT alignment
time qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment export/aligned-rep-seqs.qza \
--o-masked-alignment export/masked-aligned-rep-seqs.qza \
--o-tree export/unrooted-tree.qza \
--o-rooted-tree export/rooted-tree.qza

#export the ASV table into biom format
qiime tools export \
--input-path table.qza \
--output-path export

#convert the biom format into tsv
biom convert \
-i export/feature-table.biom \
-o export/table.tsv --to-tsv

#export representative sequences
qiime tools export \
--input-path rep-seqs.qza \
--output-path export

#export taxonomy file
qiime tools export \
--input-path taxonomy.qza \
--output-path export

#export the phylogenetic tree
qiime tools export \
--input-path rooted-tree.qza \
--output-path export

#deactivate conda environment
conda deactivate