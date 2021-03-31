#Bash script for downloading raw sequences and calling QIIME2/DADA2 to process them.
#Amy Solman
#04-03-2021

mkdir Data
cd Data
mkdir Sommers2020
cd Sommers2020

conda install -c bioconda sra-tools

for (( i = 327; i <= 408; i++ ))
  do
  prefetch SRR12805$i
done

for (( i = 327; i <= 408; i++ ))
  do
  cd SRR12805$i
  fastq-dump SRR12805$i.sra --split-files --origfmt
  cd ..
done

cd Documents/GitHub/PhD/Projects/Practise/Data-wrangling/
mkdir fastq-sommers2020
cp Sommers2020/SRR12805*/*.fastq fastq-sommers2020/
cp Sommers2020/SraRunTable.txt fastq-sommers2020/

cd fastq-sommers2020
mkdir qiime_import
mv *.fastq qiime_import
cd qiime_import
rename 's/_/_00_L001_R/' *
rename 's/.fastq/_001.fastq/' *
gzip *.fastq

cd../..
Rscript 1b-metadata-wrangle.R

mkdir work

cd ..
conda activate qiime2-2021.2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path fastq-sommers2020/qiime_import_16S \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path work/demux-paired-end.qza
  
  qiime cutadapt trim-paired \
--i-demultiplexed-sequences work/demux-paired-end.qza \
--p-cores 16 \
--p-front-f GTGCCAGCMGCCGCGGTAA \
--p-adapter-f ATTAGAWACCCBDGTAGTCC \
--p-front-r GGACTACHVGGGTWTCTAAT \
--p-adapter-r TTACCGCGGCKGCTGGCAC \
--p-minimum-length 160 \
--o-trimmed-sequences work/demux-paired-end-trimmed.qza \
--verbose 

qiime demux summarize \
--i-data work/demux-paired-end-trimmed.qza \
--o-visualization work/demux-paired-end-trimmed.qzv

###VISUALIZE###
qiime tools view work/demux-paired-end-trimmed.qzv

###CALL DADA2###

qiime dada2 denoise-paired \
--i-demultiplexed-seqs work/demux-paired-end-trimmed.qza \
--p-trunc-len-f 200 \
--p-trunc-len-r 160 \
--p-n-threads 0 \
--output-dir work/DADA2_denoising_output \
--verbose

###LETS LOOK AT THE RESULTS###
qiime metadata tabulate \
  --m-input-file work/DADA2_denoising_output/denoising_stats.qza \
  --o-visualization work/DADA2_denoising_output/denoising_stats.qzv
  
  qiime tools view work/DADA2_denoising_output/denoising_stats.qzv
  
###SUMMARY OF THE REMAINING SEQUENCE READS###

  qiime feature-table tabulate-seqs \
--i-data work/DADA2_denoising_output/representative_sequences.qza \
--o-visualization work/DADA2_denoising_output/rep_seqs.qzv 

###Through the sequence table generated below you can click on a representative sequence and BLAST it
qiime tools view work/DADA2_denoising_output/rep_seqs.qzv

###FEATURE TABLE###

qiime feature-table summarize \
--i-table work/DADA2_denoising_output/table.qza \
--o-visualization work/DADA2_denoising_output/table.qzv 

qiime tools view work/DADA2_denoising_output/table.qzv

###ASSIGN TAXONOMY###

###TO DO THIS USING SILVA WE HAVE TO INSTALL AND ACTIVATE A DIFFERENT VERSION OF QIIME2

conda deactivate 

wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-osx-conda.yml
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-osx-conda.yml
rm qiime2-2019.10-py36-osx-conda.yml
conda activate qiime2-2019.10
conda deactivate 

 wget https://data.qiime2.org/2019.10/common/silva-132-99-515-806-nb-classifier.qza
mv silva-132-99-515-806-nb-classifier.qza work/

conda activate qiime2-2019.10

####WARNING: THE BELOW COMMAND TAKES AN HOUR TO RUN#### 

qiime feature-classifier classify-sklearn \
--i-classifier work/silva-132-99-515-806-nb-classifier.qza \
--i-reads work/DADA2_denoising_output/representative_sequences.qza \
--output-dir work/classified_sequences \
--verbose

# and visualize
qiime metadata tabulate \
  --m-input-file work/classified_sequences/classification.qza \
  --o-visualization work/classified_sequences/taxonomy.qzv
  
  qiime tools view work/classified_sequences/taxonomy.qzv
  
###LETS MAKE A PHYLOGENETIC TREE!###

mkdir work/phylogeny

qiime alignment mafft \
  --i-sequences work/DADA2_denoising_output/representative_sequences.qza \
  --o-alignment work/phylogeny/aligned-rep-seqs.qza
  
###MASK THE ALIGNMENT###

qiime alignment mask \
  --i-alignment work/phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment work/phylogeny/masked-aligned-rep-seqs.qza
  
###CONSTRUCT THE PHYLOGENY###

qiime phylogeny fasttree \
  --i-alignment work/phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree work/phylogeny/fasttree-tree.qza
  
###ROOT THE TREE SO IT CAN BE USED IN UNIFRAC###

qiime phylogeny midpoint-root \
  --i-tree work/phylogeny/fasttree-tree.qza \
  --o-rooted-tree work/phylogeny/fasttree-tree-rooted.qza

###LET'S EXPORT OUR WORK

mkdir work/export

qiime tools export \
  --input-path work/DADA2_denoising_output/table.qza \
  --output-path work/export/table
  
biom convert \
-i work/export/table/feature-table.biom \
-o work/export/table/table.tsv --to-tsv

qiime tools export \
  --input-path work/DADA2_denoising_output/representative_sequences.qza \
  --output-path work/export/rep-seqs.fasta
  
qiime tools export \
  --input-path work/classified_sequences/classification.qza \
  --output-path work/export/taxonomy
  
qiime tools export \
  --input-path work/phylogeny/fasttree-tree-rooted.qza \
  --output-path work/export/exported-tree
  
conda deactivate