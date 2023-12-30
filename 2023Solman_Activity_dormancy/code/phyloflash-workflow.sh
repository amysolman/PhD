#!/bin/bash

# Phyloflash workflow for totalRNA

#install phyloflash via mamba
mamba create -n pf phyloflash
 
#activate the environment
mamba activate pf 

#check prerequisits
phyloFlash.pl -check_env #all found 

#navigate to our metatranscriptomics work folder
cd mbl_2022-mt

#make and enter a folder for phyloflash analysis
mkdir phyloflash
cd phyloflash

#Download pre-formated database - SILVA 138.1
wget https://zenodo.org/record/7892522/files/138.1.tar.gz

#unpack the database
tar -xzf 138.1.tar.gz

#run phyloflash on trimmed files
#readlength = expected read length (default is 100)

for file in ../Fastq-trimmed/*forward_paired.fastq.gz; 
do 
R1="${file:17:6}";
echo "$R1"
phyloFlash.pl -lib "$R1" \
-read1 ../Fastq-trimmed/"$R1"_forward_paired.fastq.gz \
-read2 ../Fastq-trimmed/"$R1"_reverse_paired.fastq.gz \
-readlength 250
done