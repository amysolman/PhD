#!/usr/bin/env python3 

#Script converting qPCR results to gene copies per mL/g sample.

#Python script to:
#1) Open qPCR results and sample data files
#2) Extract mean quantity (i.e. mean genes copies per uL)
#Exmaple: 500 gene copies per uL
#3) Calculate total ng DNA added to the reaction (DNA concentration * 4uL). Remember to divide by 10 for those that were diluted
#Exmaple: 2ng/uL * 4uL = 8 ng DNA
#4) Divide mean gene copies by total ng DNA in reaction to get gene copies per ng DNA
#Example: 500 / 8 = 62.5 gene copies per ng DNA
#5) Calculate the total ng DNA in the DNA extract elute (remember to multiply by 2 if a double extraction occured)
#Example: 2ng/uL * 100uL elute = 200 ng DNA in total extract elute
#6) Multiply gene copies per ng by total ng in the DNA extract elute
#Example: 62.5 * 2000 = 125000 gene copies in total DNA extract elute
#7) Divide by the weight or volume filtered to get gene copies per mL or g
#Example: 125000 ng DNA / 100 mL = 1250 gene copies per mL

__appname__ = 'PMA-qPCR.py'

__version__ = '0.0.1'

#import libraries 
import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt

#1) Open qPCR results and extract mean quantity (i.e. mean genes copies per uL) by removing unwanted rows/columns

#get our dataframes
files = glob.glob("../data/*results.csv")
files
list = []
for f in files:
    df = pd.read_csv(f, skiprows=7)
    list.append(df)

#function for editing dataframes
def editdfs(df):
    df = df.loc[df["Task"] == "UNKNOWN"]
    df = df[["Sample Name", "Quantity Mean"]]
    return(df)

#run function and put output into new list
list2 = [editdfs(df) for df in list]

#combine our qPCR results dataframes 
qpcr_df = pd.concat(list2)

#keep only unique rows
qpcr_df2 = qpcr_df.drop_duplicates()

#2) Get metadata, replace "Below detection" values and combine with qPCR results

#sample data file with volumes and weights
samp_data = pd.read_csv("../data/metadata.csv")

#replace "Below detection" values with half the lower limit of the HS Qubit Kit range (0.005/2 = 0.0025 ng/uL)
samp_data = samp_data.replace('Below detection', 0.005/2)

#make sure our DNA concentrations are numeric
samp_data["DNA Concentration"] = pd.to_numeric(samp_data["DNA Concentration"])
samp_data["DNA Concentration 1"] = pd.to_numeric(samp_data["DNA Concentration 1"])
samp_data["DNA Concentration 2"] = pd.to_numeric(samp_data["DNA Concentration 2"])

#make the sample name columns the same format
samp_data['Sample Name'] = samp_data['Sample Name'].str.strip()
qpcr_df2['Sample Name'] = qpcr_df2['Sample Name'].str.strip()

#merge the sample data and qPCR results into one dataframe
merged_df = samp_data.merge(qpcr_df2, how = 'inner', on = ['Sample Name'])

#3) Multiply gene copies of 10-fold diluted samples by 10

merged_df['Quantity Mean Edit'] = np.select([merged_df['Diluted']], [merged_df["Quantity Mean"] * 10], default=merged_df["Quantity Mean"])

#3) Calculate gene copiers per ng DNA

merged_df['Gene Copies Per Ng'] = merged_df['Quantity Mean Edit'] / merged_df['DNA Concentration']

#4) Calculate the total ng DNA in the DNA extract elute

#for samples that underwent a double extraction, work out the total ng DNA per extraction and add them together,
#otherwise multiply the DNA concentration by 100uL (the volume of the elute)
merged_df['Total Ng Extracted'] = np.select([merged_df['Double Extraction']], [(merged_df["DNA Concentration 1"] * 100) +  (merged_df["DNA Concentration 2"] * 100)], default=merged_df["DNA Concentration"] * 100)

#6) Multiply gene copies per ng by total ng in the DNA extract elute to get total gene copies extracted from each filter/sediment
merged_df['Gene Copies Extracted'] = merged_df['Gene Copies Per Ng'] * merged_df['Total Ng Extracted']

#7) Divide by the weight or volume filtered to get gene copies per mL or mg
merged_df['Gene Copies Per mL/mg'] = merged_df['Gene Copies Extracted'] / merged_df['Volume/Weight']

#save the edited dataframe
merged_df.to_csv("../data/meta_qpcr_res.csv", sep='\t')

#calculate the mean values per sample (as each environmental sample was processed as an experimental double)
#merge into a GroupBy object
new_df = merged_df[['Sample', 'Treatment', 'Habitat', 'Gene Copies Per mL/mg']]
grouped = new_df.groupby(['Sample', 'Treatment', 'Habitat'])
mean_genes = grouped.mean()
mean_genes

#export
mean_genes.to_csv("../data/qPCR-mean.csv", sep='\t')