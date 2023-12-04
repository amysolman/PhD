#!/usr/bin/env python3 

#Script converting qPCR results to gene copies per mL/g sample.

#Python script to:
#1) Open qPCR results and sample data files
#2) Extract mean quantity (i.e. mean genes copies per uL)
#Exmaple: 10,000 gene copies per uL
#3) Multiply estimated gene copies per uL by the number of uL added to the reaction.
#Example: 10,000 * 4 = 40,000 = total gene copies per reaction
#4) Divide by the total ng DNA added to the reaction
#Example: 40,000 / 8 = 5,000 = gene copies per ng DNA
#5) Multiple by original sample DNA concentration (remember to multiply by 10 if a 10-fold dilution occured)
#Example: 5,000 * 2 = 10,000 = gene copies per uL elute 
#6) Multiply by total elute volume (remember this was double for samples that had a double extraction)
#Example: 10,000 * 100 = 1,000,000
#7) Divide by the original sample volume or weight
#Example: 1,000,000 / 50 = 20,000 gene copies per mg/mL

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

#3) Multiply estimated gene copies per uL by the number of uL added to the reaction.
merged_df['Gene Copies Per Reaction'] = merged_df['Quantity Mean'] * 4

#4) Divide by the total ng DNA added to the reaction
merged_df['Gene Copies Per Ng DNA'] = merged_df['Gene Copies Per Reaction'] / ( merged_df['DNA Concentration'] * 4 )

#5) Multiple by original sample DNA concentration (remember to multiply by 10 if a 10-fold dilution occured)
merged_df['Gene Copies Per uL Elute'] = np.select([merged_df['Diluted']], [merged_df['Gene Copies Per Ng DNA'] * merged_df['DNA Concentration'] * 10], default=merged_df['Gene Copies Per Ng DNA'] * merged_df['DNA Concentration'])

#6) Multiply by total elute volume (remember this was double for samples that had a double extraction)
merged_df['Total Gene Copies Extracted'] = np.select([merged_df['Double Extraction']], [(merged_df['Gene Copies Per uL Elute'] * 200)], default=merged_df['Gene Copies Per uL Elute'] * 100)

#7) Divide by the original sample volume or weight
merged_df['Gene Copies Per mL/mg'] = merged_df['Total Gene Copies Extracted'] / merged_df['Volume/Weight']

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