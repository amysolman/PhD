#!/usr/bin/env python3 

#Script converting qPCR results to gene copies per mL/g sample.

#Python script to:
#1) Open qPCR results and sample data files
#2) Extract mean quantity (i.e. mean genes copies per uL)
#Exmaple: 500 gene copies per uL
#3) Calculate total ng DNA added to the reaction (DNA concentration * 4uL)
#Exmaple: 2ng/uL * 4uL = 8 ng DNA
#4) Divide mean gene copies by total ng DNA in reaction to get gene copies per ng DNA
#Example: 500 / 8 = 62.5 gene copies per ng DNA
#5) Calculate the total ng DNA in the DNA extract elute (remember to multiply by 10 if the sample was diluted or by 2 if a double extraction occured)
#Example: 2ng/uL * 10 * 100uL elute = 2000 ng DNA in total extract elute
#6) Multiply gene copies per ng by total ng in the DNA extract elute
#Example: 62.5 * 2000 = 125000 gene copies in total DNA extract elute
#7) Divide by the weight or volume filtered to get gene copies per mL or g
#Example: 125000 ng DNA / 100 mL = 1250 gene copies per mL
#8) Plot the results

__appname__ = 'PMA-qPCR.py'

__version__ = '0.0.1'

#import libraries 
import pandas as pd
import glob

#1) Open qPCR results and sample data files

#qPCR results
files = glob.glob("../data/*.csv")
files
list = []
for f in files:
    df = pd.read_csv(f, skiprows=7)
    list.append(df)

#sample data file with volumes and weights
samp_data = pd.read_csv("../data/metadata.csv")

#2) Extract mean quantity (i.e. mean genes copies per uL)
def editdfs(df):
    df = df.loc[df["Task"] == "UNKNOWN"]
    df = df[["Sample Name", "Quantity"]]
    return(df)

list2 = [editdfs(df) for df in list]