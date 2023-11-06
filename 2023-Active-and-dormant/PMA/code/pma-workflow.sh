#!/bin/bash

#prep the qPCR data for analysis
python PMA-qPCR.py

#plot and test for significant differences
Rscript qPCR-plot.R