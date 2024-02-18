#!/usr/bin/bash

pwd="$(pwd)"

# get cancer-associated breakpoints
echo "Processing cancer-associated DNA strand breaks..."
Rscript 00_GetBaseTable.R $pwd
Rscript 00_GetSVs.R $pwd

# analyse COSMIC breakpoints
echo "Analysing cancer-associated DNA strand breaks..."
Rscript 02_BreakagePersistence.R $pwd
Rscript 03_GapsBetweenBreaks.R $pwd

# extract genomic features
echo "Extracting genomic features of human genome..."
Rscript 02_GetGenomicFeatures.R $pwd

# get negative control breakpoints with kmertone
echo "Run kmertone to extract k-mer scores..."
Rscript 03_ProcessCOSMICBreaks.R $pwd

# copy files to fragility folders
mkdir -p ../../data/experiments/COSMIC/control_coordinates
cp ../data/kmertone/COSMIC/kmertone_scores/control_coordinates/* ../../data/experiments/COSMIC/control_coordinates/

# # get single nucleotide positions of human genome
echo "Extract single nucleotide positions of human genome..."
Rscript 02_Bin_breaks_full_genome.R 1 TRUE 1 FALSE $pwd