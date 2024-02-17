#!/usr/bin/bash

pwd="$(pwd)/"

# get cancer-associated breakpoints
Rscript 00_GetBaseTable.R $pwd
Rscript 00_GetSVs.R $pwd

# analyse COSMIC breakpoints
Rscript 02_BreakagePersistence.R $pwd
Rscript 03_GapsBetweenBreaks.R $pwd

# extract genomic features
Rscript 02_GetGenomicFeatures.R $pwd

# get negative control breakpoints with kmertone
Rscript 03_ProcessCOSMICBreaks.R $pwd

# copy files to fragility folders
mkdir -p ../../data/experiments/COSMIC/control_coordinates
cp ../data/kmertone/COSMIC/kmertone_scores/control_coordinates/* ../../data/experiments/COSMIC/control_coordinates/

# get single nucleotide positions of human genome
Rscript 02_Bin_breaks_full_genome.R 1 TRUE 1 FALSE $pwd