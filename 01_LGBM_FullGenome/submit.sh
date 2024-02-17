#!/usr/bin/bash

pwd="$(pwd)/"
RNAfold_path=$1

# for training process
Rscript 01_Extract_Features_for_Model.R TRUE 1 $pwd $RNAfold_path

# optimise lightGBM model
source activate myenv
python run_optimal_models.py
conda deactivate

# extract feature matrix of the full human genome
for chr in {1..22}
do
    Rscript 01_Extract_Features_for_Model.R FALSE $chr $pwd $RNAfold_path
done

# for testing on full genome
source activate myenv
for chr in {1..22}
do
    python3 predict_genome.py -chr $chr
done
conda deactivate

# Get genomic features
Rscript ../../05_Cosmic/scripts/02_GetGenomicFeatures.R $pwd

# Run all remaining analyses
Rscript 02_Overlap_with_genic_feats.R $pwd
Rscript 03_Random_sample_LMH_bins.R $pwd
Rscript 04_Overlap_with_TFBS.R $pwd
Rscript 05_Bins_Overlap_vs_NoOverlap.R $pwd
Rscript 06_Finding_optimal_bin_size.R $pwd
Rscript ../03_Chromothripsis.R $pwd

for bw in 1960 10000 20000
do
    Rscript 03_LMH_binned_breaks.R $bw $pwd
done