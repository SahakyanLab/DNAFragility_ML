#!/usr/bin/bash

# for training process
Rscript 01_Extract_Features_for_Model.R TRUE 1

# optimise lightGBM model
source /home/imm/hert6114/anaconda3/bin/activate myenv
python run_optimal_models.py
conda deactivate

# extract feature matrix of the full human genome
for chr in {1..22}
do
    Rscript 01_Extract_Features_for_Model.R FALSE $chr
done

# for testing on full genome
source /home/imm/hert6114/anaconda3/bin/activate myenv
for chr in {1..22}
do
    python3 predict_genome.py -chr $chr
done
conda deactivate

# Get genomic features
Rscript ../../05_Cosmic/scripts/02_GetGenomicFeatures.R

# Run all remaining analyses
Rscript 02_Overlap_with_genic_feats.R
Rscript 03_Random_sample_LMH_bins.R
Rscript 04_Overlap_with_TFBS.R
Rscript 05_Bins_Overlap_vs_NoOverlap.R
Rscript 06_Finding_optimal_bin_size.R
Rscript ../03_Chromothripsis.R

for bw in 1960 10000 20000
do
    Rscript 03_LMH_binned_breaks.R $bw
done