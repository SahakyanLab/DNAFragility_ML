#!/usr/bin/bash

pwd="$(pwd)/"
RNAfold_path=$1

# get 3-mer based crash test data
Rscript 01_model_feat_demo.R K562_DMSO_DSBs 3 TRUE TRUE $pwd $RNAfold_path
Rscript 01_model_feat_demo.R K562_Top2_mediated_DSBs 3 TRUE TRUE $pwd $RNAfold_path

# get 5-mer based crash test data
Rscript 01_model_feat_demo.R K562_Top2_mediated_DSBs 5 TRUE TRUE $pwd $RNAfold_path

# get 3-mer based complete data
Rscript 01_model_feat_demo.R K562_DMSO_DSBs 3 FALSE FALSE $pwd $RNAfold_path
Rscript 01_model_feat_demo.R K562_Top2_mediated_DSBs 3 FALSE FALSE $pwd $RNAfold_path

# get 5-mer based complete data
Rscript 01_model_feat_demo.R K562_Top2_mediated_DSBs 5 FALSE FALSE $pwd $RNAfold_path

# copy data over as example
for file in ../data/ML_demo/K562_DMSO_DSBs_kmerwindow-3crash_test-FALSEonly_breaks-FALSEkmer-8_seed-41964_zscore_{features.parquet,colnames.csv,train.parquet,test.parquet}
do
    cp "$file" "${file//DSBs/DSBs_noK562rm}"
done