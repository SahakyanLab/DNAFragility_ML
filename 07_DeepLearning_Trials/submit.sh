#!/usr/bin/bash

pwd="$(pwd)/"
bin_width=20000

# for classification with features
Rscript 00_Extract_Features_for_Model.R $pwd FALSE $bin_width

# for regression with long-context range
Rscript 00_Extract_Features_for_Model.R $pwd TRUE $bin_width

# run various configurations of the deep DNA fragility model
source activate hyena-dna

# regression
epochs=100
python run_deep_fragility.py --regression --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --regression --with_LSTM --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --regression --with_featmatrix --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --regression --with_featmatrix --with_LSTM --epochs $epochs --bin_width $bin_width

# classification
epochs=20
python run_deep_fragility.py --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --with_LSTM --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --with_featmatrix --epochs $epochs --bin_width $bin_width
python run_deep_fragility.py --with_featmatrix --with_LSTM --epochs $epochs --bin_width $bin_width

conda deactivate