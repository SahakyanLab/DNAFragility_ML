#!/usr/bin/bash

RNAfold_path=$1
fast_matrix=$2

# get features
bash run_getfeatures.sh $RNAfold_path $fast_matrix

# run crash tests
bash run_crashtests.sh

# optimise models
bash run_optimal_models.sh