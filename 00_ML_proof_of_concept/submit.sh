#!/usr/bin/bash

RNAfold_path=$1

# get features
bash run_getfeatures.sh $RNAfold_path

# run crash tests
bash run_crashtests.sh

# optimise models
bash run_optimal_models.sh