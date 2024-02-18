#!/usr/bin/bash

cd ./setup/

# install relevant software
bash install_packages.sh

# download datasets for ML demonstrations
bash get_MLdemo_datasets.sh

# download and concatenate DNAfrAIlib feature library
bash get_feature_lib.sh