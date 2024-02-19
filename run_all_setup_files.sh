#!/usr/bin/bash

# directory placehodlers
mkdir -p ./COSMIC/data/{annotations,BaseTable,COSMIC,experiments,kmertone}
mkdir -p ./data/{experiments,human_viruses,kmertone,liftover,models,parameters,range_effects,ranges,ref,setup,TFBS}

cd ./setup/

# install relevant software
bash install_packages.sh

# download datasets for ML demonstrations
bash get_MLdemo_datasets.sh

# download and concatenate DNAfrAIlib feature library
bash get_feature_lib.sh