#!/usr/bin/bash

# install R packages
Rscript install_packages.R

# install python packages
conda create --name fragility_model --file requirements.txt -y