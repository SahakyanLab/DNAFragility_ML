#!/usr/bin/bash

# install R packages
Rscript install_packages.R

# install python packages
conda create --name myenv --file requirements.txt -y