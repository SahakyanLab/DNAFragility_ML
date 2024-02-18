#!/usr/bin/bash

pwd="$(pwd)/"

# setup folder structure
mkdir -p ../data/setup/{K562_DMSO_DSBs,K562_Top2_mediated_DSBs}/kmertone
mkdir -p ../data/experiments/{K562_DMSO_DSBs,K562_Top2_mediated_DSBs}/{breakpoint_positions,control_coordinates}

# download DMSO-treated DNA fragility dataset
wget -P ../data/setup/K562_DMSO_DSBs/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444988/suppl/GSM3444988%5FK562%5FDMSO.bed.gz
gunzip ../data/setup/K562_DMSO_DSBs/GSM3444988_K562_DMSO.bed.gz

# download etoposide-treated DNA fragility dataset
wget -P ../data/setup/K562_Top2_mediated_DSBs/ https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3444nnn/GSM3444989/suppl/GSM3444989%5FK562%5FETO.bed.gz
gunzip -f ../data/setup/K562_Top2_mediated_DSBs/GSM3444989_K562_ETO.bed.gz

# process into the appropriate format for kmertone processing and run kmertone
for file in $(find ../data/ranges -name "*.gz")
do
    echo $file
    if [ -f "$file" ]
    then
        gunzip "$file"
    fi
done
Rscript process_MLdemo_datasets.R $pwd

# copy files over to appropriate directories
for name in K562_DMSO_DSBs K562_Top2_mediated_DSBs
do
    cp ../data/setup/$name/kmertone_scores/control_coordinates/* ../data/experiments/$name/control_coordinates/
done