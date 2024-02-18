#!/usr/bin/bash

# clone repo
git clone https://github.com/SahakyanLab/DNAfrAIlib.git

# get feature libraries
cd ./DNAfrAIlib/
bash get_tables.sh
cd ../

# copy feature libraries to appropriate folders
mkdir -p ../data/kmertone/QueryTable/
cp ./DNAfrAIlib/QueryTable_kmer* ../data/kmertone/QueryTable/