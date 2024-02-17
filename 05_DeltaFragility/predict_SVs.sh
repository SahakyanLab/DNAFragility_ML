#!/usr/bin/bash

source activate myenv

# predict on SNP data
python predict_SVs.py --SV_type SNP

# predict on SV data
python predict_SVs.py --SV_type SV

conda deactivate