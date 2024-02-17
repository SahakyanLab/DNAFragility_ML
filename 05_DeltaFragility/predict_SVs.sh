#!/usr/bin/bash

source /home/imm/hert6114/anaconda3/bin/activate myenv

# predict on SNP data
python predict_SVs.py --SV_type SNP

# predict on SV data
python predict_SVs.py --SV_type SV

conda deactivate