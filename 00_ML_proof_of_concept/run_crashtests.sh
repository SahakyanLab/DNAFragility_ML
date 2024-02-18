#!/usr/bin/bash

source activate fragility_model

# optimise crash tests using 3-mer data
python run_crashtest.py --base_name K562_DMSO_DSBs --kmer_window 3
python run_crashtest.py --base_name K562_Top2_mediated_DSBs --kmer_window 3

# optimise crash tests using 5-mer data
python run_crashtest.py --base_name K562_Top2_mediated_DSBs --kmer_window 5

conda deactivate