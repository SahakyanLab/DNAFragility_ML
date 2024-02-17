#!/usr/bin/bash

source /home/imm/hert6114/anaconda3/bin/activate myenv

# optimise logistic classifier using 3-mer data
export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_DMSO_DSBs --kmer_window 3 --model_type LM --rm_cols K562_cells_Top2 --rm_cols K562_cells_Top2

export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_DMSO_DSBs_noK562rm --kmer_window 3 --model_type LM

export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_Top2_mediated_DSBs --kmer_window 3 --model_type LM --rm_cols K562_cells_Top2

# optimise LightGBM using 3-mer data
export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_DMSO_DSBs --kmer_window 3 --model_type LGBM --rm_cols K562_cells_Top2

export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_DMSO_DSBs_noK562rm --kmer_window 3 --model_type LGBM

export PYTHONHASHSEED=0
python run_optimal_models.py --base_name K562_Top2_mediated_DSBs --kmer_window 3 --model_type LGBM --rm_cols K562_cells_Top2

# # optimise logistic classifier using 5-mer data
# export PYTHONHASHSEED=0
# python run_optimal_models.py --base_name K562_Top2_mediated_DSBs --kmer_window 5 --model_type LM --rm_cols K562_cells_Top2

# # optimise LightGBM using 5-mer data
# export PYTHONHASHSEED=0
# python run_optimal_models.py --base_name K562_Top2_mediated_DSBs --kmer_window 5 --model_type LGBM --rm_cols K562_cells_Top2

conda deactivate