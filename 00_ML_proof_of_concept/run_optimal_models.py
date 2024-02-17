import argparse
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from optimise_model import OptunaModel

parser = argparse.ArgumentParser(description="Run ML demo with Optuna.")
parser.add_argument(
    "-name", "--base_name", required=True, 
    type=str, help="Long range kmer size"
)
parser.add_argument(
    "-k", "--kmer_window", required=True, 
    type=int, help="Long range kmer size"
)
parser.add_argument(
    "-model", "--model_type", required=True, 
    type=str, help="'LM' for Logistic Classifier or 'LGBM' for LightGBM"
)
parser.add_argument(
    "-rc", "--rm_cols", required=False, default=None,
    type=str, help="Cols to remove to avoid potential data leakage. E.g.: 'K562_cells_Top2'"
)
args = parser.parse_args()

base_name = args.base_name
only_breaks = False
crash_test = False
cols_to_remove = args.rm_cols
model_type = args.model_type
num_features = 50

model = OptunaModel(
    base_name=base_name, 
    base_dir='ML_demo',
    kmer_window=args.kmer_window, 
    crash_test=crash_test, 
    only_breaks=only_breaks, 
    model_type=model_type,
    n_trials=50,
    final_model=False,
    cols_to_remove=cols_to_remove
)
model.run_optuna_study()
model.plot_results(num_features=num_features)