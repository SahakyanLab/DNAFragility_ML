import os
import sys
from optimise_model import OptunaModel

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

base_name = 'models/python/lightgbm'
only_breaks = False
crash_test = False
cols_to_remove = None
model_type = 'LGBM'
num_features = 25

model = OptunaModel(
    base_name=base_name, 
    base_dir=base_name,
    kmer_window=3, 
    crash_test=crash_test, 
    only_breaks=only_breaks, 
    model_type=model_type,
    n_trials=100,
    final_model=True,
    cols_to_remove=cols_to_remove
)
model.run_optuna_study()
model.plot_results_of_final_model(num_features=num_features)

'../data/models/python/lightgbm/best_LGBM_model_studytrials.csv'