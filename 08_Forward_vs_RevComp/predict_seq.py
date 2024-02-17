import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from fragility_model import DNAFragilityModel

# global variables
dir_path = '../data/Forward_vs_RevComp/'
custom_result_prefix = None
chr_num = None

path_to_data = []
for d in os.listdir(dir_path):
    full_path = os.path.join(dir_path, d)
    if os.path.isdir(full_path):
        path_to_data.append(full_path)

for i, dest_path in enumerate(path_to_data):    
    print(f"Processing sequence {i+1}/{len(path_to_data)}: {os.path.basename(dest_path)}")
    
    model = DNAFragilityModel(
        path_to_data=dest_path,
        path_to_save=dest_path,
        custom_result_prefix=custom_result_prefix,
        chr_num=chr_num,
        feature_matrix_pattern='feature_matrix',
        positions_pattern='positions',
        model_file='best_LGBM_model.txt'
    )
    model.run_model()