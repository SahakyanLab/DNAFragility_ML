import argparse
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from fragility_model import DNAFragilityModel

parser = argparse.ArgumentParser(
    description="Predict DNA fragility at single-base resolution."
)
parser.add_argument(
    "-SV", "--SV_type", required=True, 
    type=str, help="Processing 'SV' or 'SNP'"
)
args = parser.parse_args()

# global variables
dir_path = '/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility/data'
if args.SV_type == 'SV':
    dir_path = dir_path + '/deltafragility'
elif args.SV_type == 'SNP':
    dir_path = dir_path + '/deltafragility_SNPs'

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
        custom_result_prefix='before',
        chr_num=None,
        feature_matrix_pattern='before_feature_matrix',
        positions_pattern='before_positions',
        model_file='best_LGBM_model.txt'
    )
    model.run_model()
    
    model = DNAFragilityModel(
        path_to_data=dest_path,
        path_to_save=dest_path,
        custom_result_prefix='after',
        chr_num=None,
        feature_matrix_pattern='after_feature_matrix',
        positions_pattern='after_positions',
        model_file='best_LGBM_model.txt'
    )
    model.run_model()