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
    type=str, help="Processing 'SNP' or anything else"
)
args = parser.parse_args()

# global variables
chr_num = None
if args.SV_type == 'SNP':
    dir_path = '../data/SNP_nullomer_sequence_SNPs'
elif args.SV_type == 'Not_SNP':
    dir_path = '../data/nullomer_and_shuffle'

path_to_data = []
for d in os.listdir(dir_path):
    full_path = os.path.join(dir_path, d)
    if os.path.isdir(full_path):
        path_to_data.append(full_path)

for i, dest_path in enumerate(path_to_data):    
    print(f"Processing sequence {i+1}/{len(path_to_data)}: {os.path.basename(dest_path)}")
    if args.SV_type == 'SNP':
        model = DNAFragilityModel(
            path_to_data=dest_path,
            path_to_save=dest_path,
            custom_result_prefix='Nullomer',
            chr_num=None,
            feature_matrix_pattern='Nullomer_feature_matrix',
            positions_pattern='Nullomer_positions',
            model_file='best_LGBM_model.txt'
        )
        model.run_model()
        
        model = DNAFragilityModel(
            path_to_data=dest_path,
            path_to_save=dest_path,
            custom_result_prefix='Present',
            chr_num=None,
            feature_matrix_pattern='Present_feature_matrix',
            positions_pattern='Present_positions',
            model_file='best_LGBM_model.txt'
        )
    elif args.SV_type == 'Not_SNP':
        model = DNAFragilityModel(
            path_to_data=dest_path,
            path_to_save=dest_path,
            custom_result_prefix=None,
            chr_num=chr_num,
            feature_matrix_pattern='feature_matrix',
            positions_pattern='positions',
            model_file='best_LGBM_model.txt'
        )
        
    model.run_model()