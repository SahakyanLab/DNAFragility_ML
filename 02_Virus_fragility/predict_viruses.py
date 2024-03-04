import argparse
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from fragility_model import DNAFragilityModel

parser = argparse.ArgumentParser(
    description="Predict DNA fragility at single-base resolution"
)
parser.add_argument(
    "--random", 
    action=argparse.BooleanOptionalAction, 
    default=False,
    help="Processing true or randomly generated DNA virus sequences"
)
args = parser.parse_args()

# global variables
dir_path = '../data/human_viruses_shuffle/' if args.random else '../data/human_viruses/'
custom_result_prefix = None
chr_num = None

path_to_data = []
for d in os.listdir(dir_path):
    full_path = os.path.join(dir_path, d)
    if os.path.isdir(full_path):
        path_to_data.append(full_path)

for i in range(len(path_to_data)):
    print(f"Processing {i+1}/{len(path_to_data)+1} sequence: {path_to_data[i].split('/')[-1]}")
    
    model = DNAFragilityModel(
        path_to_data=path_to_data[i],
        path_to_save=path_to_data[i],
        custom_result_prefix=custom_result_prefix,
        chr_num=chr_num,
        feature_matrix_pattern='feature_matrix',
        positions_pattern='positions',
        model_file='best_LGBM_model.txt'
    )
    model.run_model()