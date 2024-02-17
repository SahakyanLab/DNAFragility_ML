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
    "-chr", "--chromosome", required=True, 
    type=str, help="Processing chromosome number"
)
args = parser.parse_args()

# global variables
path_to_data = '../data/models/python/lightgbm/chr' + args.chromosome
dest_path = (
    '/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility' + 
    '/data/models/python/lightgbm/chr' + args.chromosome
)
custom_result_prefix = None

model = DNAFragilityModel(
    path_to_data=path_to_data,
    path_to_save=dest_path,
    custom_result_prefix=custom_result_prefix,
    chr_num=args.chromosome,
    feature_matrix_pattern='FullGenome_FeatureMatrix_1_file_chunk_',
    positions_pattern='FullGenome_Ranges_1_file_chunk_',
    model_file='best_LGBM_model.txt'
)
model.run_model()