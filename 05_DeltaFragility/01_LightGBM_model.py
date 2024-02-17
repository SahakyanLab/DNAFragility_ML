import argparse
import os
import sys
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.metrics import roc_curve
import time
import warnings
import concurrent.futures
warnings.filterwarnings("ignore")

class DNAFragilityModel:
    def __init__(
            self, pwd, directory_path, subdir_name, 
            feature_matrix_pattern, positions_pattern,
            before_or_after,
            model_file='best_LGBM_model.txt', 
            chr_num=None
        ):
        self.chr_num = chr_num
        self.before_or_after = before_or_after
        self.features = None
        self.df_positions = None
        self.pwd = pwd
        self.directory_path = directory_path
        self.subdir_name = subdir_name
        self.feature_matrix_pattern = feature_matrix_pattern
        self.positions_pattern = positions_pattern
        
        path_to_model = '/media/hert6114/Paddy_5TB/ProjectBoard_Patrick/04_DNAFragility/data/models/python/lightgbm'
        self.model = lgb.Booster(
            model_file=f"{path_to_model}/{model_file}"
        )
        self.thresholds_df = pd.read_csv(f"{path_to_model}/best_LGBM_model_thresholds.csv")
    
    def _load_parquet(self, file_path):
        return pd.read_parquet(file_path, use_threads=True)

    def _load_csv(self, file_path):
        return pd.read_csv(file_path)
    
    def get_total_chunks(self):
        file_chunks = [
            f for f in os.listdir(self.subdir_name)
            if f.startswith(self.feature_matrix_pattern) and f.endswith('.parquet')
        ]
        return len(file_chunks)

    def get_data(self, chunk_range):
        t1 = time.time()
        cur_msg = "Loading in all data files"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        # Find all file chunks for feature matrix
        file_chunks = sorted([
            os.path.join(self.subdir_name, f) 
            for f in os.listdir(self.subdir_name) 
            if f.startswith(self.feature_matrix_pattern) and f.endswith('.parquet')
        ])[chunk_range[0]:chunk_range[1]]
        
        # Load and concatenate file chunks in parallel for feature matrix
        with concurrent.futures.ThreadPoolExecutor() as executor:
            data_frames = list(executor.map(self._load_parquet, file_chunks))
        df = pd.concat(data_frames, ignore_index=True)

        # Find all file chunks for positions data
        position_file_chunks = sorted([
            os.path.join(self.subdir_name, f) 
            for f in os.listdir(self.subdir_name) 
            if f.startswith(self.positions_pattern) and f.endswith('.parquet')
        ])[chunk_range[0]:chunk_range[1]]
        
        # Load and concatenate file chunks in parallel for positions data
        with concurrent.futures.ThreadPoolExecutor() as executor:
            position_frames = list(executor.map(self._load_parquet, position_file_chunks))
        self.df_positions = pd.concat(position_frames, ignore_index=True)
        
        self.features = df.columns.tolist()
        X = df.to_numpy(dtype=np.float32)
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
        return X
    
    @staticmethod
    def find_threshold_for_fpr(fpr, tpr, thresholds, target_fpr=0.2):
        idx = np.argmin(np.abs(fpr - target_fpr))
        return thresholds[idx], fpr[idx], tpr[idx]
    
    def process_data(self):
        self.df_positions.drop(columns=['seqnames', 'width'], inplace=True)
        
    def predict(self, X, chunk_idx, target_fpr=[0.001, 0.005, 0.01, 0.05]):
        t1 = time.time()
        cur_msg = "Classifying the presence or absence of breaks"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        y_probs = self.model.predict(X, device='gpu')

        for target_rate in target_fpr:
            threshold = self.thresholds_df[self.thresholds_df['FPR'] == target_rate]['threshold'].iloc[0]
            y_pred = np.where(y_probs > threshold, 1, 0)
            self.df_positions[f'Pred_{target_rate}'] = y_pred
            
        self.df_positions.to_parquet(
            f'{self.subdir_name}/{self.before_or_after}_result_chunk_{chunk_idx}.parquet',
            engine='pyarrow'
        )
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
    def save_results(self):
        t1 = time.time()
        cur_msg = "Saving all prediction results"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        result_chunks = sorted([
            os.path.join(self.subdir_name, f)
            for f in os.listdir(self.subdir_name)
            if f.startswith(f'{self.before_or_after}_result_chunk_') and f.endswith('.parquet')
        ])

        with concurrent.futures.ThreadPoolExecutor() as executor:
            result_frames = list(executor.map(self._load_parquet, result_chunks))
        df_results = pd.concat(result_frames, ignore_index=True)
        
        df_results.to_parquet(
            # f'{self.path_to_data}/final_pred_chr_{self.chr_num}.parquet',
            f'{self.subdir_name}/{self.before_or_after}_final_pred.parquet',
            engine='pyarrow'
        )
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
    def run_model(self):
        total_chunks = self.get_total_chunks()
        chunk_size = 2
        total_iters = range(0, total_chunks, chunk_size)
        num_iters = len(total_iters)

        for idx, i in enumerate(total_iters, start=1):
            print(f"Chunk {idx} of {num_iters}.")
            chunk_range = (i, min(i + chunk_size, total_chunks))
            X_test = self.get_data(chunk_range=chunk_range)
            self.process_data()
            self.predict(X_test, chunk_idx=i // chunk_size)

        self.save_results()
    
if __name__ == "__main__":
    # Iterate over all subdirectories in the specified path
    # directory_path = "/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility/data/deltafragility/"
    directory_path = "/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility/data/deltafragility_SNPs/"
    subdirectories = [os.path.join(directory_path, d) for d in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, d))]

    for i, subdir in enumerate(subdirectories):
        print(f"Processing sequence {i+1}/{len(subdirectories)}: {os.path.basename(subdir)}")
        before_model = DNAFragilityModel(
            pwd="/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility/",
            directory_path=directory_path,
            subdir_name=subdir,
            feature_matrix_pattern='before_feature_matrix',
            positions_pattern='before_positions',
            before_or_after='before',
            model_file='best_LGBM_model.txt'
        )
        before_model.run_model()

        after_model = DNAFragilityModel(
            pwd="/media/hert6114/Paddy_6TB/ProjectBoard_Patrick/04_DNAFragility/",
            directory_path=directory_path,
            subdir_name=subdir,
            feature_matrix_pattern='after_feature_matrix',
            positions_pattern='after_positions',
            before_or_after='after',
            model_file='best_LGBM_model.txt'
        )
        after_model.run_model()