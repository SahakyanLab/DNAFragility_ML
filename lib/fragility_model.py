import os
import sys
import pandas as pd
import numpy as np
import lightgbm as lgb
import time
import warnings
import concurrent.futures
warnings.filterwarnings("ignore")

class DNAFragilityModel:
    def __init__(
            self, custom_result_prefix, 
            path_to_data, path_to_save,
            feature_matrix_pattern, positions_pattern,
            model_file='best_LGBM_model.txt', 
            chr_num=None
        ):
        self.chr_num = chr_num
        self.path_to_data = path_to_data
        self.path_to_save = path_to_save
        self.feature_matrix_pattern = feature_matrix_pattern
        self.positions_pattern = positions_pattern
        self.custom_result_prefix = custom_result_prefix
        
        path_to_model = '../data/models/python/lightgbm'
        self.model = lgb.Booster(
            model_file= path_to_model + '/' + model_file
        )
        self.thresholds_df = pd.read_csv(
            path_to_model + '/best_LGBM_model_thresholds.csv'
        )
        
        if not os.path.exists(self.path_to_data):
            os.makedirs(self.path_to_data, exist_ok=True)
            
        if not os.path.exists(self.path_to_save):
            os.makedirs(self.path_to_save, exist_ok=True)

    def _load_parquet(self, file_chunks):
        return pd.read_parquet(file_chunks, use_threads=True)
    
    def get_total_chunks(self):
        file_chunks = [
            f for f in os.listdir(self.path_to_data)
            if f.startswith(self.feature_matrix_pattern) and f.endswith('.parquet')
        ]
        return len(file_chunks)

    def get_data(self, chunk_range):
        t1 = time.time()
        cur_msg = "Importing all datasets"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        # Find all file chunks for feature matrix
        file_chunks = sorted([
            os.path.join(self.path_to_data, f) 
            for f in os.listdir(self.path_to_data) 
            if f.startswith(self.feature_matrix_pattern) and f.endswith('.parquet')
        ])[chunk_range[0]:chunk_range[1]]
        
        # Load and concatenate file chunks in parallel for feature matrix
        with concurrent.futures.ThreadPoolExecutor() as executor:
            data_frames = list(executor.map(self._load_parquet, file_chunks))
        df = pd.concat(data_frames, ignore_index=True)

        # Find all file chunks for positions data
        position_file_chunks = sorted([
            os.path.join(self.path_to_data, f) 
            for f in os.listdir(self.path_to_data) 
            if f.startswith(self.positions_pattern) and f.endswith('.parquet')
        ])[chunk_range[0]:chunk_range[1]]
        
        # Load and concatenate file chunks in parallel for positions data
        with concurrent.futures.ThreadPoolExecutor() as executor:
            position_frames = list(executor.map(self._load_parquet, position_file_chunks))
        self.df_positions = pd.concat(position_frames, ignore_index=True)
        
        df.drop(columns=['predictor'], errors='ignore', inplace=True)
        self.features = df.columns.tolist()
        X = df.to_numpy(dtype=np.float32)
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
        return X
    
    @staticmethod
    def find_threshold_for_fpr(fpr, tpr, thresholds, target_fpr=0.001):
        idx = np.argmin(np.abs(fpr - target_fpr))
        return thresholds[idx], fpr[idx], tpr[idx]
    
    def process_target_col(self):
        if 'Breaks' in self.df_positions.columns:
            self.df_positions = self.df_positions[['Breaks']]
            self.df_positions['Breaks'] = self.df_positions['Breaks'].map({'YES': 1, 'NO': 0})        
            self.df_positions.rename(columns={'Breaks': 'True'}, inplace=True)
        else:
            self.df_positions['True'] = 0
            self.df_positions = self.df_positions.filter(['True'])
        
    def predict(self, X, chunk_idx, target_fpr=[0.001, 0.005, 0.01, 0.05]):
        t1 = time.time()
        cur_msg = "Classifying DNA fragility"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        y_probs = self.model.predict(X, device='gpu')

        for target_rate in target_fpr:
            threshold = self.thresholds_df[self.thresholds_df['FPR'] == target_rate]['threshold'].iloc[0]
            y_pred = np.where(y_probs > threshold, 1, 0)
            self.df_positions[f'Pred_{target_rate}'] = y_pred
            
        file_name = (
            self.path_to_save + '/' + 
            (self.custom_result_prefix + '_' if self.custom_result_prefix is not None else '') + 
            'result_chunk_' + str(chunk_idx) + 
            '.parquet'
        )
        self.df_positions.to_parquet(file_name, engine='pyarrow')
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
    def save_results(self):
        t1 = time.time()
        cur_msg = "Saving all prediction results"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        file_name = (
            (self.custom_result_prefix + '_' if self.custom_result_prefix is not None else '') + 
            'result_chunk_'
        )
        result_chunks = sorted([
            os.path.join(self.path_to_save, f)
            for f in os.listdir(self.path_to_save)
            if f.startswith(file_name) and f.endswith('.parquet')
        ])

        with concurrent.futures.ThreadPoolExecutor() as executor:
            result_frames = list(executor.map(self._load_parquet, result_chunks))
        df_results = pd.concat(result_frames, ignore_index=True)
        
        file_name = (
            self.path_to_save + '/' + 
            (self.custom_result_prefix + '_' if self.custom_result_prefix is not None else '') + 
            'final_pred' + 
            ('_chr' + str(self.chr_num) if self.chr_num is not None else '') + 
            '.parquet'
        )
        df_results.to_parquet(file_name, engine='pyarrow')
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
    def run_model(self):
        total_chunks = self.get_total_chunks()
        chunk_size = 2
        total_iters = range(0, total_chunks, chunk_size)
        num_iters = len(total_iters)

        for idx, i in enumerate(total_iters, start=1):
            if num_iters > 1:
                msg = 'Chunk ' + str(idx) + ' of ' + str(num_iters) + '.'
                if self.chr_num:
                    msg = 'Chromosome ' + str(self.chr_num) + '. ' + msg
                print(msg)
                
            chunk_range = (i, min(i + chunk_size, total_chunks))
            X_test = self.get_data(chunk_range=chunk_range)
            self.process_target_col()
            self.predict(X_test, chunk_idx=i // chunk_size)

        self.save_results()