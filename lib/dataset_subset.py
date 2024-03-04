import pandas as pd
import numpy as np
from pyfaidx import Fasta
import torch
from torch.utils.data import Dataset
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)
from utils import torch_fromstring

from utils import *

class DatasetSubset(Dataset):
    def __init__(self, ranges, sequence, features_path=None, device='cpu'):
        self.ranges = pd.read_csv(ranges)
        # self.sequence = py2bit.open(sequence)
        self.sequence = Fasta(sequence)
        self.features = None
        self.ranges['Breaks'] = np.log2(1 + self.ranges['Breaks'])
        self.range_width = self.ranges.iloc[0]['end'] - self.ranges.iloc[0]['start']
        self.device = device
        
        if features_path is not None:
            self.features = pd.read_parquet(features_path, use_threads=True)
            self.features.drop(columns=['predictor'], inplace=True)

    # extract sequence from 2bit file
    def extract_sequence(self, seqname, start, end):
        # return self.sequence.sequence(seqname, int(start), int(end)) # 2bit
        return self.sequence[seqname][start:end].seq # fasta

    def one_hot_encode(self, sequence, fixed_length):
        # pad or truncate the sequence to the fixed length
        sequence = sequence.ljust(fixed_length, '.')[:fixed_length]
        # convert sequence to one-hot encoding using pre-computed embeddings
        seq_chrs = torch_fromstring(sequence)
        return one_hot_embed[seq_chrs.long()]

    def __len__(self):
        return len(self.ranges)
    
    @property
    def seq_length(self):
        return self.range_width
    
    def __getitem__(self, idx):
        row = self.ranges.iloc[idx]
        sequence = self.extract_sequence(row['seqnames'], row['start'], row['end'])
        one_hot = self.one_hot_encode(sequence, fixed_length=self.range_width)
        breaks = torch.tensor(row['Breaks'], dtype=torch.float32)
        
        # Move data to GPU if available
        one_hot = one_hot.to(self.device)
        breaks = breaks.to(self.device)
        
        # get engineered feature matrix if available
        if self.features is not None:
            feature_matrix = self.features.iloc[idx]
            feature_matrix = torch.tensor(
                feature_matrix.values, 
                dtype=torch.float32
            ).to(self.device)
            return one_hot, breaks, feature_matrix
        else:
            return one_hot, breaks