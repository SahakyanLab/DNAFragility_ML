import torch
import torch.nn as nn
import torch.nn.functional as F

from equirc.pytorch_rclayers import RegToRegConv

class CNN(nn.Module):
    def __init__(self, seq_length, use_feature_matrix, with_LSTM, dropout, regression, verbose=False):
        super(CNN, self).__init__()
        
        self.seq_length = seq_length
        self.use_feature_matrix = use_feature_matrix
        self.with_LSTM = with_LSTM
        self.dropout = nn.Dropout(p=dropout)
        self.regression = regression
        self.verbose = verbose
        
        # hidden neurons in the LSTM 
        self.lstm_size = 4
        
        # Define convolutional layers
        self.layers = nn.Sequential(
            RegToRegConv(reg_in=2, reg_out=32, kernel_size=7, padding=3), # 2 * reg_out = 64
            # nn.Conv1d(4, 64, kernel_size=7, stride=1, padding=3),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2),
            
            RegToRegConv(reg_in=32, reg_out=16, kernel_size=7, padding=3),
            # nn.Conv1d(64, 32, kernel_size=7, stride=1, padding=3),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2),
            
            RegToRegConv(reg_in=16, reg_out=8, kernel_size=7, padding=3),
            # nn.Conv1d(32, 16, kernel_size=7, stride=1, padding=3),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2),
                        
            RegToRegConv(reg_in=8, reg_out=4, kernel_size=7, padding=3),
            # nn.Conv1d(16, 8, kernel_size=7, stride=1, padding=3),
            nn.BatchNorm1d(8),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2),
            
            RegToRegConv(reg_in=4, reg_out=2, kernel_size=7, padding=3),
            # nn.Conv1d(8, 4, kernel_size=7, stride=1, padding=3),
            nn.BatchNorm1d(4),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2, stride=2)
        )
        
        # Compute the length of the sequence after the convolutional and pooling layers
        # 5 max pooling layers
        self.conv_outlen = self.seq_length // (2 ** 5)  # 5 max pooling layers
        self.fc_dim = self.lstm_size * self.conv_outlen
        
        # Bi-directional LSTM layer to process the output of the convolutional
        # layers before passing it to the self-attention layer.
        if self.with_LSTM:
            self.lstm = nn.LSTM(
                input_size=4, 
                hidden_size=self.lstm_size, 
                num_layers=1, 
                batch_first=True, 
                bidirectional=True
            )
            self.fc_dim = self.fc_dim - 2 * 2
        
        self.fc_layers = nn.Sequential(
            nn.Linear(self.fc_dim, 64),
            nn.ReLU(),
            self.dropout,
            
            nn.Linear(64, 32),
            nn.ReLU(),
            self.dropout,
        
            nn.Linear(32, 1)
        )
        
    def forward(self, x):
        if self.verbose: 
            print(f'x shape: {x.shape}')
        
        # Pass through convolutional, BN, and pooling layers
        x = self.layers(x)
        if self.verbose: 
            print(f'x shape after layers: {x.shape}')
   
        # Reshape the tensor for LSTM 
        # (batch_size, seq_len, input_size)
        x = x.view(-1, self.conv_outlen, 4)
        if self.verbose: 
            print(f'x shape reshaping for LSTM:  {x.shape}')

        # Pass through LSTM
        if self.with_LSTM:
            x, (hn, cn) = self.lstm(x)
            
            # switch dimensions i maxpool along seq_len dimension
            x = x.permute(0, 2, 1)
            x = F.max_pool1d(x, kernel_size=2, stride=2)
        
        # switch dimensions back
        x = x.permute(0, 2, 1)
        if self.verbose: 
            print(f'x shape after passing LSTM: {x.shape}')
        
        # Reshape back to shape expected by self-attention layer
        # (batch_size, input_size, seq_len)
        x = x.permute(0, 2, 1)
                
        # Determine the current batch size
        current_batch_size = x.size(0)
        
        # Flatten the tensor
        x = x.reshape(current_batch_size, -1)
        if self.verbose: 
            print(f'x shape after flattening tensor: {x.shape}')

        # Pass through fully connected layers        
        x = self.fc_layers(x)
        if self.verbose: 
            print(f'x shape after fully connected layers: {x.shape}')
    
        return x
    
class CNNandFeaturesModel(nn.Module):
    def __init__(self, data, with_LSTM, dropout, regression, verbose=False):
        super(CNNandFeaturesModel, self).__init__()
        
        self.data = data
        self.with_LSTM = with_LSTM
        self.use_feature_matrix = self.data.features is not None
        self.fc_dim = 1 + 32 if self.use_feature_matrix else 1
        self.dropout = nn.Dropout(p=dropout)
        self.regression = regression
        self.verbose = verbose
        
        self.cnn = CNN(
            seq_length=self.data.seq_length, 
            use_feature_matrix=self.use_feature_matrix, 
            with_LSTM=self.with_LSTM,
            dropout=dropout,
            regression=regression
        )
        
        # Define fully connected layers for additional feature matrix
        if self.use_feature_matrix:
            self.num_additional_features = len(self.data.features.columns)
            self.fc_additional1 = nn.Linear(self.num_additional_features, 128)
            self.fc_additional2 = nn.Linear(128, 64)
            self.fc_additional3 = nn.Linear(64, 32)

            # Define fully connected layers post-concatenation
            # 1 from CNN and 32 from additional feature matrix
            self.fc1 = nn.Linear(self.fc_dim, 32)
            self.fc2 = nn.Linear(32, 16)
            self.fc3 = nn.Linear(16, 1)
        
    def forward(self, x, feature_matrix):
        if self.verbose: 
            print(f'X shape: {x.shape}')
        
        # Reshape x to have dimensions (batch_size, 4, sequence_length)
        x = x.view(-1, 4, self.cnn.seq_length)
        if self.verbose: 
            print(f'X shape after reshape: {x.shape}')
        
        final_output = cnn_output = self.cnn(x)
        if self.verbose: 
            print(f'cnn_output shape: {cnn_output.shape}')
        
        if self.use_feature_matrix:
            additional_output = F.relu(self.fc_additional1(feature_matrix))
            additional_output = F.relu(self.fc_additional2(additional_output))
            additional_output = F.relu(self.fc_additional3(additional_output))
            
            combined = torch.cat((cnn_output, additional_output), dim=1)
            if self.verbose: 
                print(f'Combined shape: {combined.shape}')
        
            x = F.relu(self.fc1(combined))
            if self.verbose: 
                print(f'x shape of combined: {x.shape}')
            
            x = self.dropout(x)
            
            x = F.relu(self.fc2(x))
            x = self.dropout(x)
            
            final_output = self.fc3(x)
        
        if not self.regression:
            final_output = F.sigmoid(final_output)
            
        return final_output