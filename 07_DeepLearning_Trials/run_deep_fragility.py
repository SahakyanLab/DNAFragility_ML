import argparse
import os
import sys
import pandas as pd
import numpy as np

from tqdm import tqdm 

from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, gaussian_kde
from sklearn.metrics import roc_auc_score, confusion_matrix, ConfusionMatrixDisplay, accuracy_score, precision_recall_curve, auc

import torch
from torch.utils.data import DataLoader
import torch.optim as optim
import equirc as eq

import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from dataset_subset import DatasetSubset
from cnn_features_model import CNNandFeaturesModel
from utils import seed_everything
from losses import RMSELoss

parser = argparse.ArgumentParser(description="Trial deep fragility models.")
parser.add_argument(
    "--regression", 
    action=argparse.BooleanOptionalAction, 
    default=False,
    help="If provided, runs regression, else classification task"
)
parser.add_argument(
    "--with_featmatrix", 
    action=argparse.BooleanOptionalAction, 
    default=False,
    help="If provided, with feature matrix, else only one-hot-encoding"
)
parser.add_argument(
    "--with_LSTM", 
    action=argparse.BooleanOptionalAction, 
    default=False,
    help="If provided, model architecture with LSTM, else just CNN"
)
parser.add_argument(
    "-epochs", "--epochs", required=True, 
    type=int, help="Number of epochs"
)
parser.add_argument(
    "-bw", "--bin_width", required=True, 
    type=int, help="Size of each bin width"
)
args = parser.parse_args()

# path to folders
regression = args.regression
with_additional_feat = args.with_featmatrix
with_LSTM = args.with_LSTM
epochs = args.epochs
bw = args.bin_width
verbose = False

# path to folders
train_ranges_path = f"./data/train_regression-{regression}.csv"
test_ranges_path = f"./data/test_regression-{regression}.csv"
sequence = "../data/ref/chm13v2.0.fa"

# model with engineered features?
train_features_path = f"./data/train_regression-{regression}.parquet" if with_additional_feat else None
test_features_path = f"./data/test_regression-{regression}.parquet" if with_additional_feat else None

# modelling parameters
batch_size = 64 if regression else 256
num_workers = 0
learning_rate = 6e-6
dropout = 0
weight_decay = 1e-4
threshold = 0.5
seed = 1234
seed_everything(seed=seed)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

train_data = DatasetSubset(
    ranges=train_ranges_path, 
    sequence=sequence, 
    features_path=train_features_path,
    device=device
)
test_data = DatasetSubset(
    ranges=test_ranges_path, 
    sequence=sequence, 
    features_path=test_features_path,
    device=device
)

if with_additional_feat:
    scaler = StandardScaler()
    
    # transform train data
    train_features = train_data.features.loc[:, ].values
    train_features_scaled = scaler.fit_transform(train_features)
    train_data.features.loc[:, ] = train_features_scaled
    
    # apply transformation to test data
    test_features = test_data.features.loc[:, ].values
    test_features_scaled = scaler.transform(test_features)
    test_data.features.loc[:, ] = test_features_scaled

train_dataloader = DataLoader(
    train_data, 
    batch_size=batch_size, 
    num_workers=num_workers
)
test_dataloader = DataLoader(
    test_data, 
    batch_size=batch_size, 
    num_workers=num_workers
)

# if GPU is available
model = CNNandFeaturesModel(
    data=train_data, 
    with_LSTM=with_LSTM,
    dropout=dropout,
    regression=regression,
    verbose=verbose
).to(device)

optim = torch.optim.AdamW(
    model.parameters(),
    lr=learning_rate,
    weight_decay=weight_decay
)
loss_fn = RMSELoss() if regression else torch.nn.BCELoss()

save_loss, save_corr_coef, save_accuracy = [], [], []
for epoch in range(1, epochs+1):
    model.train()

    total_loss, correct = 0, 0  
    outputs_list, targets_list = [], []

    train_pb = tqdm(train_dataloader, total=len(train_dataloader))

    for batch in train_pb:
        # zero gradients for every batch
        optim.zero_grad()
        
        if with_additional_feat:
            inputs, targets, feature_matrix = batch
            outputs = model(inputs, feature_matrix)
        else:
            inputs, targets = batch            
            outputs = model(inputs, None)

        # loss and gradients
        loss = loss_fn(outputs.squeeze(), targets)
        loss.backward()
        
        # adjust learning weights
        optim.step()

        # update loss
        total_loss += loss.item()
        
        # accumulate outputs and targets
        outputs_list.extend(outputs.tolist())
        targets_list.extend(targets.tolist())
        
        if not regression:
            outputs_tensor_conv = torch.where(outputs.squeeze() > threshold, 1, 0)
            accuracy = float((outputs_tensor_conv == targets).float().mean()) * 100

        train_pb.set_description(
            "Epoch " + str(epoch) + "/" + str(epochs) + ", " + 
            "Loss: " + format(loss.item(), '.3f') + 
            (", Accuracy: " + format(accuracy, '.3f') if not regression else "")
        )
        
    # average results per epoch
    avg_loss = total_loss / len(train_dataloader)
    outputs_tensor = torch.tensor(outputs_list)
    targets_tensor = torch.tensor(targets_list)
    
    if regression:
        outputs_tensor = (2 ** outputs_tensor.view(-1, 1)) - 1
        targets_tensor = (2 ** targets_tensor.view(-1, 1)) - 1
    
        stack_tensor = torch.cat((outputs_tensor, targets_tensor), dim=1)
        corr_coef_matrix = torch.corrcoef(stack_tensor.t())
        corr_coef = corr_coef_matrix[0, 1].item()
    else:
        outputs_tensor_conv = torch.where(outputs_tensor > threshold, 1, 0).reshape(-1)
        accuracy = float((outputs_tensor_conv == targets_tensor).float().mean()) * 100

    print(
        "Epoch: " + str(epoch) + "/" + str(epochs) + ", " + 
        "Avg Loss: " + format(avg_loss, '.3f') + ", " + 
        (
            "Pearson R: " + format(corr_coef, '.3f') if regression else "Accuracy: " + 
            format(accuracy, '.3f')
        )
    )
    
    save_loss.append(avg_loss)
    if regression:
        save_corr_coef.append(corr_coef)
    else:
        save_accuracy.append(accuracy)
    
save_loss = np.array(save_loss)
if regression:
    save_corr_coef = np.array(save_corr_coef)
    training_metrics = pd.DataFrame({
        'Epoch': range(1,epochs+1),
        'Training_Loss': save_loss,
        'Training_PearsonR': save_corr_coef
    })
else:
    save_accuracy = np.array(save_accuracy)
    training_metrics = pd.DataFrame({
        'Epoch': range(1,epochs+1),
        'Training_Loss': save_loss,
        'Training_Accuracy': save_accuracy
    })
    
if not os.path.exists('../figures/DeepLearning_Trials/'):
    os.makedirs('../figures/DeepLearning_Trials/', exist_ok=True)
        
training_metrics.to_csv(
    './data/TrainingMetrics_' + 
    'Regression_' + str(regression) + '_' + 
    'FeatMatrix_' + str(with_additional_feat) + '_' + 
    'LSTM_' + str(with_LSTM) + '_' + 
    'BinWidth_' + str(bw) + '.csv', 
    index=False
)
        
all_outputs, all_targets = [], []
with torch.no_grad():
    model.eval()
    correct, total, total_loss = 0, 0, 0
    outputs_list, targets_list = [], []
    
    test_pb = tqdm(test_dataloader, total=len(test_dataloader))
    
    for batch in test_pb:
        if with_additional_feat:
            inputs, targets, feature_matrix = batch
            outputs = model(inputs, feature_matrix)
        else:
            inputs, targets = batch            
            outputs = model(inputs, None)

        # loss and gradients
        loss = loss_fn(outputs.squeeze(), targets)
        
        # update loss
        total_loss += loss.item()
        
        # accumulate outputs and targets
        outputs_list.extend(outputs.tolist())
        targets_list.extend(targets.tolist())
        
        if not regression:
            outputs_tensor_conv = torch.where(outputs.squeeze() > threshold, 1, 0)
            accuracy = float((outputs_tensor_conv == targets).float().mean()) * 100

        test_pb.set_description(
            "Loss: " + format(loss.item(), '.3f') + 
            (", Accuracy: " + format(accuracy, '.3f') if not regression else "")
        )
        
    # average results per epoch
    avg_loss = total_loss / len(test_dataloader)
    outputs_tensor = torch.tensor(outputs_list)
    targets_tensor = torch.tensor(targets_list)
    
    if regression:
        outputs_tensor = outputs_tensor.view(-1, 1)
        targets_tensor = targets_tensor.view(-1, 1)
        stack_tensor = torch.cat((outputs_tensor, targets_tensor), dim=1)
        corr_coef_matrix = torch.corrcoef(stack_tensor.t())
        corr_coef_logscale = corr_coef_matrix[0, 1].item()
        
        outputs_tensor = (2 ** outputs_tensor) - 1
        targets_tensor = (2 ** targets_tensor) - 1
        stack_tensor = torch.cat((outputs_tensor, targets_tensor), dim=1)
        corr_coef_matrix = torch.corrcoef(stack_tensor.t())
        corr_coef_truescale = corr_coef_matrix[0, 1].item()
    else:
        outputs_tensor_conv = torch.where(outputs_tensor > threshold, 1, 0).reshape(-1)
        accuracy = float((outputs_tensor_conv == targets_tensor).float().mean()) * 100
        auroc = roc_auc_score(outputs_tensor_conv, targets_tensor)
        
        precision, recall, _ = precision_recall_curve(targets_tensor, outputs_tensor_conv)
        auprc = auc(recall, precision)

    print(
        "Avg Loss: " + format(avg_loss, '.3f') + 
        (
            ", Pearson R: " + format(corr_coef_truescale, '.3f') if regression else 
            ", Accuracy: " + format(accuracy, '.3f') + 
            ", AUROC: " + format(auroc, '.3f') + 
            ", AUPRC: " + format(auprc, '.3f')
        )
    )
    
all_outputs = np.array(outputs_list).reshape(-1)
all_targets = np.array(targets_list)

# Convert the lists or numpy arrays to a DataFrame
model_output = pd.DataFrame({
    'Outputs': all_outputs,
    'Targets': all_targets
})

model_output.to_csv(
    './data/ModelOutput_' + 
    'Regression_' + str(regression) + '_' + 
    'FeatMatrix_' + str(with_additional_feat) + '_' + 
    'LSTM_' + str(with_LSTM) + '_' + 
    'BinWidth_' + str(bw) + '.csv', 
    index=False
)

def plot_confusion_matrix(y_true, y_pred, auroc=None, auprc=None):
    fontsize = 12
    
    # Evaluate the model
    accuracy = accuracy_score(y_true, y_pred)
    conf_matrix = confusion_matrix(y_true, y_pred) / len(y_pred)

    print(f'Accuracy: {accuracy:.2f}')

    # Plot the confusion matrix
    tick_labels = ['False', 'True']
    
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        conf_matrix, annot=True, fmt=".2%", 
        annot_kws={"size": fontsize},
        xticklabels=tick_labels, yticklabels=tick_labels
    )
    
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=fontsize)
    
    plt.xlabel('Predicted', fontsize=fontsize)
    plt.ylabel('True', fontsize=fontsize)
    
    title = (
        "Accuracy: " + format(accuracy, '.3f') + 
        ", AUROC: " + format(auroc, '.3f') + 
        ", AUPRC: " + format(auprc, '.3f')
    )
    plt.title(title, fontsize=fontsize)
    plt.savefig(
        '../figures/DeepLearning_Trials/ConfusionMatrix_' + 
        'Regression_' + str(regression) + '_' + 
        'FeatMatrix_' + str(with_additional_feat) + '_' + 
        'LSTM_' + str(with_LSTM) + '_' + 
        'BinWidth_' + str(bw) + '.pdf', 
        dpi=300
    )

def plot_predictions_vs_true(
        yhat, y, bw, corr_coef=[corr_coef_logscale, corr_coef_truescale], 
        log_space=True, y_range=None, fix_yrange=False
    ):
    labelsize = 12
    fontsize = labelsize + 2
    
    if not log_space:
        y = 2**y - 1
        yhat = 2**yhat - 1
    
    # Estimate density
    xy = np.vstack([y, yhat])
    z = gaussian_kde(xy)(xy)
    
    # new figure with side-by-side plots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # Create the first subplot on the left side
    ax_1 = axs[0]
    sns.scatterplot(
        x=y, y=yhat, 
        hue=z, 
        size=z, 
        alpha=0.2, 
        c=z,
        cmap="jet",
        legend=False, 
        ax=ax_1
    )
    sns.lineplot(
        x=[-0.05 * max(y), max(y)], 
        y=[-0.05 * max(y), max(y)], 
        color='grey',
        ax=ax_1
    )
    # line regression plot
    sns.regplot(
        x=y, y=yhat, 
        scatter=False, 
        color='darkred',
        ax=ax_1
    )
    # display pearson correlation coefficient
    if fix_yrange:
        ax_1.set_xlim([-0.05 * max(y), max(y)])
        ax_1.set_ylim([-0.05 * max(y), max(y)])
    
    ax_1.set_xlabel('True (log)' if log_space else 'True', fontsize=fontsize)
    ax_1.set_ylabel('Predicted (log)' if log_space else 'Predicted', fontsize=fontsize)

    # Create the second subplot on the right side
    ax_2 = axs[1]
    sns.kdeplot(
        y, label='True', alpha=0.8, 
        lw=3, fill=False, ax=ax_2, 
        bw_adjust=1, color='#2166ac'
    )
    sns.kdeplot(
        yhat, label='Prediction', alpha=0.8, 
        lw=3, fill=False, ax=ax_2, 
        bw_adjust=1, color='#b2182b'
    )
    ax_2.set_xlabel('Number of breaks (log)' if log_space else 'Number of breaks', fontsize=fontsize)
    ax_2.set_ylabel('Density', fontsize=fontsize)
    
    for i in range(len(axs)):
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        
        # increase tick label sizes
        axs[i].tick_params(axis='both', which='major', labelsize=labelsize)

    title = "Pearson R: " + (format(corr_coef[0], '.3f') if log_space else format(corr_coef[1], '.3f'))
    ax_1.set_title(title, fontsize=labelsize)
    ax_2.legend(fontsize=fontsize)
    
    plt.tight_layout()
    fig.savefig(
        '../figures/DeepLearning_Trials/Regression_' + 
        'FeatMatrix_' + str(with_additional_feat) + '_' + 
        'LSTM_' + str(with_LSTM) + '_' + 
        'BinWidth_' + str(bw) + '_' + 
        'LogScale_' + str(log_space) + '.pdf', 
        dpi=300
    )
    
if regression:
    plot_predictions_vs_true(
        y=all_outputs, yhat=all_targets, 
        bw=bw, fix_yrange=True, log_space=True
    )

    plot_predictions_vs_true(
        y=all_outputs, yhat=all_targets, 
        bw=bw, fix_yrange=True, log_space=False
    )
else:
    plot_confusion_matrix(
        y_true=targets_tensor.numpy(), 
        y_pred=outputs_tensor_conv.numpy(),
        auroc=auroc,
        auprc=auprc
    )