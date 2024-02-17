import argparse
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_dir, '../lib'))

from crashtest import CrashTest

parser = argparse.ArgumentParser(description="Run ML crash tests.")
parser.add_argument(
    "-name", "--base_name", required=True, 
    type=str, help="Long range kmer size"
)
parser.add_argument(
    "-k", "--kmer_window", required=True, 
    type=int, help="Long range kmer size"
)
args = parser.parse_args()

base_name = args.base_name
only_breaks = True
crash_test = True
cols_to_remove = 'K562_cells_Top2'

model = CrashTest(
    base_name=base_name, 
    kmer_window=args.kmer_window, 
    crash_test=crash_test, 
    only_breaks=only_breaks, 
    cols_to_remove=cols_to_remove,
    default=True
)
model.run_models()

if args.kmer_window == 3:
    if base_name == "K562_DMSO_DSBs":
        ymin_auc, ymax_auc = 0.74, 0.83
        ymin_acc, ymax_acc = 0.71, 0.76
    else:
        ymin_auc, ymax_auc = 0.70, 0.82
        ymin_acc, ymax_acc = 0.65, 0.75
else:
    ymin_auc, ymax_auc = 0.65, 0.83
    if base_name == "K562_DMSO_DSBs":
        ymin_acc, ymax_acc = 0.58, 0.78
    else:
        ymin_acc, ymax_acc = 0.58, 0.78
    
model.plot_results_auc(ymin=ymin_auc, ymax=ymax_auc)
model.plot_results_acc(ymin=ymin_acc, ymax=ymax_acc)