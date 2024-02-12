import pandas as pd
import numpy as np
from gerrychain.random import random
import lightgbm as lgb
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.model_selection import train_test_split
import optuna
from optuna.samplers import TPESampler
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

class CrashTest:
    def __init__(self, base_name, kmer_window, crash_test, only_breaks, cols_to_remove=None, default=True):
        self.base_name = base_name
        self.kmer_window = kmer_window
        self.crash_test = crash_test
        self.only_breaks = only_breaks
        self.full_name = (
            self.base_name + '_' + 
            'kmerwindow-' + str(self.kmer_window) +
            'crash_test-' + str(self.crash_test).upper() + 
            'only_breaks-' + str(self.only_breaks).upper() +
            'kmer-8_seed-41964_zscore'
        )
        self.cols_to_remove = cols_to_remove
        self.default = default
        self.cpu = 2
        self.seed = 1234
        self.threshold = 0.5
        self.n_trials = 10
        self.n_train_iter = self.n_trials * 4
        self.df_colnames = pd.read_csv(f'../data/ML_demo/{self.full_name}_colnames.csv')
        self.columns_breaks = ['BREAKS', 'LONGRANGE_PREPARAMS']
        self.columns_kmers = ['KMER_COUNTS']
        self.exp_labels = ['Breakage', 'Kmer_counts']
        
        if not os.path.exists('../figures/ML_demo'):
            os.makedirs('../figures/ML_demo', exist_ok=True)
            
        if not os.path.exists('../data/ML_demo'):
            os.makedirs('../data/ML_demo', exist_ok=True)
        
        np.random.seed(self.seed)
        random.seed(self.seed)
        
    def get_data(self, filename):        
        df = pd.read_parquet(
            '../data/ML_demo/' + 
            self.full_name + '_' + filename + '.parquet', 
            use_threads=True
        )
        target_col = 'predictor'
        features = [col for col in df.columns if col != target_col]
        
        if self.cols_to_remove is not None:
            features = list(filter(lambda x: self.cols_to_remove not in x, features))
        
        y = df[target_col].map({'YES': 1, 'NO': 0}).to_numpy(dtype=np.float32)
        X = df[features]
        return X, y
    
    # optuna objective function for LightGBM
    def objective_lgb(self, trial, X_train, X_test, y_train, y_test):
        try:
            params = {
                "objective": "binary",
                "metric": "binary_logloss",
                "verbosity": -1,
                "boosting_type": "gbdt",
                "device_type": "gpu",
                "num_threads": self.cpu,
                "reg_alpha": trial.suggest_float("reg_alpha", 1e-8, 10.0, log=True),
                "reg_lambda": trial.suggest_float("reg_lambda", 1e-8, 10.0, log=True),
                "num_leaves": trial.suggest_int("num_leaves", 2, 100),
                "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.1),
                "n_estimators": trial.suggest_int("n_estimators", 1000, 10000),
                "colsample_bytree": trial.suggest_float("colsample_bytree", 0.1, 1.0),
                "subsample": trial.suggest_float("subsample", 0.1, 1.0)
            }

            model = lgb.LGBMClassifier(**params, random_state=self.seed)
            model.fit(X_train, y_train)
            preds = model.predict_proba(X_test)[:, 1]
            fpr, tpr, _ = roc_curve(y_test, preds)
            auroc = auc(fpr, tpr)
            return auroc

        except lgb.basic.LightGBMError as e:
            return -np.inf
        
    def objective_logreg(self, trial, X_train, X_test, y_train, y_test):
        C_float = trial.suggest_loguniform("C", 1e-2, 1)

        model = LogisticRegression(C=C_float, random_state=self.seed)
        model.fit(X_train, y_train)
        y_probs = model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_probs)
        auroc = auc(fpr, tpr)
        return auroc

    def train_and_evaluate_models(self):    
        X_train, y_train = self.get_data(filename='train')
        X_test, y_test = self.get_data(filename='test')

        downsampling = [1000, 3000, 5000, 10000, 100000, 500000]
        full_train_size = len(X_train)
        downsampling.append(full_train_size)

        scaler = StandardScaler()
        results = []
        
        # downsample the training data
        for size in downsampling:
            if size != full_train_size:            
                X_train_downsampled, _, y_train_downsampled, _ = train_test_split(
                    X_train, y_train, stratify=y_train, train_size=size, random_state=self.seed
                )
                X_train_downsampled = X_train_downsampled.reset_index()
            else:
                X_train_downsampled = X_train
                y_train_downsampled = y_train

            # train and test for each subset of columns
            for i, subset_columns in enumerate([self.columns_breaks, self.columns_kmers]):
                print(f"Fitting models to columns: {self.exp_labels[i]} with {size} samples.")
                
                cols_of_int = self.df_colnames[self.df_colnames['group'].isin(subset_columns)]['columns'].to_list()
                X_train_np = X_train_downsampled[X_train_downsampled.columns.intersection(cols_of_int)].to_numpy(dtype=np.float32)
                X_test_np = X_test[X_test.columns.intersection(cols_of_int)].to_numpy(dtype=np.float32)

                # train the LightGBM model with the best parameters found
                if self.default:
                    best_params = dict()
                    best_params["objective"]="binary"
                    best_params["metric"]="binary_logloss"
                    best_params["boosting_type"]="gbdt"
                    best_params["device_type"]="gpu"
                    best_params["num_threads"]=self.cpu
                    best_params['verbose']=-1
                    lgbc = lgb.LGBMClassifier(**best_params, random_state=self.seed)
                else:
                    sampler_lgb = TPESampler(seed=self.seed)
                    study_lgb = optuna.create_study(
                        direction="maximize",
                        sampler=sampler_lgb,
                        pruner=optuna.pruners.HyperbandPruner(
                            min_resource=1, 
                            max_resource=self.n_train_iter, 
                            reduction_factor=3
                        )
                    )
                    study_lgb.optimize(lambda trial: self.objective_lgb(
                        trial=trial, 
                        X_train=X_train, X_test=X_test_np,
                        y_train=y_train_downsampled, y_test=y_test
                    ), n_trials=self.n_trials)
                    best_params = study_lgb.best_params
                    best_params['verbose'] = -1
                    lgbc = lgb.LGBMClassifier(**best_params, random_state=self.seed)

                lgbc.fit(X_train_np, y_train_downsampled)
                lgbc_y_pred = lgbc.predict_proba(X_test_np)[:, 1]
                lgbc_fpr, lgbc_tpr, _ = roc_curve(y_test, lgbc_y_pred)
                lgbc_auroc = auc(lgbc_fpr, lgbc_tpr)
                lgbc_acc = accuracy_score(y_test, np.where(lgbc_y_pred > self.threshold, 1, 0))
                
                # train the linear model with the best parameters found
                X_train_np = scaler.fit_transform(X_train_np)
                X_test_np = scaler.transform(X_test_np)
                
                if self.default:
                    lmc = LogisticRegression(random_state=self.seed)
                else:
                    sampler_logreg = TPESampler(seed=self.seed)
                    study_logreg = optuna.create_study(
                        direction="maximize", 
                        sampler=sampler_logreg,
                        pruner=optuna.pruners.HyperbandPruner(
                            min_resource=1, 
                            max_resource=self.n_train_iter, 
                            reduction_factor=3
                        )
                    )
                    study_logreg.optimize(lambda trial: self.objective_logreg(
                        trial=trial, 
                        X_train=X_train, X_test=X_test_np,
                        y_train=y_train_downsampled, y_test=y_test
                    ), n_trials=self.n_trials)
                    best_logreg_params = study_logreg.best_params
                    lmc = LogisticRegression(**best_logreg_params, random_state=self.seed)

                lmc.fit(X_train_np, y_train_downsampled)
                lmc_y_probs = lmc.predict_proba(X_test_np)[:, 1]
                lmc_fpr, lmc_tpr, _ = roc_curve(y_test, lmc_y_probs)
                lmc_auroc = auc(lmc_fpr, lmc_tpr)
                lmc_acc = accuracy_score(y_test, np.where(lmc_y_probs > self.threshold, 1, 0))

                if self.default:
                    results.append({
                        'downsampling_size': size,
                        'feature_group': self.exp_labels[i],
                        'lgb_accuracy': lgbc_acc,
                        'lgb_auroc': lgbc_auroc,
                        'lmc_accuracy': lmc_acc,
                        'lmc_auroc': lmc_auroc
                    })
                else:
                    results.append({
                        'downsampling_size': size,
                        'feature_group': self.exp_labels[i],
                        'lgb_accuracy': lgbc_acc,
                        'lgb_auroc': lgbc_auroc,
                        'lmc_accuracy': lmc_acc,
                        'lmc_auroc': lmc_auroc,
                        'lgb_best_params': [list(best_params.items())],
                        'lmc_best_params': [list(best_logreg_params.items())]
                    })

        return pd.DataFrame(results)
    def run_models(self):
        df_res = self.train_and_evaluate_models()
        self.save_results(df_res)

    def save_results(self, df_res):
        df_res.to_csv(
            '../data/ML_demo/' + self.full_name + '_' +
            'LM_vs_LightGBM_default-' + 
            str(self.default).upper() + '.csv',
            index=False
        )
        
        self.df_res_long = pd.melt(
            df_res,
            id_vars=['downsampling_size', 'feature_group'],
            var_name=['metric_model']
        ).assign(
            model=lambda df: df['metric_model'].str.split('_').str[0].replace({
                'lmc': 'Linear', 'lgb': 'LightGBM'
            }),
            metric=lambda df: df['metric_model'].str.split('_').str[1]
        ).drop(
            columns=['metric_model']
        )

    def plot_results_auc(self, ymin, ymax):
        df_res_long_auroc = self.df_res_long[self.df_res_long['metric'] == 'auroc']

        ax = sns.FacetGrid(
            df_res_long_auroc, 
            col="model", 
            col_order=['Linear', 'LightGBM'],
            hue="feature_group", 
            height=5, 
            aspect=1, 
            palette=['#2166ac', '#b2182b']
        )
        ax = ax.map(sns.scatterplot, "downsampling_size", "value", alpha=0.75)
        ax = ax.map(sns.pointplot, "downsampling_size", "value", alpha=0.75, ci=None)

        ax.set_axis_labels('Number of Samples', 'AUROC score')
        ax.set_titles('{col_name} Model')
        ax.add_legend(
            title='Feature Group', 
            loc='upper center', 
            bbox_to_anchor=(0.5, -0.1), 
            ncol=10
        )
        ax.set(ylim=(ymin, ymax))
        
        plt.gcf().axes[0].tick_params(labelsize=12)
        plt.tight_layout()
        plt.savefig(
            '../figures/ML_demo/' + self.full_name + '_' +
            'LM_vs_LightGBM_default_auroc_default-' + 
            str(self.default).upper() + '.pdf',
            dpi=600, bbox_inches='tight'
        )
        
    def plot_results_acc(self, ymin, ymax):
        df_res_long_acc = self.df_res_long[self.df_res_long['metric'] == 'accuracy']

        ax = sns.FacetGrid(
            df_res_long_acc, 
            col="model", 
            col_order=['Linear', 'LightGBM'],
            hue="feature_group", 
            height=5, 
            aspect=1, 
            palette=['#2166ac', '#b2182b']
        )
        ax = ax.map(sns.scatterplot, "downsampling_size", "value", alpha=0.75)
        ax = ax.map(sns.pointplot, "downsampling_size", "value", alpha=0.75, ci=None)

        ax.set_axis_labels('Number of Samples', 'Accuracy score')
        ax.set_titles('{col_name} Model')
        ax.add_legend(
            title='Feature Group', 
            loc='upper center', 
            bbox_to_anchor=(0.5, -0.1), 
            ncol=10
        )
        ax.set(ylim=(ymin, ymax))

        plt.gcf().axes[0].tick_params(labelsize=12)
        plt.tight_layout()
        plt.savefig(
            '../figures/ML_demo/' + self.full_name + '_' +
            'LM_vs_LightGBM_default_acc_default-' + 
            str(self.default).upper() + '.pdf',
            dpi=600, bbox_inches='tight'
        )