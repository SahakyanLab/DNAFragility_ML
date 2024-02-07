import pandas as pd
import numpy as np
import pickle
from gerrychain.random import random
import lightgbm as lgb
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, accuracy_score
import optuna
from optuna.samplers import TPESampler
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os
import sys
import warnings
warnings.filterwarnings('ignore')

class OptunaModel:
    def __init__(
        self, base_name, base_dir, 
        kmer_window, crash_test, only_breaks, model_type, 
        n_trials=50, final_model=False, cols_to_remove=None
    ):
        self.base_name = base_name
        self.base_dir = base_dir
        self.kmer_window = kmer_window
        self.crash_test = crash_test
        self.only_breaks = only_breaks
        self.final_model = final_model
        
        if self.final_model:
            self.full_name = 'best_LGBM_model'
        else:
            self.full_name = (
                self.base_name + '_' + 
                'kmerwindow-' + str(self.kmer_window) +
                'crash_test-' + str(self.crash_test).upper() + 
                'only_breaks-' + str(self.only_breaks).upper() +
                'kmer-8_seed-41964_zscore'
            )
        self.cols_to_remove = cols_to_remove
        self.n_trials = n_trials
        self.cpu = 2
        self.seed = 1234
        self.threshold = 0.5
        self.n_train_iter = self.n_trials * 4
        self.model_type = model_type
        
        self.group_colors = {
            'QM_PARAMETERS': 'salmon',
            'G4MAP': 'blue',
            'BREAKS': 'purple',
            'LONGRANGE_PREPARAMS': 'purple',
            'G4_REGEX': 'darkorange',
            'SINGLETON': 'brown',
            'GC_CONTENT': 'pink',
            'GC_COUNT': 'gray',
            'KMER_COUNTS': 'olive',
            'VIENNA_RNA': 'cornflowerblue',
            'DNA_SHAPE': 'deeppink'
        }
        if self.final_model:
            self.df_colnames = pd.read_csv(
                "../data/ML_demo/K562_Top2_mediated_DSBs_" + 
                "kmerwindow-5crash_test-TRUE" + 
                "only_breaks-TRUEkmer-8_" + 
                "seed-41964_zscore_colnames.csv"
            )
        else:
            self.df_colnames = pd.read_csv(f'../data/{self.base_dir}/{self.full_name}_colnames.csv')
        self.df_colnames['hex'] = self.df_colnames['group'].map(self.group_colors)
        
        self.columns_breaks = ['BREAKS', 'LONGRANGE_PREPARAMS']
        self.columns_kmers = ['KMER_COUNTS']
        self.exp_labels = ['Breakage', 'Kmer_counts']
        
        if not os.path.exists(f'../figures/{self.base_dir}'):
            os.makedirs(f'../figures/{self.base_dir}', exist_ok=True)
            
        if not os.path.exists(f'../data/{self.base_dir}'):
            os.makedirs(f'../data/{self.base_dir}', exist_ok=True)
            
        np.random.seed(self.seed)
        random.seed(self.seed)
        
    def get_data(self, filename, to_numpy=True, verbose=True):   
        if verbose:
            t1 = time.time()
            cur_msg = "Importing " + filename + " dataset"
            l = ''.join(['.' for _ in range(60 - len(cur_msg))])
            sys.stdout.write(f"{cur_msg}{l}")
            sys.stdout.flush()
             
        path_to_file_name = '../data/' + self.base_dir + '/'
        if not self.final_model:
            path_to_file_name = path_to_file_name + self.full_name + '_'
            
        df = pd.read_parquet(
            path_to_file_name + filename + '.parquet', 
            use_threads=True
        )
        target_col = 'predictor'
        features = [col for col in df.columns if col != target_col]
        
        if self.cols_to_remove is not None:
            features = list(filter(
                lambda x: self.filter_features(x, self.cols_to_remove), features
            ))
        
        y = df[target_col].map({'YES': 1, 'NO': 0})
        X = df[features]
        
        if to_numpy:
            y = y.to_numpy(dtype=np.float32)
            X = X.to_numpy(dtype=np.float32)
            
        if verbose:
            total_time = time.time() - t1
            sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
            sys.stdout.flush()
        
        return X, y, features
    
    def filter_features(self, feature, cols_to_remove):
        feature_lower = feature.lower()
        cols_to_remove_lower = cols_to_remove.lower()
        prefix = cols_to_remove_lower.split('_')[0]

        if feature_lower == cols_to_remove_lower:
            return False
        if any(seq.lower() in feature_lower for seq in ['ATACseq', 'Chipseq', 'Dnaseseq', 'FAIREseq', 'Epigenome']):
            return prefix in feature_lower
        return True
    
    # optuna objective function for LightGBM
    def objective_lgb(self, trial, X_train, X_test, y_train, y_test):
        try:
            params = {
                "objective": "binary",
                "metric": "binary_logloss",
                "verbosity": -1,
                "boosting_type": "gbdt",
                "device": "gpu",
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
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        
        params = {
            # Regularisation strength
            "C": trial.suggest_loguniform("C", 1e-2, 1),
            
            # Tolerance for stopping criteria
            "tol": trial.suggest_uniform("tol", 1e-6, 1e-3),
            
            # Maximum number of iterations for solver to converge
            "max_iter": trial.suggest_int("max_iter", 100, 1000)
        }
        
        try:
            model = LogisticRegression(**params, penalty='l2', random_state=self.seed)
            model.fit(X_train_scaled, y_train)
            y_probs = model.predict_proba(X_test_scaled)[:, 1]
            
            fpr, tpr, _ = roc_curve(y_test, y_probs)
            auroc = auc(fpr, tpr)
            
            return auroc
        except Exception as e:
            return -np.inf
        
    def run_optuna_trials(self):
        X_train, y_train, _ = self.get_data(filename='train', to_numpy=True)
        X_test, y_test, _ = self.get_data(filename='test', to_numpy=True)
        
        t1 = time.time()
        cur_msg = "Starting Optuna hyperparameter tuning"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}\n")
        sys.stdout.flush()
                
        # Make the sampler behave in a deterministic way.
        sampler = TPESampler(seed=self.seed)
        self.study = optuna.create_study(
            direction="maximize",
            sampler=sampler,
            pruner=optuna.pruners.HyperbandPruner(
                min_resource=1, 
                max_resource=self.n_train_iter, 
                reduction_factor=3
            )
        )

        if self.model_type == 'LGBM':
            self.study.optimize(lambda trial: self.objective_lgb(
                trial=trial, 
                X_train=X_train, X_test=X_test, 
                y_train=y_train, y_test=y_test
            ), n_trials=self.n_trials)
        elif self.model_type == 'LM':
            self.study.optimize(lambda trial: self.objective_logreg(
                trial=trial, 
                X_train=X_train, X_test=X_test, 
                y_train=y_train, y_test=y_test
            ), n_trials=self.n_trials)
                
        df_studytrials = self.study.trials_dataframe()
        df_studytrials.rename(columns={'value': 'AUROC'}, inplace=True)
        df_studytrials.to_csv(
            '../data/' + self.base_dir + '/' + 
            self.full_name + '_' +
            (self.model_type + '_model_' if not self.final_model else '') + 
            'studytrials.csv',
            index=False
        )
        
        print(f"Number of finished trials: {len(self.study.trials)}")
        print("Best trial:")
        print(f"  AUROC: {self.study.best_trial.value:.3f}")
        print("  Params: ")
        for key, value in self.study.best_trial.params.items():
            print(f"    {key}: {value:.3f}")
            
        # importance of each hyper-parameter
        hyp_imp = optuna.importance.get_param_importances(self.study)
        hyp_imp = pd.DataFrame(hyp_imp.items(), columns=['feature', 'importance'])
        hyp_imp.sort_values('importance', ascending=False, inplace=True)

        fig, ax = plt.subplots(figsize=(6,4))
        y_pos = range(len(hyp_imp))
        ax.hlines(y=y_pos, xmin=0, xmax=hyp_imp['importance'], color='grey')
        ax.scatter(hyp_imp['importance'], y_pos, color='darkorange', s=100)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(hyp_imp['feature'])

        # Highest importance at the top
        ax.invert_yaxis()

        ax.set_xlabel('Hyperparameter Importance')
        ax.set_ylabel('Hyperparameter')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        fig.savefig(
            '../figures/' + self.base_dir + '/' + 
            self.full_name + '_' +
            ('best_' + self.model_type + '_' if not self.final_model else '') + 
            'HyperparameterImportance.pdf',
            dpi=600, bbox_inches='tight'
        )
        
        if self.final_model:
            all_auroc_values, all_trees = [], []
            for trial in self.study.trials:
                    if trial.state == optuna.trial.TrialState.COMPLETE:
                        all_auroc_values.append(trial.values[0])
                        all_trees.append(trial.params['n_estimators'])
                        
            model_complexity = pd.DataFrame({
                'Number of trees': all_trees, 
                'AUROC': all_auroc_values
            })

            # visualise aucroc vs. complexity
            plt.figure(figsize=(7,5))
            sns.scatterplot(
                data=model_complexity, 
                x='Number of trees', y='AUROC', 
                color='navy', s=50
            )

            plt.title('Performance vs. complexity of model')
            plt.tight_layout()
            plt.savefig(
                '../figures/' + self.base_dir + '/' + 
                self.full_name + '_Classifier_AUROC-vs-Trees.pdf',
                dpi=600, bbox_inches='tight'
            )
            
    def run_optim_model(self):
        X_train, y_train, _ = self.get_data(
            filename='train', to_numpy=True, verbose = False
        )
        X_test, y_test, self.feature_names = self.get_data(
            filename='test', to_numpy = False, verbose = False
        )
            
        t1 = time.time()
        cur_msg = "Running optimised model"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        if self.model_type == 'LGBM': 
            self.study.best_trial.params["objective"]="binary"
            self.study.best_trial.params["metric"]="binary_logloss"
            self.study.best_trial.params["verbosity"]=1
            self.study.best_trial.params["boosting_type"]="gbdt"
            self.study.best_trial.params["device"]="gpu"
            self.study.best_trial.params["num_threads"]=self.cpu
            
            self.best_model = lgb.LGBMClassifier(
                **self.study.best_trial.params, 
                random_state=self.seed
            )
        elif self.model_type == 'LM':
            self.best_model = LogisticRegression(
                **self.study.best_trial.params, 
                random_state=self.seed
            )
        self.best_model.fit(X_train, y_train)
                
        # save model
        base_name = '../data/' + self.base_dir + '/' + self.full_name + '_' 
        if self.model_type == 'LGBM':
            self.best_model.booster_.save_model(
                base_name + 'best_LGBM_model.txt'
            )
        elif self.model_type == 'LM':
            file_name = base_name + 'best_LM_model.pkl'
            with open(file_name, 'wb') as file:
                pickle.dump(self.best_model, file)
            
        # Make predictions on the test set
        y_probs = self.best_model.predict_proba(X_test)[:, 1]
        self.fpr, self.tpr, thresholds = roc_curve(y_test, y_probs)
        auroc = auc(self.fpr, self.tpr)
        acc = accuracy_score(np.where(y_probs > 0.5, 1, 0), y_test)
        
        if self.final_model:
            def find_threshold_for_fpr(fpr, tpr, thresholds, target_fpr=0.2):
                idx = np.argmin(np.abs(fpr - target_fpr))
                return thresholds[idx], fpr[idx], tpr[idx]
            
            thresholds_df = pd.DataFrame(columns=["threshold", "FPR", "TPR"])
            target_fprs = [0.001, 0.005, 0.01, 0.05]

            for t_fpr in target_fprs:
                threshold, _, TPR = find_threshold_for_fpr(
                    self.fpr, self.tpr, thresholds, t_fpr
                )
                thresholds_df = thresholds_df.append({
                    "threshold": threshold, 
                    "FPR": t_fpr, 
                    "TPR": TPR,
                    'TPR_FPR_ratio': TPR / t_fpr
                }, ignore_index=True)
                
            thresholds_df.to_csv(
                '../data/' + self.base_dir + '/' + 
                self.full_name + '_thresholds.csv',
                index=False
            )
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()

        print(f"AUC on test set: {auroc:.3f}")
        print(f"Accuracy:        {acc:.3f}")
        
    def get_feat_imp(self):
        t1 = time.time()
        cur_msg = "Extracting most important features"
        l = ''.join(['.' for _ in range(60 - len(cur_msg))])
        sys.stdout.write(f"{cur_msg}{l}")
        sys.stdout.flush()
        
        if self.model_type == 'LGBM':
            feat_imp = self.best_model.feature_importances_
            feat_imp = feat_imp / feat_imp.max()
            self.feat_imp = pd.DataFrame({
                'Feature': self.feature_names, 
                'Importance_abs': feat_imp
            })
        elif self.model_type == 'LM':
            self.feat_imp = pd.DataFrame(self.feature_names, columns=['Feature'])
            self.feat_imp['Importance'] = self.best_model.coef_.T
            self.feat_imp['Sign'] = np.where(self.feat_imp['Importance'] > 0, "+", "-")
            self.feat_imp['Importance_abs'] = self.feat_imp['Importance'].abs()
            
        self.feat_imp.sort_values('Importance_abs', ascending=False, inplace=True)

        # match hex codes
        cols_hex_dict = pd.Series(
            self.df_colnames['hex'].values, 
            index=self.df_colnames['columns']
        ).to_dict()
        self.feat_imp['hex'] = self.feat_imp['Feature'].map(cols_hex_dict)

        self.feat_imp.to_csv(
            '../data/' + self.base_dir + '/' + 
            self.full_name + '_' +
            ('best_' + self.model_type + '_model_' if not self.final_model else '') + 
            'feature_importance.csv',
            index=False
        )
        pd.DataFrame({'TPR': self.tpr, 'FPR': self.fpr}).to_csv(
            '../data/' + self.base_dir + '/' + 
            self.full_name + '_' +
            ('best_' + self.model_type + '_model_' if not self.final_model else '') + 
            'ROC_curve.csv',
            index=False
        )
        
        total_time = time.time() - t1
        sys.stdout.write(f"DONE! -- {total_time:.2f} seconds\n")
        sys.stdout.flush()
        
    def get_feature_contributions(self):
        # find contribution of breakage scores in full table of features
        cols_group_dict = pd.Series(
            self.df_colnames['group'].values, 
            index=self.df_colnames['columns']
        ).to_dict()
        self.feat_imp['group'] = self.feat_imp['Feature'].map(cols_group_dict)

        all_coefs = self.feat_imp['Importance_abs'].sum()
        group_contr = self.feat_imp.groupby('group')['Importance_abs'].sum() / all_coefs * 100

        self.breakscore_contr = group_contr[['BREAKS', 'LONGRANGE_PREPARAMS']].sum()
        self.kmercount_contr = group_contr[['KMER_COUNTS']].sum()
        
        self.feat_imp.fillna('#787878', inplace=True)
        self.feat_imp.reset_index()
        
    def plot_results(self, num_features):
        feat_imp_subset = self.feat_imp.head(num_features)

        plt.figure(figsize=(8, 10))

        # Iterating over rows to plot each feature with its corresponding sign color
        if self.model_type == 'LGBM':
            palette = ['grey']
        elif self.model_type == 'LM':
            palette = ['#2166ac' if sign == '+' else '#b2182b' for sign in feat_imp_subset['Sign']]
        else:
            palette = None
        
        ax = sns.barplot(
            x='Importance_abs', 
            y='Feature', 
            data=feat_imp_subset, 
            palette=palette,
            alpha=0.8
        )

        # Coloring y-axis labels based on hex color in feat_imp_subset
        plt.yticks(range(len(feat_imp_subset)), feat_imp_subset['Feature'])

        # Coloring y-axis labels based on hex color in the table
        for index, label in enumerate(plt.gca().get_yticklabels()):
            hex_color = feat_imp_subset.iloc[index]['hex']
            label.set_color(hex_color)

        plt.xlabel('Coefficient (Absolute Value)')
        plt.ylabel('')
        plt.title(
            '''
                Feature importance from top {:.0f} features \n
                Contributions: Break score: {:.1f}%, kmer counts: {:.1f}%
            '''.format(num_features, self.breakscore_contr, self.kmercount_contr),
            loc='left',
            fontsize=10
        )

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        plt.savefig(
            '../figures/' + self.base_dir + '/' + 
            self.full_name + '_' +
            'best_' + self.model_type + '_' + 
            'FeatureImportance.pdf',
            dpi=600, bbox_inches='tight'
        )
        
    def plot_results_of_final_model(self, num_features):
        feat_imp_subset = self.feat_imp.head(num_features)

        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        # Plot AUROC Curve
        ax1 = axs[0]
        ax1.plot(self.fpr, self.tpr, color='darkorange', lw=2, label='ROC curve (area = 0.90)')
        ax1.plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlabel('False Positive Rate')
        ax1.set_ylabel('True Positive Rate')
        ax1.set_title('Receiver Operating Characteristic')
        ax1.legend(loc="lower right")

        # Plot Feature Importance (Horizontal Lollipop plot)
        ax2 = axs[1]
        y_positions = range(len(feat_imp_subset))
        ax2.hlines(y=y_positions, xmin=0, xmax=feat_imp_subset['importance'], color='grey')
        ax2.scatter(feat_imp_subset['importance'], y_positions, color='#787878', s=100)
        ax2.set_yticks(y_positions)
        ax2.set_yticklabels(feat_imp_subset['feature'])

        # Highest importance at the top
        ax2.invert_yaxis()
        ax2.set_title(f'Top {num_features} Scaled LightGBM Feature Importances')
        ax2.set_xlabel('Feature Importance Score (Scaled)')

        # Coloring y-axis labels based on hex color in the table
        for index, label in enumerate(plt.gca().get_yticklabels()):
            hex_color = feat_imp_subset.iloc[index]['hex']
            label.set_color(hex_color)

        for i in range(len(axs)):
            axs[i].spines['top'].set_visible(False)
            axs[i].spines['right'].set_visible(False)


        plt.tight_layout()
        fig.savefig(
            '../figures/' + self.base_dir + '/' + 
            self.full_name + '_Classifier.pdf',
            dpi=600, bbox_inches='tight'
        )
        
    def run_optuna_study(self):
        self.run_optuna_trials()
        self.run_optim_model()
        self.get_feat_imp()
        self.get_feature_contributions()