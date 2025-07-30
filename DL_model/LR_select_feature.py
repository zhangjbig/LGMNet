import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import math, joblib
import os
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from skopt import BayesSearchCV
from skopt.space import Integer, Categorical, Real

torch.manual_seed(42)
np.random.seed(42)

result_dir = f"./results"
os.makedirs(result_dir, exist_ok=True)


def load_data_label(train_csv_path, test_csv_path):
    # 加载训练集和测试集
    train_data = pd.read_csv(train_csv_path)
    test_data = pd.read_csv(test_csv_path)
    # 分离特征和标签
    X_train = train_data.iloc[:, 3:]  # 第4列至最后一列是特征
    y_train = train_data['label']  # 第3列是标签
    X_test = test_data.iloc[:, 3:]  # 测试集特征
    y_test = test_data['label']  # 测试集标签

    # 特征标准化
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    joblib.dump(scaler, "trained_scaler1.pkl")

    return X_train, y_train, X_test, y_test


def select_top_features(X_train, y_train, X_test):
    lr_model = LogisticRegression(
        penalty='l2',
        C=1.31,
        solver='lbfgs',
        max_iter=5000,
        random_state=42
    )
    lr_model.fit(X_train, y_train)
    importances = np.abs(lr_model.coef_[0])
    sorted_indices = np.argsort(importances)[::-1]
    top_n_indices = sorted_indices[:32]
    print(top_n_indices)
    X_train_selected = X_train[:, top_n_indices]
    X_test_selected = X_test[:, top_n_indices]
    return X_train_selected, X_test_selected

def save_csv(X_train_selected,y_train,X_test_selected,y_test):
    # 保存选出的特征和标签到 CSV 文件
    train_df = pd.DataFrame(X_train_selected)
    train_df.insert(0, 'label', y_train)
    train_df.to_csv(os.path.join('data/selected_train167.csv'), index=True)

    test_df = pd.DataFrame(X_test_selected)
    test_df.insert(0, 'label', y_test)
    test_df.to_csv(os.path.join('data/selected_test131.csv'), index=True)

'''def create_data_loaders(X_train_selected, y_train, X_test_selected, y_test, test_size=0.2, random_state=42):
    X_train_selected = torch.FloatTensor(X_train_selected)
    y_train = torch.FloatTensor(y_train)
    X_test_selected = torch.FloatTensor(X_test_selected)
    y_test = torch.FloatTensor(y_test)

    X_train, X_val, y_train, y_val = train_test_split(X_train_selected, y_train, test_size=test_size,
                                                      random_state=random_state)

    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=32, shuffle=True)
    val_loader = DataLoader(TensorDataset(X_val, y_val), batch_size=32, shuffle=False)
    test_loader = DataLoader(TensorDataset(X_test_selected, y_test), batch_size=32, shuffle=False)

    return train_loader, val_loader, test_loader
'''


train_path = './data/train167.csv'
test_path = './data/test131.csv'

X_train, y_train, X_test, y_test = load_data_label(train_path, test_path)
X_train_selected, X_test_selected = select_top_features(X_train, y_train, X_test)
save_csv(X_train_selected,y_train,X_test_selected,y_test)
