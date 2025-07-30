import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
from imblearn.over_sampling import SMOTE
from sklearn.preprocessing import PolynomialFeatures
import os
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
from torch.optim.swa_utils import AveragedModel, SWALR, update_bn
import datetime, time,random

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier

torch.manual_seed(42)
np.random.seed(42)

result_dir = ("./compare")
os.makedirs(result_dir, exist_ok=True)



class InceptionV2_1D(nn.Module):
    def __init__(self, in_channels, ch1x1, ch3x3red, ch3x3, ch5x5red, ch5x5, pool_proj):
        super().__init__()

        # Branch 1: 1x1 convolution
        self.branch1 = nn.Sequential(
            nn.Conv1d(in_channels, ch1x1, 1),
            nn.BatchNorm1d(ch1x1),
            nn.ReLU()
        )

        # Branch 2: 1x1 conv + 3x3 conv
        self.branch2 = nn.Sequential(
            nn.Conv1d(in_channels, ch3x3red, 1),
            nn.BatchNorm1d(ch3x3red),
            nn.ReLU(),
            nn.Conv1d(ch3x3red, ch3x3, 3, padding=1),
            nn.BatchNorm1d(ch3x3),
            nn.ReLU()
        )

        # Branch 3: 1x1 conv + 3x3 conv + 3x3 conv (factorized 5x5)
        self.branch3 = nn.Sequential(
            nn.Conv1d(in_channels, ch5x5red, 1),
            nn.BatchNorm1d(ch5x5red),
            nn.ReLU(),
            nn.Conv1d(ch5x5red, ch5x5red, 3, padding=1),
            nn.BatchNorm1d(ch5x5red),
            nn.ReLU(),
            nn.Conv1d(ch5x5red, ch5x5, 3, padding=1),
            nn.BatchNorm1d(ch5x5),
            nn.ReLU()
        )

        # Branch 4: MaxPool + 1x1 conv
        self.branch4 = nn.Sequential(
            nn.MaxPool1d(3, stride=1, padding=1),
            nn.Conv1d(in_channels, pool_proj, 1),
            nn.BatchNorm1d(pool_proj),
            nn.ReLU()
        )

    def forward(self, x):
        return torch.cat([
            self.branch1(x),
            self.branch2(x),
            self.branch3(x),
            self.branch4(x)
        ], 1)


# 使用Inception-v2的GoogleNet
class LR_GoogleV2_MLP(nn.Module):
    def __init__(self, input_dim, mlp_dims=[256, 128], dropout=0.4):
        super().__init__()

        # GoogleNet 1D with Inception-v2 modules
        self.cnn = nn.Sequential(
            nn.Conv1d(1, 64, 7, stride=2, padding=3),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(3, stride=2, padding=1),

            InceptionV2_1D(64, 64, 96, 128, 16, 32, 32),
            InceptionV2_1D(256, 128, 128, 192, 32, 96, 64),
            nn.MaxPool1d(3, stride=2, padding=1),

            InceptionV2_1D(480, 192, 96, 208, 16, 48, 64),
            InceptionV2_1D(512, 160, 112, 224, 24, 64, 64),
            InceptionV2_1D(512, 128, 128, 256, 24, 64, 64),
            InceptionV2_1D(512, 112, 144, 288, 32, 64, 64),
            InceptionV2_1D(528, 256, 160, 320, 32, 128, 128),
            nn.MaxPool1d(3, stride=2, padding=1)
        )

        self._get_cnn_output_dim(input_dim)

        # 增强的MLP分类器
        self.classifier = nn.Sequential(
            nn.Linear(self.cnn_out_features, mlp_dims[0]),
            nn.BatchNorm1d(mlp_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.BatchNorm1d(mlp_dims[1]),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.cnn(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)



class InceptionV3_1D(nn.Module):
    def __init__(self, in_channels, ch1x1, ch3x3red, ch3x3, ch3x3red2, ch3x3_1, ch3x3_2, pool_proj):
        super().__init__()

        # Branch 1: 1x1 convolution
        self.branch1 = nn.Sequential(
            nn.Conv1d(in_channels, ch1x1, 1),
            nn.BatchNorm1d(ch1x1),
            nn.ReLU()
        )

        # Branch 2: 1x1 conv + 3x3 conv
        self.branch2 = nn.Sequential(
            nn.Conv1d(in_channels, ch3x3red, 1),
            nn.BatchNorm1d(ch3x3red),
            nn.ReLU(),
            nn.Conv1d(ch3x3red, ch3x3, 3, padding=1),
            nn.BatchNorm1d(ch3x3),
            nn.ReLU()
        )

        # Branch 3: 1x1 conv + 3x3 conv + 3x3 conv (factorized)
        self.branch3 = nn.Sequential(
            nn.Conv1d(in_channels, ch3x3red2, 1),
            nn.BatchNorm1d(ch3x3red2),
            nn.ReLU(),
            nn.Conv1d(ch3x3red2, ch3x3_1, 3, padding=1),
            nn.BatchNorm1d(ch3x3_1),
            nn.ReLU(),
            nn.Conv1d(ch3x3_1, ch3x3_2, 3, padding=1),
            nn.BatchNorm1d(ch3x3_2),
            nn.ReLU()
        )

        # Branch 4: MaxPool + 1x1 conv
        self.branch4 = nn.Sequential(
            nn.MaxPool1d(3, stride=1, padding=1),
            nn.Conv1d(in_channels, pool_proj, 1),
            nn.BatchNorm1d(pool_proj),
            nn.ReLU()
        )

    def forward(self, x):
        return torch.cat([
            self.branch1(x),
            self.branch2(x),
            self.branch3(x),
            self.branch4(x)
        ], 1)


# 使用Inception-v3的GoogleNet
class LR_GoogleV3_MLP(nn.Module):
    def __init__(self, input_dim, mlp_dims=[256, 128], dropout=0.3):
        super().__init__()

        # GoogleNet 1D with Inception-v3 modules
        self.cnn = nn.Sequential(
            nn.Conv1d(1, 32, 3, stride=2, padding=1),  # Changed to 3x3 stride 2
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(32, 32, 3, padding=1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(32, 64, 3, padding=1),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(3, stride=2, padding=1),

            InceptionV3_1D(64, 64, 48, 64, 64, 96, 96, 32),
            InceptionV3_1D(256, 64, 48, 64, 64, 96, 96, 64),
            nn.MaxPool1d(3, stride=2, padding=1),

            InceptionV3_1D(288, 64, 48, 64, 64, 96, 96, 64),
            InceptionV3_1D(288, 64, 48, 64, 64, 96, 96, 64),
            InceptionV3_1D(288, 64, 48, 64, 64, 96, 96, 64),
            InceptionV3_1D(288, 64, 48, 64, 64, 96, 96, 64),
            InceptionV3_1D(288, 64, 48, 64, 64, 96, 96, 64),
            nn.MaxPool1d(3, stride=2, padding=1)
        )

        self._get_cnn_output_dim(input_dim)

        # MLP部分保持不变
        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(self.cnn_out_features, mlp_dims[0]),
            nn.Mish(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.Mish(inplace=True),
            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.cnn(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)



class CNN(nn.Module):
    def __init__(self, input_dim, num_channels=64, kernel_size=3, dropout=0.1):
        super().__init__()
        # 输入形状: (batch_size, 1, input_dim)
        self.conv_layers = nn.Sequential(
            nn.Conv1d(in_channels=1, out_channels=num_channels, kernel_size=kernel_size, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Dropout(dropout),

            nn.Conv1d(in_channels=num_channels, out_channels=num_channels * 2, kernel_size=kernel_size, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Dropout(dropout)
        )

        # 计算展平后的维度
        self._to_linear = None
        self._get_conv_output_dim(input_dim)

        self.classifier = nn.Sequential(
            nn.Linear(self._to_linear, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def _get_conv_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.conv_layers(dummy)
        self._to_linear = dummy.view(1, -1).shape[1]

    def forward(self, x):
        # 输入形状: (batch_size, input_dim)
        x = x.unsqueeze(1)  # 添加通道维度 (batch_size, 1, input_dim)
        x = self.conv_layers(x)
        x = x.view(x.size(0), -1)  # 展平
        return self.classifier(x)

class RNN(nn.Module):
    def __init__(self, input_dim, hidden_size=128, num_layers=2, dropout=0.1):
        super().__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers

        # 使用普通RNN替换LSTM
        self.rnn = nn.RNN(
            input_size=input_dim,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,
            dropout=dropout if num_layers > 1 else 0
        )
        self.dropout = nn.Dropout(dropout)

        self.classifier = nn.Sequential(
            nn.Linear(hidden_size, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        # 输入形状: (batch_size, input_dim)
        x = x.unsqueeze(1)  # (batch_size, 1, input_dim)

        # 初始化隐藏状态（仅需h0，无细胞状态c0）
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)

        # RNN前向传播（输出out和最终隐藏状态hn）
        out, hn = self.rnn(x, h0)  # out形状: (batch_size, 1, hidden_size)

        # 取最后一个时间步的输出
        out = out[:, -1, :]  # (batch_size, hidden_size)
        out = self.dropout(out)

        return self.classifier(out)

class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dims=[256,128], dropout=0.1):
        super().__init__()
        # 定义全连接层序列
        layers = []
        in_features = input_dim
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_features, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            in_features = h_dim

        # 分类头
        layers.append(nn.Linear(hidden_dims[-1], 1))
        layers.append(nn.Sigmoid())

        self.mlp = nn.Sequential(*layers)

    def forward(self, x):
        # 输入形状: (batch_size, input_dim)
        return self.mlp(x)  # 输出形状: (batch_size,1)

class CMLP(nn.Module):
    def __init__(self, input_dim, cnn_channels=64, mlp_dims=[256,128], dropout=0.1):
        super().__init__()
        # CNN部分
        self.cnn = nn.Sequential(
            nn.Conv1d(1, cnn_channels, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            nn.Conv1d(cnn_channels, cnn_channels * 2, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout)
        )

        # 计算CNN输出维度
        self._get_cnn_output_dim(input_dim)

        # MLP部分（修正输入维度）
        mlp_layers = []
        in_dim = self.cnn_out_features
        for h_dim in mlp_dims:
            mlp_layers += [
                nn.Linear(in_dim, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ]
            in_dim = h_dim
        mlp_layers.append(nn.Linear(mlp_dims[-1], 1))
        mlp_layers.append(nn.Sigmoid())

        self.mlp = nn.Sequential(*mlp_layers)

    def _get_cnn_output_dim(self, input_dim):
        # 计算实际输出维度
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)  # 正确计算展平维度

    def forward(self, x):
        # CNN处理
        x = x.unsqueeze(1)  # [B,1,D]
        x = self.cnn(x)  # [B,C,F]
        x = x.view(x.size(0), -1)  # [B, C*F]

        # MLP处理
        return self.mlp(x)

class CRNN(nn.Module):
    def __init__(self, input_dim, cnn_channels=64, rnn_hidden=128, dropout=0.2):
        super().__init__()
        # CNN部分保持不变
        self.cnn = nn.Sequential(
            nn.Conv1d(1, cnn_channels, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            nn.Conv1d(cnn_channels, cnn_channels * 2, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout)
        )

        # 修改维度计算方式
        self._get_cnn_output_dim(input_dim)

        # 修正RNN的input_size为通道数
        self.rnn = nn.RNN(
            input_size=self.cnn_out_channels,  # 使用CNN输出通道数
            hidden_size=rnn_hidden,
            num_layers=2,
            batch_first=True,
            dropout=dropout
        )

        # 分类器保持不变
        self.classifier = nn.Sequential(
            nn.Linear(rnn_hidden, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_channels = dummy.shape[1]  # 获取通道数 [B, C, F]
        self.seq_length = dummy.shape[2]  # 获取序列长度

    def forward(self, x):
        # CNN处理
        x = x.unsqueeze(1)  # [B, 1, D]
        x = self.cnn(x)  # [B, C, F]

        # 调整形状适应RNN
        x = x.permute(0, 2, 1)  # [B, F, C] (序列长度, 特征维度)

        # RNN处理
        out, _ = self.rnn(x)  # [B, F, H]
        out = out[:, -1, :]  # 取最后时间步 [B, H]

        return self.classifier(out)

class LR_CMLP(nn.Module):
    def __init__(self, input_dim, cnn_channels=64, mlp_dims=[256,128], dropout=0.1):
        super().__init__()
        # RF_CMLP更换导入数据
        # CNN部分
        self.cnn = nn.Sequential(
            nn.Conv1d(1, cnn_channels, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            nn.Conv1d(cnn_channels, cnn_channels * 2, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout)
        )

        # 计算CNN输出维度
        self._get_cnn_output_dim(input_dim)

        # MLP部分（修正输入维度）
        mlp_layers = []
        in_dim = self.cnn_out_features
        for h_dim in mlp_dims:
            mlp_layers += [
                nn.Linear(in_dim, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ]
            in_dim = h_dim
        mlp_layers.append(nn.Linear(mlp_dims[-1], 1))
        mlp_layers.append(nn.Sigmoid())

        self.mlp = nn.Sequential(*mlp_layers)

    def _get_cnn_output_dim(self, input_dim):
        # 计算实际输出维度
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)  # 正确计算展平维度

    def forward(self, x):
        # CNN处理
        x = x.unsqueeze(1)  # [B,1,D]
        x = self.cnn(x)  # [B,C,F]
        x = x.view(x.size(0), -1)  # [B, C*F]

        # MLP处理
        return self.mlp(x)




from sklearn.ensemble import RandomForestClassifier
def select_top_features(X_train, y_train, X_test):
    rf_model = RandomForestClassifier(
        n_estimators=290,  # 决策树数量
        max_depth=30,  # 树的最大深度
        min_samples_split=4,  # 分裂所需的最小样本数，默认2
        min_samples_leaf=2,  # 叶子节点所需的最小样本数，默认1
        random_state=42,
        n_jobs=1
    )
    rf_model.fit(X_train, y_train)
    feature_importances = rf_model.feature_importances_
    top_n_features = np.argsort(feature_importances)[-20:]  # 按重要性排序，取前 N 个特征
    # 筛选重要特征
    X_train_selected = X_train[:, top_n_features]
    X_test_selected = X_test[:, top_n_features]
    return X_train_selected, X_test_selected



class DetNetBottleneck(nn.Module):
    def __init__(self, in_channels, dilation=1):
        super().__init__()
        inter_channels = in_channels // 4
        self.conv1 = nn.Conv1d(in_channels, inter_channels, 1, bias=False)
        self.bn1 = nn.BatchNorm1d(inter_channels)

        self.conv2 = nn.Conv1d(
            inter_channels, inter_channels, 3,
            padding=dilation, dilation=dilation, bias=False)
        self.bn2 = nn.BatchNorm1d(inter_channels)

        self.conv3 = nn.Conv1d(inter_channels, in_channels, 1, bias=False)
        self.bn3 = nn.BatchNorm1d(in_channels)
        self.relu = nn.ReLU(inplace=True)

    def forward(self, x):
        identity = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)

        out += identity
        out = self.relu(out)
        return out


class LR_Des_CMLP(nn.Module):
    def __init__(self, input_dim, cnn_channels=32, mlp_dims=[256,128], dropout=0.2):
        super().__init__()

        # 修改后的CNN部分（加入DetNet）
        self.cnn = nn.Sequential(
            nn.Conv1d(1, cnn_channels, 3, padding='same'),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            DetNetBottleneck(cnn_channels, dilation=1),  # 添加DetNet模块
            DetNetBottleneck(cnn_channels, dilation=2),  # 添加多尺度空洞卷积

            nn.Conv1d(cnn_channels, cnn_channels * 2, 3, padding='same'),
            nn.BatchNorm1d(cnn_channels * 2),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            DetNetBottleneck(cnn_channels * 2, dilation=1),  # 第二个DetNet模块
        )

        # 计算CNN输出维度
        self._get_cnn_output_dim(input_dim)

        # MLP部分保持不变...
        mlp_layers = []
        in_dim = self.cnn_out_features
        for h_dim in mlp_dims:
            mlp_layers += [
                nn.Linear(in_dim, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ]
            in_dim = h_dim
        mlp_layers.append(nn.Linear(mlp_dims[-1], 1))
        mlp_layers.append(nn.Sigmoid())
        self.mlp = nn.Sequential(*mlp_layers)

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)  # [B,1,D]
        x = self.cnn(x)  # 经过包含DetNet的CNN
        x = x.view(x.size(0), -1)

        #MLP处理
        return self.mlp(x)

class Autoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim=64):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )

        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 256),
            nn.ReLU(),
            nn.Linear(256, input_dim)
        )

    def forward(self, x):
        z = self.encoder(x)
        reconstructed = self.decoder(z)
        return reconstructed

# 串联结构
class LR_AE_CMLP(nn.Module):
    def __init__(self, input_dim, cnn_channels=64, mlp_dims=[256,128], dropout=0.1, latent_dim=32):
        super().__init__()

        # 自编码器部分
        self.autoencoder = Autoencoder(input_dim, latent_dim)

        # CNN部分
        self.cnn = nn.Sequential(
            nn.Conv1d(1, cnn_channels, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            nn.Conv1d(cnn_channels, cnn_channels * 2, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout)
        )

        # 计算CNN输出维度
        self._get_cnn_output_dim(input_dim)

        # 联合特征维度
        combined_dim = self.cnn_out_features + latent_dim

        # MLP分类器
        mlp_layers = []
        in_dim = combined_dim
        for h_dim in mlp_dims:
            mlp_layers += [
                nn.Linear(in_dim, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ]
            in_dim = h_dim
        mlp_layers.append(nn.Linear(mlp_dims[-1], 1))
        mlp_layers.append(nn.Sigmoid())

        self.mlp = nn.Sequential(*mlp_layers)

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        # 自编码器重建
        reconstructed = self.autoencoder(x)

        # CNN处理
        x_cnn = x.unsqueeze(1)
        cnn_features = self.cnn(x_cnn).view(x.size(0), -1)

        # 自编码器潜在特征
        latent_features = self.autoencoder.encoder(x)

        # 特征融合
        combined = torch.cat([cnn_features, latent_features], dim=1)

        # 分类预测
        return self.mlp(combined), reconstructed



class LR_VGG_MLP(nn.Module):
    def __init__(self, input_dim, conv_channels=64, mlp_dims=[256,128], dropout=0.2):
        super().__init__()

        # VGG风格卷积块
        self.features = nn.Sequential(
            # Block 1: 2个卷积层 + 池化
            nn.Conv1d(1, conv_channels, 3, padding='same'),
            #nn.BatchNorm1d(conv_channels),
            nn.ReLU(),
            nn.Conv1d(conv_channels, conv_channels, 3, padding='same'),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

            # Block 2: 2个卷积层 + 池化
            nn.Conv1d(conv_channels, conv_channels * 2, 3, padding='same'),
            #nn.BatchNorm1d(conv_channels * 2),
            nn.ReLU(),
            nn.Conv1d(conv_channels * 2, conv_channels * 2, 3, padding='same'),
            #nn.BatchNorm1d(conv_channels * 2),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout),

        )

        # 计算输出维度
        self._get_cnn_output_dim(input_dim)

        # MLP部分
        mlp_layers = []
        in_dim = self.cnn_out_features
        for h_dim in mlp_dims:
            mlp_layers += [
                nn.Linear(in_dim, h_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ]
            in_dim = h_dim
        mlp_layers.append(nn.Linear(mlp_dims[-1], 1))
        mlp_layers.append(nn.Sigmoid())
        self.classifier = nn.Sequential(*mlp_layers)

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.features(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.features(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)

class PatchEmbedding(nn.Module):
    def __init__(self, input_dim=167, patch_size=16, embed_dim=128):
        super().__init__()
        self.num_patches = (input_dim + patch_size - 1) // patch_size
        self.patch_size = patch_size
        self.proj = nn.Linear(patch_size, embed_dim)
        self.cls_token = nn.Parameter(torch.randn(1, 1, embed_dim))
        self.pos_embed = nn.Parameter(torch.randn(1, self.num_patches + 1, embed_dim))

    def forward(self, x):
        # x shape: [B, 1, D]
        x = x.squeeze(1)  # [B, D]

        # Zero-padding if needed
        residual = x.size(1) % self.patch_size
        if residual != 0:
            x = F.pad(x, (0, self.patch_size - residual))

        # Create patches
        x = x.view(x.size(0), -1, self.patch_size)  # [B, num_patches, patch_size]
        x = self.proj(x)  # [B, num_patches, embed_dim]

        # Add cls token
        cls_tokens = self.cls_token.expand(x.size(0), -1, -1)
        x = torch.cat((cls_tokens, x), dim=1)  # [B, num_patches+1, embed_dim]

        # Add positional embedding
        x += self.pos_embed
        return x


class TransformerBlock(nn.Module):
    def __init__(self, embed_dim=128, num_heads=4, mlp_ratio=2):
        super().__init__()
        self.norm1 = nn.LayerNorm(embed_dim)
        self.attn = nn.MultiheadAttention(embed_dim, num_heads)
        self.norm2 = nn.LayerNorm(embed_dim)
        self.mlp = nn.Sequential(
            nn.Linear(embed_dim, embed_dim * mlp_ratio),
            nn.ReLU(),
            nn.Linear(embed_dim * mlp_ratio, embed_dim)
        )

    def forward(self, x):
        # x shape: [B, N+1, D]
        x = x.permute(1, 0, 2)  # [N+1, B, D]
        attn_out, _ = self.attn(x, x, x)
        x = x + attn_out
        x = self.norm1(x)

        mlp_out = self.mlp(x)
        x = x + mlp_out
        x = self.norm2(x)
        return x.permute(1, 0, 2)


import math
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(max_len, d_model)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe.unsqueeze(0))

    def forward(self, x):
        x = x + self.pe[:, :x.size(1), :]
        return self.dropout(x)


class LR_Transformer(nn.Module):
    def __init__(self, input_dim, d_model=128, nhead=8, num_layers=3, dim_feedforward=256, dropout=0.1):
        super().__init__()
        self.d_model = d_model
        self.input_embedding = nn.Linear(input_dim, d_model)
        self.pos_encoder = PositionalEncoding(d_model, dropout)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model, nhead=nhead, dim_feedforward=dim_feedforward,
            dropout=dropout, batch_first=True
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

        self.classifier = nn.Sequential(
            nn.Linear(d_model, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.input_embedding(x) * math.sqrt(self.d_model)
        x = self.pos_encoder(x)
        x = self.transformer_encoder(x)
        x = x[:, 0, :]
        return self.classifier(x)



class LR_ViT(nn.Module):
    def __init__(self, input_dim, mlp_dims=[256,128],
                 embed_dim=128, num_layers=4, num_heads=4,
                 patch_size=64, dropout=0.1):
        super().__init__()

        # Vision Transformer
        self.patch_embed = PatchEmbedding(input_dim, patch_size, embed_dim)

        self.transformer = nn.Sequential(*[
            TransformerBlock(embed_dim, num_heads)
            for _ in range(num_layers)
        ])

        # MLP Head
        self.mlp_head = nn.Sequential(
            nn.LayerNorm(embed_dim),
            nn.Linear(embed_dim, mlp_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        # Input shape: [B, 1, D]
        x = self.patch_embed(x)

        # Transformer
        x = self.transformer(x)  # [B, N+1, D]

        # Use cls token for classification
        cls_token = x[:, 0, :]  # [B, D]

        # MLP Head
        return self.mlp_head(cls_token)



class LR_Alex_MLP(nn.Module):
    def __init__(self, input_dim, mlp_dims=[256, 128], dropout=0.5):
        super().__init__()

        # AlexNet 1D版本（调整参数适配短序列）
        self.cnn = nn.Sequential(
            nn.Conv1d(1, 64, 5, stride=2, padding=2),  # 输出尺寸: (L-5+4)//2 +1
            nn.ReLU(),
            nn.MaxPool1d(3, stride=2, padding=1),  # 输出尺寸: (L-3+2)//2 +1

            nn.Conv1d(64, 192, 3, padding=1),
            nn.ReLU(),
            nn.MaxPool1d(3, stride=2, padding=1),

            nn.Conv1d(192, 384, 3, padding=1),
            nn.ReLU(),

            nn.Conv1d(384, 256, 3, padding=1),
            nn.ReLU(),

            nn.Conv1d(256, 256, 3, padding=1),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1)  # 自适应池化避免尺寸错误
        )

        self._get_cnn_output_dim(input_dim)

        # MLP部分
        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(256, mlp_dims[0]),  # 直接使用256作为输入维度
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.ReLU(),
            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.cnn(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)

class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, downsample=None):
        super().__init__()
        self.conv1 = nn.Conv1d(in_channels, out_channels, kernel_size=3,
                               stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv1d(out_channels, out_channels, kernel_size=3,
                               stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm1d(out_channels)
        self.downsample = downsample
        self.dropout = nn.Dropout(0.3)  # 增加dropout防止过拟合

    def forward(self, x):
        identity = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)
        out = self.dropout(out)

        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            identity = self.downsample(x)

        out += identity
        out = self.relu(out)
        return out


class LR_Res_MLP(nn.Module):
    def __init__(self, input_dim, init_channels=32, mlp_dims=[128, 64], dropout=0.3):
        super().__init__()

        # 初始卷积层
        self.features = nn.Sequential(
            nn.Conv1d(1, init_channels, 3, padding='same', bias=False),
            nn.BatchNorm1d(init_channels),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(dropout)
        )

        # 残差块配置（根据小样本调整结构）
        self.layer1 = self._make_layer(init_channels, init_channels, stride=1)
        self.layer2 = self._make_layer(init_channels, init_channels * 2, stride=1)

        # 自适应池化替代固定池化
        self.avgpool = nn.AdaptiveAvgPool1d(4)

        # 动态计算特征维度
        self._get_cnn_output_dim(input_dim)

        # 精简MLP结构
        self.classifier = nn.Sequential(
            nn.Linear(self.cnn_out_features, mlp_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def _make_layer(self, in_channels, out_channels, stride):
        downsample = None
        if stride != 1 or in_channels != out_channels:
            downsample = nn.Sequential(
                nn.Conv1d(in_channels, out_channels, 1, stride=stride, bias=False),
                nn.BatchNorm1d(out_channels)
            )
        return nn.Sequential(
            ResidualBlock(in_channels, out_channels, stride, downsample),
            ResidualBlock(out_channels, out_channels)  # 每个阶段两个残差块
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.features(dummy)
        dummy = self.layer1(dummy)
        dummy = self.layer2(dummy)
        dummy = self.avgpool(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.features(x)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.avgpool(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)



class Inception1D(nn.Module):
    def __init__(self, in_channels, ch1x1, ch3x3red, ch3x3, ch5x5red, ch5x5, pool_proj,
                 use_attention=True,
                 reduction=16,  # 通道注意力降维比例
                 kernel_size=7
                 ):
        super().__init__()

        #self.use_attention = use_attention

        self.branch1 = nn.Sequential(
            nn.Conv1d(in_channels, ch1x1, 1),
            nn.BatchNorm1d(ch1x1),
            nn.ReLU()
        )

        self.branch2 = nn.Sequential(
            nn.Conv1d(in_channels, ch3x3red, 1),
            nn.BatchNorm1d(ch3x3red),
            nn.ReLU(),
            nn.Conv1d(ch3x3red, ch3x3, 3, padding=1),
            nn.BatchNorm1d(ch3x3),
            nn.ReLU()
        )

        self.branch3 = nn.Sequential(
            nn.Conv1d(in_channels, ch5x5red, 1),
            nn.BatchNorm1d(ch5x5red),
            nn.ReLU(),
            nn.Conv1d(ch5x5red, ch5x5, 5, padding=2),
            nn.BatchNorm1d(ch5x5),
            nn.ReLU()
        )

        self.branch4 = nn.Sequential(
            nn.MaxPool1d(3, stride=1, padding=1),
            nn.Conv1d(in_channels, pool_proj, 1),
            nn.BatchNorm1d(pool_proj),
            nn.ReLU()
        )



    def forward(self, x):
        return torch.cat([
            self.branch1(x),
            self.branch2(x),
            self.branch3(x),
            self.branch4(x)
        ], 1)

class LR_Google_mlp(nn.Module):
    def __init__(self, input_dim, mlp_dims=[256,128], dropout=0.2):
        super().__init__()

        # GoogleNet 1D版本
        self.cnn = nn.Sequential(
            nn.Conv1d(1, 64, 7, stride=2, padding=3),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(3, stride=2, padding=1),

            Inception1D(64, 64, 96, 128, 16, 32, 32,reduction=16),
            Inception1D(256, 128, 128, 192, 32, 96, 64,reduction=16),
            nn.MaxPool1d(3, stride=2, padding=1),

            Inception1D(480, 192, 96, 208, 16, 48, 64),
            Inception1D(512, 160, 112, 224, 24, 64, 64),
            Inception1D(512, 128, 128, 256, 24, 64, 64),
            Inception1D(512, 112, 144, 288, 32, 64, 64),
            Inception1D(528, 256, 160, 320, 32, 128, 128),
            nn.MaxPool1d(3, stride=2, padding=1)
        )

        self._get_cnn_output_dim(input_dim)

        # MLP部分
        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(self.cnn_out_features, mlp_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(mlp_dims[0], mlp_dims[1]),
            nn.ReLU(),
            nn.Linear(mlp_dims[1], 1),
            nn.Sigmoid()
        )

    def _get_cnn_output_dim(self, input_dim):
        dummy = torch.zeros(1, 1, input_dim)
        dummy = self.cnn(dummy)
        self.cnn_out_features = dummy.view(1, -1).size(1)

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.cnn(x)
        x = x.view(x.size(0), -1)
        return self.classifier(x)


#ML
rf_model = RandomForestClassifier(
    n_estimators=290,  # 决策树数量
    max_depth=30,  # 树的最大深度
    min_samples_split=4,  # 分裂所需的最小样本数，默认2
    min_samples_leaf=2,  # 叶子节点所需的最小样本数，默认1
    random_state=42,
    n_jobs=1
)

svm_model = SVC(
    C=0.03,  # 正则化参数，越小正则化越强，默认1.0
    kernel='linear',  # 核函数，'linear'、'rbf'、'poly'等，默认'rbf'
    gamma='scale',  # 核函数的系数，'scale'或'auto'，默认'scale'
    probability=True,  # 是否启用概率估计，默认False
    random_state=42  # 确保结果可重复
)

xgb_model = XGBClassifier(
    n_estimators=290,
    max_depth=30,
    learning_rate=0.1,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    n_jobs=1,
    eval_metric='logloss'
)


def calculate_metrics(y_true, y_pred, y_prob):
    return {
        'accuracy': accuracy_score(y_true, y_pred),
        'f1': f1_score(y_true, y_pred),
        'auc': roc_auc_score(y_true, y_prob),
        'cm': confusion_matrix(y_true, y_pred)
    }


def train_model(model, train_loader, val_loader, criterion, optimizer, scheduler, epochs=100):
    history = {'train_loss': [], 'val_loss': [], 'train_metrics': [], 'val_metrics': []}
    best_auc = 0
    best_model = None

    for epoch in range(epochs):
        model.train()
        epoch_train_loss = 0
        all_preds, all_probs, all_labels = [], [], []

        for inputs, labels in train_loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels.unsqueeze(1))
            loss.backward()
            optimizer.step()
            epoch_train_loss += loss.item()

            probs = outputs.squeeze().detach().numpy()
            preds = (probs > 0.5).astype(int)
            if np.isscalar(probs):
                probs = [probs]
                preds = [preds]  # 这里也要添加
            all_preds.extend(preds)
            all_probs.extend(probs)
            all_labels.extend(labels.numpy())

        train_metrics = calculate_metrics(all_labels, all_preds, all_probs)

        # 修改这里：只取前两个返回值
        val_loss, val_metrics, _, _ = evaluate_loader(model, val_loader, criterion)

        history['train_loss'].append(epoch_train_loss / len(train_loader))
        history['val_loss'].append(val_loss)
        history['train_metrics'].append(train_metrics)
        history['val_metrics'].append(val_metrics)

        # 更新调度器
        # scheduler.step(val_metrics['auc'])

        if val_metrics['auc'] > best_auc:
            best_auc = val_metrics['auc']
            best_model = model.state_dict().copy()

    if best_model is not None:
        model.load_state_dict(best_model)
    return model, history


# 修改evaluate_loader函数
def evaluate_loader(model, loader, criterion):
    model.eval()
    loss = 0
    all_preds, all_probs, all_labels = [], [], []

    with torch.no_grad():
        for inputs, labels in loader:
            outputs = model(inputs)
            loss += criterion(outputs, labels.unsqueeze(1)).item()

            probs = outputs.squeeze().numpy()
            preds = (probs > 0.5).astype(int)
            if np.isscalar(probs):
                probs = [probs]
                preds = [preds]  # 这里也要添加
            all_preds.extend(preds)
            all_probs.extend(probs)
            all_labels.extend(labels.numpy())

    metrics = calculate_metrics(all_labels, all_preds, all_probs)
    return loss / len(loader), metrics, np.array(all_labels), np.array(all_probs)

def plot_loss_curves(history, fold, result_dir):
    plt.figure(figsize=(12, 6))
    plt.plot(history['train_loss'], label='Train Loss')
    plt.plot(history['val_loss'], label='Validation Loss')
    plt.title(f'Fold {fold} Training/Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(os.path.join(result_dir, f'loss_curve_fold_{fold}.png'))
    plt.close()


def plot_roc_curves(all_val_true, all_val_probs, all_test_true, all_test_probs, result_dir, png_name):
    """
    绘制验证集和测试集的ROC曲线
    """
    # 创建图形
    plt.figure(figsize=(14, 6))

    # 验证集ROC曲线
    plt.subplot(1, 2, 1)
    for i, (val_true, val_probs) in enumerate(zip(all_val_true, all_val_probs)):
        fpr, tpr, _ = roc_curve(val_true, val_probs)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=3, alpha=0.5, label=f'Run {i + 1} (AUC={roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Validation Set ROC Curves (30 Runs)')


    # 测试集ROC曲线
    plt.subplot(1, 2, 2)
    for i, (test_true, test_probs) in enumerate(zip(all_test_true, all_test_probs)):
        fpr, tpr, _ = roc_curve(test_true, test_probs)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=3, alpha=0.5, label=f'Run {i + 1} (AUC={roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Test Set ROC Curves (30 Runs)')


    plt.tight_layout()
    plt.savefig(os.path.join(result_dir, png_name))
    plt.close()


def load_data_label(train_csv_path, test_csv_path):
    train_data = pd.read_csv(train_csv_path)
    test_data = pd.read_csv(test_csv_path)

    X_train = train_data.iloc[:, 3:].values.astype(np.float32)
    y_train = train_data['label'].values.astype(np.float32)
    X_test = test_data.iloc[:, 3:].values.astype(np.float32)
    y_test = test_data['label'].values.astype(np.float32)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)


    return X_train, y_train, X_test, y_test




def main():
    start_time = time.time()  # 记录程序开始时间

    train_path = './data/train167.csv'
    test_path = './data/test131.csv'

    X_train_orig, y_train_orig, X_test, y_test = load_data_label(train_path, test_path)
    X_train_orig, X_test = select_top_features(X_train_orig, y_train_orig, X_test)


    # 存储所有结果
    results = []
    # 存储所有运行的验证集和测试集的真实标签和预测概率
    all_val_true = []
    all_val_probs = []
    all_test_true = []
    all_test_probs = []
    '''
    model_names=['LR','RF','SVM','XGBoost','CNN','RNN','MLP','CMLP','RCNN',
                 'LR_CMLP','RF_CMLP','LR_AE_CMLP','LR_Des_CMLP','LR_VGG_MLP',
                 'LR_Transformer','LR_ViT','LR_Alex_MLP','LR_Res_MLP','LR_Google_MLP',
                 'LR_GoogleV2_MLP','LR_AttenGoogleV2_MLP','LR_GoogleV3_MLP']'''

    for seed_index in range(0,30):  # 随机种子
        seed = 9*seed_index + 480
        # print(f"\n=== Run {seed_index+1}/10 (Seed {seed}) ===")

        # 重置随机数种子
        torch.manual_seed(seed)
        np.random.seed(seed)

        # 应用SMOTE上采样
        smote = SMOTE(sampling_strategy='auto', k_neighbors=15, random_state=42)
        X_train_res, y_train_res = smote.fit_resample(X_train_orig, y_train_orig)

        # 创建原始样本标记（前N个为True）
        n_original = X_train_orig.shape[0]
        is_original = np.array([True] * n_original + [False] * (len(X_train_res) - n_original))

        # 8:2划分训练集和验证集
        train_orig_idx, val_orig_idx = train_test_split(
            np.arange(n_original),
            test_size=0.2,
            stratify=y_train_orig,
            random_state=seed
        )

        # 构建增强后的训练集（原始训练样本 + 所有虚拟样本）
        virtual_idx = np.where(is_original == False)[0]
        train_res_idx = np.concatenate([train_orig_idx, virtual_idx])

        # 验证集使用原始验证样本
        X_train_fold = X_train_res[train_res_idx]
        y_train_fold = y_train_res[train_res_idx]
        X_val_fold = X_train_res[val_orig_idx]  # 原始验证样本
        y_val_fold = y_train_res[val_orig_idx]

        # 标准化处理（在划分后执行）
        scaler = StandardScaler()
        X_train_fold = scaler.fit_transform(X_train_fold)
        X_val_fold = scaler.transform(X_val_fold)
        X_test_scaled = scaler.transform(X_test)
        X_train_orig_scaled = scaler.transform(X_train_orig)  # 使用相同的scaler

        # 创建数据加载器
        train_loader = DataLoader(
            TensorDataset(torch.FloatTensor(X_train_fold), torch.FloatTensor(y_train_fold)),
            batch_size=64, shuffle=True
        )
        val_loader = DataLoader(
            TensorDataset(torch.FloatTensor(X_val_fold), torch.FloatTensor(y_val_fold)),
            batch_size=64
        )
        test_loader = DataLoader(
            TensorDataset(torch.FloatTensor(X_test_scaled), torch.FloatTensor(y_test)),
            batch_size=64
        )
        # 新增：整个训练集的加载器（仅原始样本）
        train_orig_loader = DataLoader(
            TensorDataset(torch.FloatTensor(X_train_orig_scaled), torch.FloatTensor(y_train_orig)),
            batch_size=64
        )

        # 初始化模型
        # 测试不同模型时替换模型名
        model = LR_CMLP(input_dim=X_train_orig.shape[1])
        criterion = nn.BCELoss()
        optimizer = optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, mode='max', factor=0.5, patience=5
        )

        # 训练模型
        model, history = train_model(
            model, train_loader, val_loader,
            criterion, optimizer, scheduler, epochs=200
        )


        # 验证集评估（仅原始样本）
        _, val_metrics, val_true, val_probs = evaluate_loader(model, val_loader, criterion)
        best_val_auc = max([m['auc'] for m in history['val_metrics']])

        # 测试集评估
        _, test_metrics, test_true, test_probs = evaluate_loader(model, test_loader, criterion)
        test_auc = test_metrics['auc']

        # 新增：整个训练集评估（仅原始样本）
        _, train_metrics, _, _ = evaluate_loader(model, train_orig_loader, criterion)
        train_auc = train_metrics['auc']

        # 保存本次运行的验证集和测试集结果
        all_val_true.append(val_true)
        all_val_probs.append(val_probs)
        all_test_true.append(test_true)
        all_test_probs.append(test_probs)

        # 存储结果
        results.append({
            'seed': seed,
            'train_auc': train_auc,
            'val_auc': best_val_auc,
            'test_auc': test_auc
        })

        #print(f"Run {seed_index + 1} Results:")
        print(f"Train AUC: {train_auc:.4f}")
        print(f"Best Val AUC: {best_val_auc:.4f}")
        print(f"Test AUC: {test_auc:.4f}")



    pd.DataFrame(results).to_excel(os.path.join(result_dir, "seed5_results.xlsx"), index=False)
    # 新增：绘制ROC曲线
    plot_roc_curves(all_val_true, all_val_probs, all_test_true, all_test_probs, result_dir,png_name="roc_curves_5.png")


    # 在程序结束时打印时间
    end_time = time.time()  # 记录程序结束时间
    run_time = end_time - start_time  # 计算运行时长

    # 打印结束时间和运行时长
    current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
    print(f"程序结束时间: {current_time}")
    print(f"运行时长: {run_time:.2f} 秒")





if __name__ == "__main__":
    main()




