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
import datetime, time,random
import joblib

torch.manual_seed(42)
np.random.seed(42)

result_dir = ("./results")
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
class GoogleNetV2(nn.Module):
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
            all_preds.extend(preds)
            all_probs.extend(probs)
            all_labels.extend(labels.numpy())

        train_metrics = calculate_metrics(all_labels, all_preds, all_probs)
        val_loss, val_metrics = evaluate_loader(model, val_loader, criterion)

        history['train_loss'].append(epoch_train_loss / len(train_loader))
        history['val_loss'].append(val_loss)
        history['train_metrics'].append(train_metrics)
        history['val_metrics'].append(val_metrics)

        # 更新调度器（根据验证集的 AUC 调整学习率）
        #scheduler.step(val_metrics['auc'])  # 若用验证损失，改为 scheduler.step(val_loss)

        if val_metrics['auc'] > best_auc:
            best_auc = val_metrics['auc']
            best_model = model.state_dict().copy()

    if best_model is not None:
        model.load_state_dict(best_model)
    return model, history


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
            all_preds.extend(preds)
            all_probs.extend(probs)
            all_labels.extend(labels.numpy())

    metrics = calculate_metrics(all_labels, all_preds, all_probs)
    return loss / len(loader), metrics


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


def load_data_label(train_csv_path, test_csv_path):
    train_data = pd.read_csv(train_csv_path)
    test_data = pd.read_csv(test_csv_path)

    X_train = train_data.iloc[:, 2:].values.astype(np.float32)
    y_train = train_data['label'].values.astype(np.float32)
    X_test = test_data.iloc[:, 2:].values.astype(np.float32)
    y_test = test_data['label'].values.astype(np.float32)

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    #joblib.dump(scaler, "trained_scaler2.pkl")
    X_test = scaler.transform(X_test)


    return X_train, y_train, X_test, y_test





def main():
    start_time = time.time()  # 记录程序开始时间

    train_path = './data/selected_train167.csv'
    test_path = './data/selected_test131.csv'

    X_train_orig, y_train_orig, X_test, y_test = load_data_label(train_path, test_path)


    # 存储所有结果
    results = []

    for seed in range(0,30):  # 随机种子
        seed = 9*seed + 480
        print(f"\n=== Random Seed {seed} ===")

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
        #joblib.dump(scaler, "trained_scaler3.pkl")
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
        model = GoogleNetV2(input_dim=X_train_orig.shape[1])
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

        torch.save(model.state_dict(), f"best_model_{seed}.pth")

        # 验证集评估（仅原始样本）
        _, val_metrics = evaluate_loader(model, val_loader, criterion)
        best_val_auc = max([m['auc'] for m in history['val_metrics']])

        # 测试集评估
        _, test_metrics = evaluate_loader(model, test_loader, criterion)
        test_auc = test_metrics['auc']

        # 整个训练集评估（仅原始样本）
        _, train_metrics = evaluate_loader(model, train_orig_loader, criterion)
        train_auc = train_metrics['auc']

        # 存储结果（新增train_auc）
        results.append({
            'seed': seed,
            'train_auc': train_auc,  # 新增：整个训练集结果
            'val_auc': best_val_auc,  # 内部验证集结果
            'test_auc': test_auc  # 外部测试集结果
        })

        print(f"Seed {seed} Results:")
        print(f"Train AUC: {train_auc:.4f}")  # 新增打印
        print(f"Best Val AUC: {best_val_auc:.4f}")
        print(f"Test AUC: {test_metrics['auc']:.4f}")

    pd.DataFrame(results).to_excel(os.path.join(result_dir, "seed30_results.xlsx"), index=False)


    # 在程序结束时打印时间
    end_time = time.time()  # 记录程序结束时间
    run_time = end_time - start_time  # 计算运行时长

    # 打印结束时间和运行时长
    current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
    print(f"程序结束时间: {current_time}")
    print(f"运行时长: {run_time:.2f} 秒")





if __name__ == "__main__":
    main()




