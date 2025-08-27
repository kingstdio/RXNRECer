"""
Author: Zhenkun Shi
Date: 2023-10-05 10:54:13
LastEditors: Zhenkun Shi
LastEditTime: 2023-10-06 14:38:55
FilePath: /preaction/methods/alfp_end2end.py
Description:

Copyright (c) 2023 by tibd, All Rights Reserved.
"""


from cProfile import label
import os
import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import pandas as pd
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
import esm
from sklearn.model_selection import train_test_split
from tkinter import _flatten  # type: ignore
sys.path.insert(0, os.path.dirname(os.path.realpath("__file__")))
sys.path.insert(1, "../")



# Step 1: Define the Dataset and DataLoader
class ProteinDataset(torch.utils.data.Dataset):
    def __init__(self, sequences, labels):
        self.sequences = sequences
        self.labels = labels

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        # Convert sequences and labels to PyTorch tensors
        sequence = self.sequences[idx]
        label = self.labels[idx]
        label = torch.from_numpy(label).float()
        #sequence  = torch.from_numpy(sequence )
        return sequence, label

import time
Global_File_Name = "time_consuming.txt"
# 定义一个装饰器来测试方法执行时间并记录到文件
def measure_execution_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        with open(Global_File_Name, "a") as file:
            file.write(f"Method '{func.__name__}' execution time: {execution_time} seconds\n")
        return result
    return wrapper

# region Step 2: Define the EsmEmbedding with layer freezing
class EsmEmbedding(nn.Module):
    def __init__(self, device, freeze_esm_layers):
        super(EsmEmbedding, self).__init__()
        self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.freeze_esm_layers = freeze_esm_layers
        self.device = device

        # Freeze the specified layers
        for name, param in self.esm_model.named_parameters():
            if any(layer_name in name for layer_name in [f"layer.{i}." for i in range(self.freeze_esm_layers)]):
                param.requires_grad = False
        self.esm_model = self.esm_model.to(self.device)

    def forward(self, x):
        batch_converter = self.esm_alphabet.get_batch_converter()
        seq_idx = [f"seq_{i}" for i in range(len(x))]
        
        x_esm_input = list(zip(seq_idx, x))
        batch_labels, batch_strs, batch_tokens = batch_converter(x_esm_input)
        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)

        batch_tokens = batch_tokens.to(self.device)

        embd_results = self.esm_model(batch_tokens, repr_layers=[33], return_contacts=False)

        token_representations = embd_results["representations"][33]
        sequence_representations = []
        for i, tokens_len in enumerate(batch_lens):
            sequence_representations.append(token_representations[i, 1 : tokens_len - 1].mean(0))
        
        sequence_representations = torch.stack(sequence_representations, dim=0)
        return sequence_representations.to(self.device)

# endregion


# region self-attention
class SelfAttention(nn.Module):
    def __init__(self, attention_size, input_dimensions):
        super(SelfAttention, self).__init__()
        self.attention_size = attention_size

        self.W = nn.Parameter(torch.Tensor(input_dimensions, self.attention_size))
        self.b = nn.Parameter(torch.Tensor(1, self.attention_size))
        self.u = nn.Parameter(torch.Tensor(self.attention_size, 1))

        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_normal_(self.W)
        nn.init.zeros_(self.b)
        nn.init.xavier_normal_(self.u)

    def forward(self, x, mask=None):
        # x: (BATCH_SIZE, MAX_TIMESTEPS, EMBED_SIZE)
        # et: (BATCH_SIZE, MAX_TIMESTEPS, ATTENTION_SIZE)
        et = torch.tanh(torch.matmul(x, self.W) + self.b)

        # at: (BATCH_SIZE, MAX_TIMESTEPS)
        at = torch.softmax(torch.squeeze(torch.matmul(et, self.u), dim=-1), dim=-1)

        if mask is not None:
            at *= mask.float()

        # ot: (BATCH_SIZE, MAX_TIMESTEPS, EMBED_SIZE)
        atx = at.unsqueeze(dim=-1)
        ot = atx * x

        # output: (BATCH_SIZE, EMBED_SIZE)
        output = torch.sum(ot, dim=1)

        return output
# endregion


# region Step 4: Define the Model
class BGRU(nn.Module):
    def __init__(self, input_dimensions, device, gru_h_size=512, attention_size=32, dropout=0.2, output_dimensions=300, freeze_esm_layers=23):
        super(BGRU, self).__init__()
        self.device = device
        # Create an EsmEmbedding instance
        self.embedding = EsmEmbedding(freeze_esm_layers=freeze_esm_layers, device=self.device)
        # Define the GRU, attention, and output layers
        self.gru = nn.GRU(input_dimensions, gru_h_size, batch_first=True, bidirectional=True)
        self.attention = SelfAttention(attention_size, (gru_h_size * 2))
        self.output_layer = nn.Linear((gru_h_size * 2), output_dimensions)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        esm_out = self.embedding(x)  # ESM embedding
        esm_out = esm_out.unsqueeze(1)
        gru_output, _ = self.gru(esm_out)
        attention_out = self.attention(gru_output)
        # print(torch.unique(torch.eq(gru_output.squeeze(1), attention_out)))
        drop_out_x = self.dropout(attention_out)  # Apply dropout
        output = self.output_layer(drop_out_x)
        # output1 = attention_out - gru_output
        output1 = self.output_layer(gru_output.squeeze(1))
        return output,output1

@measure_execution_time
def select_samples(train_sequences,train_labels,model,k,batch_size):
    if train_sequences.shape[0]<k:
        return list(range(train_sequences.shape[0]))
    train_dataset = ProteinDataset(train_sequences, train_labels)
    loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=False)
    sums = []
    model.eval()
    with torch.no_grad():
        for batch_x_squences, batch_y_labels in loader:
            outputs = model(batch_x_squences)[0]
            probabilities = torch.sigmoid(outputs).clamp(min=1e-10, max=1-1e-10)
            # row_sum = torch.sum(torch.abs(outputs[1] - outputs[0]), dim=1)
            sums.append(torch.abs(probabilities - 0.5).sum(dim=1))
    # 将列表中的张量拼接成一个张量
    result = torch.cat(sums, dim=0)
    # 使用 torch.topk 函数获取最大值及其索引
    topk_values, topk_indices = torch.topk(result, k)
    return topk_indices

import random
@measure_execution_time
def select_random(train_sequences,train_labels,model,k,batch_size):
    if train_sequences.shape[0]<k:
        return list(range(train_sequences.shape[0]))
    # 假设 all_indices 是包含所有索引的列表
    all_indices = list(range(train_sequences.shape[0]))

    # 从 all_indices 中随机抽取 k 个索引
    k = 100  # 假设要抽取 100 个索引
    random_indices = random.sample(all_indices, k)
    return random_indices

#测试方法
from sklearn import metrics
def caculateMetrix(groundtruth, predict, baselineName, type='binary'):
    if type == 'binary':
        acc = metrics.accuracy_score(groundtruth, predict)
        precision = metrics.precision_score(groundtruth, predict, zero_division=True )
        recall = metrics.recall_score(groundtruth, predict, zero_division=True)
        f1 = metrics.f1_score(groundtruth, predict, zero_division=True)
        tn, fp, fn, tp = metrics.confusion_matrix(groundtruth, predict).ravel()
        npv = tn/(fn+tn+1.4E-45)
        # print('%12s'%baselineName, '\t\t%.6f' %acc,'\t%.6f'% precision,'\t%.6f'%npv,'\t%.6f'% recall,'\t%.6f'% f1, '\t', 'tp:',tp,'fp:',fp,'fn:',fn,'tn:',tn)
        # print('{:<24} {:<14.6f} {:<25.6f} {:<15.6f} {:<15.6f} {:<20.6f} tp: {} fp: {} fn: {} tn: {}'.format(baselineName, acc, precision, npv, recall, f1, tp, fp, fn, tn))
        print('{:<24} {:<14.6f} {:<25.6f} {:<15.6f} {:<15.6f} {:<20.6f} {:<15} {:<15} {:<15} {}'.format(baselineName, acc, precision, npv, recall, f1, tp, fp, fn, tn))

    if type =='include_unfind':
        evadf = pd.DataFrame()
        evadf['g'] = groundtruth
        evadf['p'] = predict

        evadf_hot = evadf[~evadf.p.isnull()]
        evadf_cold = evadf[evadf.p.isnull()]

        tp = len(evadf_hot[(evadf_hot.g.astype('int')==1) & (evadf_hot.p.astype('int')==1)])
        fp = len(evadf_hot[(evadf_hot.g.astype('int')==0) & (evadf_hot.p.astype('int')==1)])        
        tn = len(evadf_hot[(evadf_hot.g.astype('int')==0) & (evadf_hot.p.astype('int')==0)])
        fn = len(evadf_hot[(evadf_hot.g.astype('int')==1) & (evadf_hot.p.astype('int')==0)])
        up = len(evadf_cold[evadf_cold.g==1])
        un = len(evadf_cold[evadf_cold.g==0])
        acc = (tp+tn)/(tp+fp+tn+fn+up+un)
        precision = tp/(tp+fp)
        npv = tn/(tn+fn)
        recall = tp/(tp+fn+up)
        f1=(2*precision*recall)/(precision+recall)
        print( baselineName, 
                '\t%.6f' %acc,
                '\t%.6f'% precision,
                '\t%.6f'%npv,
                '\t%.6f'% recall,
                '\t%.6f'% f1, '\t', 
                'tp:',tp,'fp:',fp,'fn:',fn,'tn:',tn, 'up:',up, 'un:',un)

    if type == 'multi':
        acc = metrics.accuracy_score(groundtruth, predict)
        precision = metrics.precision_score(groundtruth, predict, average='macro', zero_division=True )
        recall = metrics.recall_score(groundtruth, predict, average='macro', zero_division=True)
        f1 = metrics.f1_score(groundtruth, predict, average='macro', zero_division=True)
        return acc,precision,recall,f1
        print('%12s'%baselineName, ' \t%.6f '%acc,'\t%.6f'% precision, '\t%.6f'% recall,'\t%.6f'% f1)

@measure_execution_time
def val_model(val_loader,model,criterion,device):
    model.eval()
    total_loss = 0.0
    # 创建空张量用于存储模型输出和标签
    all_outputs = torch.tensor([])
    all_labels = torch.tensor([])
    with torch.no_grad():
        for batch_x_squences, batch_y_labels in val_loader:
            batch_x_squences = list(batch_x_squences)
            outputs = model(batch_x_squences)[0]
            batch_y_labels = batch_y_labels.to(device)
            loss = criterion(outputs, batch_y_labels)
            total_loss += loss.item()
            # 将当前批次的模型输出和标签添加到张量中
            all_outputs = torch.cat((all_outputs, outputs.cpu()), dim=0)
            all_labels = torch.cat((all_labels, batch_y_labels.cpu()), dim=0)
    # 计算准确率
    y_true = (all_labels==1)
    y_pred = (torch.sigmoid(all_outputs)>0.5)
    accuracy = caculateMetrix(y_true, y_pred, baselineName="baselineName", type='multi')
    return total_loss/len(val_loader),accuracy

@measure_execution_time
#region train epoch
def train_epoch(model, train_loader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0
    # 创建空张量用于存储模型输出和标签
    all_outputs = torch.tensor([])
    all_labels = torch.tensor([])
    for batch_x_squences, batch_y_labels in train_loader:
        optimizer.zero_grad()
        batch_x_squences = list(batch_x_squences)
        outputs = model(batch_x_squences)[0]
        batch_y_labels = batch_y_labels.to(device)
        loss = criterion(outputs, batch_y_labels)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
        # 将当前批次的模型输出和标签添加到张量中
        all_outputs = torch.cat((all_outputs, outputs.cpu()), dim=0)#移到cpu上减少显存占用增加
        all_labels = torch.cat((all_labels, batch_y_labels.cpu()), dim=0)

    avg_loss = total_loss / len(train_loader)
    # 计算准确率
    y_true = (all_labels==1)
    y_pred = (torch.sigmoid(all_outputs)>0.5)
    accuracy = caculateMetrix(y_true, y_pred, baselineName="baselineName", type='multi')
    return avg_loss,accuracy

#endregion
import subprocess
def print_gpu_memory():
    command = "nvidia-smi --query-gpu=memory.total,memory.used --format=csv,noheader,nounits"
    memory_info = subprocess.check_output(command, shell=True).decode()
    memory_info_list = memory_info.strip().split('\n')

    for i, memory in enumerate(memory_info_list):
        total_memory, used_memory = memory.split(',')
        print(f"GPU {i} - 已用显存: {used_memory.strip()} MB, 总显存: {total_memory.strip()} MB")

from torch.nn.parallel._functions import Scatter
from torch.nn.parallel import DataParallel
import torch
# This code was copied from torch.nn.parallel and adapted for DataParallel to chunk lists instead of duplicating them
# (this is really all this code is here for)
def scatter(inputs, target_gpus, dim=0):
    def scatter_map(obj):
        if isinstance(obj, torch.Tensor):
            return Scatter.apply(target_gpus, None, dim, obj)
        if isinstance(obj, tuple) and len(obj) > 0:
            return list(zip(*map(scatter_map, obj)))
        if isinstance(obj, list) and len(obj) > 0:
            if not isinstance(obj, list) or len(obj) == 0 or len(target_gpus) == 0:
                return []
            num_gpus = len(target_gpus)
            num_samples = len(obj)
            samples_per_gpu = num_samples // num_gpus
            remaining_samples = num_samples % num_gpus

            distributed_samples = []
            start_idx = 0
            for gpu_idx in range(num_gpus):
                gpu_samples_count = samples_per_gpu + (1 if gpu_idx < remaining_samples else 0)
                gpu_samples = obj[start_idx:start_idx + gpu_samples_count]
                distributed_samples.append(gpu_samples)
                start_idx += gpu_samples_count

            return distributed_samples
        if isinstance(obj, dict) and len(obj) > 0:
            return list(map(type(obj), zip(*map(scatter_map, obj.items()))))
        return [obj for _ in target_gpus]

    # After scatter_map is called, a scatter_map cell will exist. This cell
    # has a reference to the actual function scatter_map, which has references
    # to a closure that has a reference to the scatter_map cell (because the
    # fn is recursive). To avoid this reference cycle, we set the function to
    # None, clearing the cell
    try:
        return scatter_map(inputs)
    finally:
        scatter_map = None


def scatter_kwargs(inputs, kwargs, target_gpus, dim=0):
    inputs = scatter(inputs, target_gpus, dim) if inputs else []
    kwargs = scatter(kwargs, target_gpus, dim) if kwargs else []
    if len(inputs) < len(kwargs):
        inputs.extend([() for _ in range(len(kwargs) - len(inputs))])
    elif len(kwargs) < len(inputs):
        kwargs.extend([{} for _ in range(len(inputs) - len(kwargs))])
    inputs = tuple(inputs)
    kwargs = tuple(kwargs)
    return inputs, kwargs


class DataParallelV2(DataParallel):
    def scatter(self, inputs, kwargs, device_ids):
        return scatter_kwargs(inputs, kwargs, device_ids, dim=self.dim)

def train(train_path, test_path,batch_size=8, esm_out_dim=1280, gru_h_dim=256, att_dim=32, dropout_rate=0.2, freeze_esm_layers = 23,epoches=5):
    # loading data
    #Train data   
    ds_train = pd.read_feather(train_path)
    cls = get_cdhit_results(cdhit_clstr_file="/hpcfs/fhome/zhujun1/project/protein/alfp/taskModel240524/temp/cdhit_results_2025-03-22_16_57_54_4SOLE5wlJQfjC6GZ.clstr")
    cluster_size = cls.cluster_id.value_counts()
    cluster_size = pd.DataFrame({'cluster_id':cluster_size.index,'cluster_size':cluster_size.values})
    cls = cls.merge(cluster_size, on='cluster_id', how='left')
    merged_data = pd.merge(cls[['cluster_id', 'uniprot_id','is_representative']], ds_train, on='uniprot_id')
    train_data = merged_data[['cluster_id', 'uniprot_id', 'is_representative','seq','label']]
    train_data = train_data[train_data['is_representative'] == True]
    X = train_data["seq"].to_numpy()
    Y = np.vstack(train_data["label"])
    # 训练集数据处理
    print("训练集",X.shape)
    X,Y = data_process_drop(X,Y)

    #Test data  
    test_data = pd.read_feather(test_path)
    test_sequences = test_data["seq"].to_numpy()
    test_labels = np.stack(test_data["label"])
    print("测试集",test_labels.shape)
    # 验证集数据处理，过长的蛋白质序列
    test_sequences,test_labels = data_process_drop(test_sequences,test_labels)
    
    #构建测试集dataloader
    test_dataset = ProteinDataset(test_sequences, test_labels)
    test_loader = DataLoader(test_dataset, batch_size = batch_size, shuffle=False)
    # 从训练集抽取一定比例的数据作为初始化集，其余作为样本池
    train_sequences, val_sequences, train_labels, val_labels = train_test_split(X, Y, test_size = 0.2, random_state=42)
    #构建训练集loader
    train_dataset = ProteinDataset(train_sequences, train_labels)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True,drop_last=True)
    labeled_number = len(train_dataset)
    print("loader_finish训练集大小==============",labeled_number)
    #构建验证集dataloader
    val_dataset = ProteinDataset(val_sequences, val_labels)
    val_loader = DataLoader(val_dataset, batch_size = 2*batch_size, shuffle=False)
    freeze_esm_layers_list = [23,28,31,32,33]
    for freeze_esm_layers in freeze_esm_layers_list:


        output_dimensions = Y.shape[1]  # Number of output classes
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = BGRU(
            input_dimensions=esm_out_dim,
            device = device,
            gru_h_size=gru_h_dim,
            attention_size=att_dim,
            dropout=dropout_rate,
            output_dimensions=output_dimensions,
            freeze_esm_layers=freeze_esm_layers,
        )
        model.to(device)
        # 获取可用的GPU数量list(range(args.devices))
        device_ids = list(range(args.devices))
        model = DataParallelV2(model, device_ids=device_ids)
        criterion = nn.BCEWithLogitsLoss()
        optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=0.0001)
        logname = args.taskname + str(freeze_esm_layers) + ".txt"
        # 以追加模式打开文件并写入数据
        with open(logname, 'a') as f:
            f.write(f"冻结层数为：{freeze_esm_layers}\n")
        cycle_loss=0.0
        cycle_acc=0.0
        best_accuracy=0.0
        for epoch in range(20):
            cycle_loss ,cycle_acc = train_epoch(model=model, train_loader=train_loader, optimizer=optimizer, criterion=criterion, device= device)
            # 写入日志文件
            with open(logname, "a") as file:
                file.write(f"Epoches [{epoch + 1}/{epoches}] - Train Loss: {cycle_loss:.6f} - Train Acc: {cycle_acc}\n")
            print("开始测试")
            val_loss,val_accuracy = val_model(val_loader, model,criterion,device)
            test_loss,test_accuracy = val_model(test_loader, model,criterion,device)
            print(f"Epoch [{freeze_esm_layers+1}/{33}] - Train Avg Loss: {cycle_loss:.6f} - Train Avg Acc: {cycle_acc}- Val Avg Loss: {val_loss:.6f} - Val Avg Acc: {val_accuracy}- Test Avg Loss: {test_loss:.6f} - Test Avg Acc: {test_accuracy}")
            # 写入日志文件
            with open(logname, "a") as file:
                file.write(f"Epoch [{freeze_esm_layers+1}/{33}] - Train Avg Loss: {cycle_loss:.6f} - Train Avg Acc: {cycle_acc}- Val Avg Loss: {val_loss:.6f} - Val Avg Acc: {val_accuracy}- Test Avg Loss: {test_loss:.6f} - Test Avg Acc: {test_accuracy}\n")
            #保存最佳模型
            best_mode_pth = args.taskname + 'mode.pth'
            if best_accuracy<test_accuracy[0]:
                torch.save(model.state_dict(),best_mode_pth)
                best_accuracy = test_accuracy[0]


#丢弃超过1022的样本
def data_process_drop(X,Y):
    X_processed = []  # 用于存储处理后的数据
    Y_processed = []  # 用于存储处理后的标签
    for i, x in enumerate(X):
        if len(x) < 1022: 
            X_processed.append(x)
            Y_processed.append(Y[i])  # 原始数据不需要拆分，直接添加到处理后的数据中
    X_processed = np.array(X_processed)
    Y_processed = np.vstack(Y_processed)
    return X_processed, Y_processed

def get_cdhit_results(cdhit_clstr_file):
    """读取cdhit聚类结果
    Args:
        cdhit_clstr_file (string): 聚类结果文件
    Returns:
        DataFrame: ['cluster_id','uniprot_id','identity'， 'is_representative']
    """
    counter = 0
    res = []
    with open(cdhit_clstr_file,'r') as f:
        for line in f:
            if 'Cluster' in line:
                cluster_id = line.replace('>Cluster','').replace('\n', '').strip()
                continue
            str_uids= line.replace('\n','').split('>')[1].replace('at ','').split('... ')
                        
            if '*' in str_uids[1]:
                identity = 1
                isrep = True
            else:
                identity = float(str_uids[1].strip('%')) /100
                isrep = False

            res = res +[[cluster_id, str_uids[0], identity, isrep ]]

    resdf = pd.DataFrame(res, columns=['cluster_id','uniprot_id','identity', 'is_representative']) #转换为DataFrame
    return resdf

def step_by_step_run():
    # loading data path
    train_path = "task240524/ds_train.feather"
    test_path = "task240524/ds_test.feather"
    # train中长度超过1022的字符串个数
    # test中长度超过1022的字符串个数
    print("----------starting-------")
    train(train_path=train_path,test_path=test_path,batch_size=args.batchsize, esm_out_dim=1280, gru_h_dim=512, att_dim=32, dropout_rate=0.2, freeze_esm_layers=32,epoches=5)

import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--taskname', type=str, default = "Freese_Test")
    parser.add_argument('--batchsize', type=int, default=64)
    parser.add_argument('--devices', type=int, default=8)
    args = parser.parse_args()
    step_by_step_run()
    print("success")