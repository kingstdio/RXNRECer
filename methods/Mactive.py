from cProfile import label
import os
import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.nn.parallel
from torch.nn.parallel._functions import Scatter
from torch.nn.parallel import DataParallel
import torch
from tqdm import tqdm
import torch.backends.cudnn as cudnn
import pandas as pd
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
import esm
import argparse
from sklearn.model_selection import train_test_split
from tkinter import _flatten  # type: ignore
from collections import OrderedDict
from sklearn import metrics
from sklearn.metrics import pairwise_distances

import json

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
    
# region Step 2: Define the EsmEmbedding with layer freezing
class EsmEmbedding(nn.Module):
    def __init__(self, device, freeze_esm_layers):
        super(EsmEmbedding, self).__init__()
        self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.freeze_esm_layers = freeze_esm_layers
        self.device = device

        # Freeze the specified layers
        for name, param in self.esm_model.named_parameters():
            if any(layer_name in name for layer_name in [f"layer.{i}." for i in range(23)]):
                param.requires_grad = False
        self.esm_model = self.esm_model.to(self.device)

    def forward(self, x):
        batch_converter = self.esm_alphabet.get_batch_converter()
        #破碎，把长度超过1022的截断，随机裁剪成1022长度的字符串
        # x = breakSquences(x)
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
        output1 = self.output_layer(gru_output.squeeze(1))
        return output,output1


#region  predict reactions with trained model using protein sequences
def predict_sequences(model, sequences, model_weight_path, dict_path, batch_size=512, device=torch.device("cpu"), seq_cut_threshold=10000):
    
    data_loader = DataLoader(sequences, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=True)
           
    # 加载已保存模型
    state_dict = torch.load(model_weight_path, map_location=device, weights_only=True)
    model.load_state_dict({k.replace("module.", ""): v for k, v in state_dict.items()})
    model.to(device)

    # 加载字典
    with open(dict_path, 'r') as f:
        rhea_dict = json.load(f)

    # 获取模型预测结果
    model.eval()
    all_reactions = []
    
    predictions = []
    with torch.no_grad():
        for batch_x in tqdm(data_loader, desc="Predicting reactions"):
            
            batch_x = [seq[:seq_cut_threshold] for seq in batch_x]
            
            outputs = model(batch_x)[0]
            sigmoid_outputs = torch.sigmoid(outputs)
            y_preds = (sigmoid_outputs > 0.5).int().cpu().tolist()

            # 更新标签：确保每个样本至少有一个标签
            for j, pred in enumerate(y_preds):
                if sum(pred) == 0:  # 如果没有大于0.5的标签
                    max_index = torch.argmax(sigmoid_outputs[j]).item()  # 找到最大值的索引
                    y_preds[j][max_index] = 1  # 将对应的标签设为1
                    
            reactions = [
                ';'.join(rhea_dict[str(index)] for index, value in enumerate(pred) if value != 0)
                for pred in y_preds
            ]
                     
            all_reactions.extend(reactions)
            
        return all_reactions
#endregion



#采样
def select_samples(train_sequences,train_labels,model,k,batch_size):
    if train_sequences.shape[0]<k:
        return list(range(train_sequences.shape[0]))
    train_dataset = ProteinDataset(train_sequences, train_labels)
    loader = DataLoader(train_dataset, batch_size=2*batch_size, shuffle=False)
    sums = []
    model.eval()
    with torch.no_grad():
        for batch_x_squences, batch_y_labels in loader:
            outputs = model(batch_x_squences)
            row_sum = torch.sum(torch.abs(outputs[0] - outputs[1]), dim=1)
            sums.append(row_sum)
    # 将列表中的张量拼接成一个张量
    result = torch.cat(sums, dim=0)
    # 使用 torch.topk 函数获取最大值及其索引
    topk_values, topk_indices = torch.topk(result, k)
    return topk_indices

#测试方法

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

#多gpu计算，数据分发方法重写

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
    train_data = pd.read_feather(train_path)
    X = train_data[train_data['isenzyme']!=False]["seq"].to_numpy()
    Y = np.vstack(train_data[train_data['isenzyme']!=False]["label"])
    # 训练集数据处理
    print("训练集",X.shape)
    X,Y = data_process_drop(X,Y)
    
    #Test data
    test_data = pd.read_feather(test_path)
    test_sequences = test_data[test_data['isenzyme']!=False]["seq"].to_numpy()
    test_labels = np.stack(test_data[test_data['isenzyme']!=False]["label"])
    print("测试集",test_sequences.shape)
    # 验证集数据处理，过长的蛋白质序列
    test_sequences,test_labels = data_process_drop(test_sequences,test_labels)
    
    #构建测试集dataloader
    test_dataset = ProteinDataset(test_sequences, test_labels)
    test_loader = DataLoader(test_dataset, batch_size = batch_size, shuffle=False)
    # 从训练集抽取一定比例的数据作为初始化集，其余作为样本池
    train_sequences, labeled_sequences, train_labels, labeled_labels = train_test_split(X, Y, test_size = 0.3, random_state=42)

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
    # 获取可用的GPU数量
    device_ids = list(range(args.devices))
    model = DataParallelV2(model, device_ids=device_ids)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=0.0001)
    logname = args.taskname+".txt"

    Cycles = 30
    poolSize = 200000
    best_accuracy=0.0
    for cycle in range(Cycles):
        #构建训练集loader
        train_dataset = ProteinDataset(labeled_sequences, labeled_labels)
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True,drop_last=True)
        labeled_number = len(train_dataset)
        print("loader_finish训练集大小==============",labeled_number)
        cycle_loss=0.0
        cycle_acc=0.0
        for epoch in range(epoches):
            cycle_loss ,cycle_acc = train_epoch(model=model, train_loader=train_loader, optimizer=optimizer, criterion=criterion, device= device)
            # 写入日志文件
            with open(logname, "a") as file:
                file.write(f"Epoches [{epoch + 1}/{epoches}] - Train Loss: {cycle_loss:.6f} - Train Acc: {cycle_acc}\n")


        test_loss,test_accuracy = val_model(test_loader, model,criterion,device)
        print(f"Epoch [{cycle+1}/{Cycles}] - Train Avg Loss: {cycle_loss:.6f} - Train Avg Acc: {cycle_acc}- Test Avg Loss: {test_loss:.6f} - Test Avg Acc: {test_accuracy}")

        # 写入日志文件
        with open(logname, "a") as file:
            file.write(f"Cycles [{cycle + 1}/{Cycles}] - Sample_Number: {labeled_number} - Train Loss: {cycle_loss:.6f} - Train Acc: {cycle_acc} - Test Loss: {test_loss:.6f} -Test Acc:{test_accuracy}\n")
        #保存最佳模型
        best_mode_pth = args.taskname + 'best_mode.pth'
        if best_accuracy<test_accuracy[0]:
            torch.save(model.state_dict(),best_mode_pth)
            best_accuracy = test_accuracy[0]
        #为了加快采样，将从样本池中随机抽出一半组成未标记样本池,（仅在全部训练数据时使用，单独酶蛋白仅22万数据直接对全部数据采样）
        indices = np.random.permutation(train_sequences.shape[0])
        train_sequences = train_sequences[indices]
        train_labels = train_labels[indices]
        # 选择前 poolSize 个样本
        X_pool = train_sequences
        Y_pool = train_labels
        #采样,全量
        index = select_samples(X_pool, Y_pool, model,k=20000,batch_size=batch_size)
        
        selected_sequences = np.array([X_pool[i] for i in index])
        selected_labels = np.array([Y_pool[i] for i in index])
        train_sequences = np.array([seq for i, seq in enumerate(train_sequences) if i not in index])
        train_labels = np.array([label for i, label in enumerate(train_labels) if i not in index])
        labeled_sequences = np.concatenate((labeled_sequences, selected_sequences))
        labeled_labels = np.concatenate((labeled_labels, selected_labels))

#丢弃超过1022的样本
def data_process_drop(X,Y):
    X_processed = []  # 用于存储处理后的数据
    Y_processed = []  # 用于存储处理后的标签
    for i, x in enumerate(X):
        if len(x) < 1022: 
            X_processed.append(x)
            Y_processed.append(Y[i])  
    X_processed = np.array(X_processed)
    Y_processed = np.vstack(Y_processed)
    return X_processed, Y_processed

def step_by_step_run():
    # loading data path
    train_path = "/hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/datasets/task240524/ds_train.feather"
    test_path = "/hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/datasets/task240524/ds_test.feather"
    # train中长度超过1022的字符串个数
    # test中长度超过1022的字符串个数
    train(train_path=train_path,test_path=test_path,batch_size=args.batchsize, esm_out_dim=1280, gru_h_dim=512, att_dim=32, dropout_rate=0.2, freeze_esm_layers=32,epoches=10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--taskname', type=str, default = "test")
    parser.add_argument('--batchsize', type=int, default=8)
    parser.add_argument('--devices', type=int, default=1)
    args = parser.parse_args()
    step_by_step_run()
    print("success")
