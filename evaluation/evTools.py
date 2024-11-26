# Standard Library Imports
import os
import sys

# Third-party Imports
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn import metrics
from concurrent.futures import ThreadPoolExecutor, as_completed
import plotly.graph_objects as go
from IPython.display import HTML
from pandarallel import pandarallel  # Importing pandarallel for parallel processing

# Setting up the path for the module
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1, '../')

# Local Imports
from config import conf as cfg
from tools import btools








# Read CSV files serially
def read_10fold_res_csv_files(file_paths):
    return [pd.read_csv(file, sep='\t') for file in file_paths]

def process_no_res(res_list, eckey, rxnkey):
    # 遍历 res_list 中的每个数据帧，计算总条目数、无预测条目数和无反应 EC 条目数
    summary = [
        [
            len(pred_detail),
            len(pred_detail[pred_detail[eckey].str.contains('NO-PREDICTION')]),
            len(pred_detail[pred_detail[rxnkey].str.contains('EC-WITHOUT-REACTION')])
        ]
        for pred_detail in res_list
    ]
    
    # 创建 DataFrame 返回结果
    return pd.DataFrame(summary, columns=['test_size', 'no_prediction', 'ec_without_rxn'])


# 并行执行标签创建的函数
def make_10folds_labels(resdf, columns_dict, rxn_label_dict, fold_num=10):
    res = []
    for i in tqdm(range(fold_num)):
        
        for src_col, lb in columns_dict.items():
            resdf[i][lb] = resdf[i][src_col].apply(lambda reaction_id: btools.make_label(reaction_id=str(reaction_id), rxn_label_dict=rxn_label_dict))
        resdf[i]['run_fold'] = i+1
        # res = res +  resdf[i]
    resdf = pd.concat(resdf, axis=0).reset_index(drop=True)
    return resdf
    


# Function to calculate metrics
def calculate_metrics(eva_df, ground_truth_col, pred_col, eva_name, avg_method='weighted'):
    res =  btools.rxn_eva_metric_with_colName(eva_df=eva_df, col_groundtruth=ground_truth_col, col_pred=pred_col, eva_name=eva_name, average_type=avg_method)
    return res

# 多线程运行评价函数
def calculate_metrics_parallel(res_df, ground_truth_col, pred_col, avg_method='weighted', max_workers=None):
    def run_metric_evaluation(index):
        return calculate_metrics(eva_df=res_df[index], ground_truth_col=ground_truth_col, pred_col=pred_col, eva_name=f'fold{index + 1}', avg_method=avg_method)
    
    results = [None] * len(res_df)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(run_metric_evaluation, i): i
            for i in range(len(res_df))
        }
        for future in as_completed(futures):
            i = futures[future]
            results[i] = future.result()
            
    results = pd.concat(results,axis=0).reset_index(drop=True)
    
    return results


#region 展开获取均值方差
def get_fold_mean_std_metrics(input_df):
    # 对数值列进行分组聚合
    res_fold_std = input_df[['baselineName', 'avgType','mAccuracy', 'mPrecision', 'mRecall', 'mF1']]
    # 对每个 baselineName 进行分组并计算均值和标准差
    res_fold_std = res_fold_std.groupby(['baselineName','avgType']).agg(['mean', 'std'])
    # 重置索引以简化处理
    res_fold_std = res_fold_std.reset_index()
    
    
    # 修改列名，将 MultiIndex 转换为单层列名
    res_fold_std.columns = ['_'.join(filter(None, col)).strip() for col in res_fold_std.columns]
    # 使用 melt 方法将列转换为行
    res_fold_std_melted = res_fold_std.melt(id_vars=['baselineName', 'avgType'], var_name='Metric_Statistic', value_name='Value')
    # 将 'Metric_Statistic' 列分割成 'Metric' 和 'Statistic'
    res_fold_std_melted[['Metric', 'Statistic']] = res_fold_std_melted['Metric_Statistic'].str.rsplit('_', n=1, expand=True)
    
    res_fold_std_melted = res_fold_std_melted.sort_values(by=['baselineName',  'Metric', 'avgType']).reset_index(drop=True)
    res_fold_std_melted = res_fold_std_melted.drop(columns=['Metric_Statistic'])
    
    # 使用 pivot 将数据转换为所需格式
    res_fold_std_pivot = res_fold_std_melted.pivot_table(index=['baselineName','avgType', 'Metric'], columns='Statistic', values='Value').reset_index()
    res_fold_std_pivot.columns.name = None
    return res_fold_std_pivot
#endregion

#region 统计没有预测结果和有EC没反应的结果
def statistic_no_res(res_df, name_col_ec, name_col_rxn):
    grouped_counts = res_df.groupby('run_fold').agg(
    test_size=('run_fold', 'count'),
    no_prediction_count=(name_col_ec, lambda x: (x == 'NO-PREDICTION').sum()),
    ec_without_reaction_count=(name_col_rxn, lambda x: (x == 'EC-WITHOUT-REACTION').sum())
    ).reset_index()
    
    return grouped_counts
#endregion



# region  计算10折交叉验证的结果评价指标
def get_eval_results_ec(baselineName, dict_rxn2id):
    
    label_file = f'{cfg.DIR_PROJECT_ROOT}/results/intermediate/ecmethods/{baselineName}_10folds_labels_res.feather'
    
    #  处理lable 数据
    if not os.path.exists(label_file):
        vali_res_file_path = [f'{cfg.DIR_RES_BASELINE}results/ec_methods/{baselineName}/fold{item}.tsv' for item in range(1, 11)]
        res =read_10fold_res_csv_files(vali_res_file_path)
        
        print('Labeling ...')
        columns_dict = {
            'rxn_groundtruth': 'lb_rxn_groundtruth',
            f'rxn_{baselineName}': f'lb_rxn_{baselineName}'
        }
        
        res = make_10folds_labels(resdf=res, columns_dict=columns_dict, rxn_label_dict=dict_rxn2id, fold_num=10)
        res.to_feather(label_file)
    
    # 计算指标
    print('Calculating metrics ...')
    metrics_flie = f'{cfg.DIR_PROJECT_ROOT}/results/intermediate/ecmethods/{baselineName}_10_folds_metrics_res.feather'
    
    if not os.path.exists(metrics_flie):
        res = pd.read_feather(label_file)
        metrics = eva_cross_validation(res_df=res, lb_groundtruth='lb_rxn_groundtruth', lb_predict=f'lb_rxn_{baselineName}', num_folds=10)
        metrics.insert(0,'baselineName', f'{baselineName}')
        metrics.to_feather(metrics_flie)

    # 计算均值方差
    print('Calculating mean and std ...')
    std_file = f'{cfg.DIR_PROJECT_ROOT}/results/intermediate/ecmethods/{baselineName}_10_folds_std_res.feather'
    if not os.path.exists(std_file):
        metrics = pd.read_feather(metrics_flie)
        std = get_fold_mean_std_metrics(input_df= metrics)
        std.to_feather(std_file)
        
    
    # 统计没有预测的rxn
    print('Statistic ec no prediction and ec with no reaction    ...')
    ec_no_rxn_file = f'{cfg.DIR_PROJECT_ROOT}/results/intermediate/ecmethods/{baselineName}_10_folds_no_rxn_res.feather'
    if not os.path.exists(ec_no_rxn_file):
        res = pd.read_feather(label_file)
        no_rxn_10_fold = statistic_no_res(res_df=res, name_col_ec=f'ec_{baselineName}', name_col_rxn=f'rxn_{baselineName}')
        no_rxn_10_fold.to_feather(ec_no_rxn_file)
        
      
    std = pd.read_feather(std_file)     
    metrics = pd.read_feather(metrics_flie)
    no_rxn_10_fold = pd.read_feather(ec_no_rxn_file)
        
    return std, metrics, no_rxn_10_fold
#endregion


#region Jpuyter 显示 10折交叉验证的结果HTML
# Function to display results as HTML
def display_html_results(metrics, std_mean, no_pred, eva_name):
    return HTML(f"""
         <div style="float:left; width:700px;">
              <h2 style='color:blue'>{eva_name} Evaluation 10 Fold Details</h2>
              {metrics.to_html()}
         </div>
         <div  style="float:left; width:600px;" >
              <h2 style='color:blue' >{eva_name} Evaluation 10 Fold std and mean Overview</h2>
                   {std_mean.to_html()}
         </div>
         <div style="float:left; width:600px;">
         <h2 style='color:blue' >{eva_name} Evaluation 10 Fold No Prediction Overview</h2>
              {no_pred.to_html()}
         </div>
         """)
#endregion


def show_ec_methods_10_eva_fig(res_metrics_data):

    # 定义颜色
    colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC', '#E7DAD2', '#999999']

    # 确保 res_fold_std 按照 F1 值排序，但保留所有度量以便绘图
    # 获取唯一的方法名称
    methods = res_metrics_data[res_metrics_data.Metric == 'mF1'].sort_values(by=['mean']).baselineName.tolist()
    
    # 创建一个绘图对象
    fig = go.Figure()

    # 为每个方法添加柱状图，并使用指定颜色
    for idx, method in enumerate(methods):
        df_method = res_metrics_data[res_metrics_data['baselineName'] == method]
        # 只在 df_method 非空的情况下添加 trace
        
        showMetics = ['mAccuracy', 'mPrecision', 'mRecall', 'mF1']
        
        if not df_method.empty:
            fig.add_trace(
                go.Bar(
                    name=method,
                    x=showMetics,
                    y=df_method['mean'],
                    width=0.11,  # 调整柱状图的宽度
                    error_y=dict(
                        type='data',
                        array=df_method['std'],
                        visible=True
                    ),
                    marker_color=colors[idx % len(colors)],  # 使用指定的颜色循环
                    marker_line=dict(color='black', width=0.5)  # 添加黑色1px外边框
                )
            )

    # 更新布局
    fig.update_layout(
        yaxis=dict(
            showline=True, linecolor='black', linewidth=1,
            showgrid=True, gridcolor='gray', gridwidth=1,
            minor=dict(
                showgrid=True,
                griddash='dash',  # 设置网格线为虚线
                gridcolor='gray',
                gridwidth=0.5,
                dtick=0.05
            ),
            dtick=0.1,  # 每隔 0.1 的主要网格线
            range=[0, 1.1]  # y 轴的刻度最大值调整为 1.2
        ),
        xaxis=dict(showline=True,linecolor='#000000',linewidth=1),       
        title=dict(
        text='Performance Comparison of EC-based Reaction Prediction Methods',
        x=0.5,
        y=0.05,
        xanchor='center',
        font=dict(
            size=20,
            weight=1000
        )
        ),
        # xaxis_title='Metric',
        yaxis_title='Mean Value',
        width=1600,
        height=600,
        barmode='group',
        bargap=0.29,  # 设置柱状图之间的间距
        template='plotly_white',
        legend=dict(
            orientation="h",  # 水平放置图例
            yanchor="bottom",
            y=0.99,  # 将图例放置在下方
            xanchor="center",
            x=0.5,
            tracegroupgap=2,  # 图例列数
            borderwidth=1,  # 添加黑色1px外边框
        ),

    )

    return fig
























def get_simi_Pred(pred_list, uniprot_rxn_dict, topk=3):
    uniprot_id_list = [item[0] for item in pred_list][:topk]
    rxn_ids = [uniprot_rxn_dict.get(uniprot_id) for uniprot_id in uniprot_id_list]
    rxn_res = (cfg.SPLITER).join(set(rxn_ids))
    return rxn_res


    
def read_h5_file(file_path):
    with pd.HDFStore(file_path, 'r') as h5:
        data = h5['data']
    return data





# 计算评估指标的函数，保持并行化以提升速度
def calculate_metrics_multi_joblib(groundtruth, predict, average_type, print_flag=False, n_jobs=4):
    metric_functions = [
        lambda gt, pr: metrics.accuracy_score(gt, pr),
        lambda gt, pr: metrics.precision_score(gt, pr, average=average_type, zero_division=True),
        lambda gt, pr: metrics.recall_score(gt, pr, average=average_type, zero_division=True),
        lambda gt, pr: metrics.f1_score(gt, pr, average=average_type, zero_division=True)
    ]

    # 使用 joblib 并行化计算
    results = Parallel(n_jobs=n_jobs)(
        delayed(metric_fn)(groundtruth, predict) for metric_fn in metric_functions
    )

    if print_flag:
        print(f'{results[0]:.6f}\t{results[1]:.6f}\t{results[2]:.6f}\t{results[3]:.6f}\t{average_type:>12s}')

    return results + [average_type]

# 评估单个折的函数
def eva_one_fold(eva_df, lb_groundtruth, lb_predict, fold_num=None, n_jobs=4):
    # 提取 groundtruth 和 predict 数据
    groundtruth = np.stack(eva_df[lb_groundtruth])
    predict = np.stack(eva_df[lb_predict])
    
    # 定义需要计算的平均类型
    average_types = ['weighted', 'micro', 'macro', 'samples']
    
    # 并行计算不同的平均类型
    results = Parallel(n_jobs=n_jobs)(
        delayed(calculate_metrics_multi_joblib)(
            groundtruth=groundtruth,
            predict=predict,
            average_type=avg_type,
            print_flag=False,
            n_jobs=1  # 内部不再嵌套并行，避免资源争用
        ) for avg_type in average_types
    )
    
    # 处理结果并创建 DataFrame
    res = pd.DataFrame(results, columns=['mAccuracy', 'mPrecision', 'mRecall', 'mF1', 'avgType'])
    if fold_num is not None:
        res.insert(0, 'runFold', fold_num)
    
    return res

# 执行多折交叉验证
def eva_cross_validation(res_df, lb_groundtruth, lb_predict, num_folds=10):

    eva_metrics = []
    for runfold in tqdm(range(1, num_folds+1)):
        res = eva_one_fold(eva_df=res_df[res_df.run_fold==runfold].reset_index(drop=True).head(100), lb_groundtruth=lb_groundtruth, lb_predict=lb_predict,fold_num=runfold)
        eva_metrics = eva_metrics + [res]
        
    eva_metrics = pd.concat(eva_metrics, axis=0).reset_index(drop=True)
    
    return eva_metrics