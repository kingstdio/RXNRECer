import sys, os
import pandas as pd
import numpy as np
import argparse
import concurrent.futures
import warnings
import matplotlib.pyplot as plt
from tqdm import tqdm

# Adding the project directory to system path for imports
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1, '../')

from config import conf as cfg
import modules.simi_caculator as simitool

# Dictionary mapping embedding methods to feature file paths
dict_feature_path = {
    'esm': cfg.FILE_EMBD_PROTEIN_ESM2_L33_650M,
    'unirep': cfg.FILE_EMBD_PROTEIN_UNIREP,
    't5': cfg.FILE_EMBD_PROTEIN_T5_SEQ
}

# Function to load features from a given method
def load_feature(method):
    print(f'Loading {method} features from {dict_feature_path.get(method)}')
    feature = pd.read_feather(dict_feature_path.get(method))
    feature = feature.rename(columns={
        'rep33': 'features' if method == 'esm' else None,
        'unirep': 'features' if method == 'unirep' else None,
        't5_per_protein': 'features' if method == 't5' else None
    }).dropna(axis=1, how='all')
    print(f'{method} features loaded, shape: {feature.shape}')
    return feature[['uniprot_id', 'features']]

# Function to merge features with training and testing data
def merge_features(embdfeatures, data_train, data_test):
    print('Adding features to training and test data')
    res_train = data_train.merge(embdfeatures, on='uniprot_id', how='left')
    res_test = data_test.merge(embdfeatures, on='uniprot_id', how='left')
    return res_train, res_test

# Function to load dataset from the given path
def load_data(path):
    data = pd.read_feather(path)
    return data[['uniprot_id', 'reaction_id', 'ec_number']]

# Function to load training and testing datasets for a specific fold
def load_dataset(fold_num):
    train_path = f'{cfg.DIR_DATASET}validation/fold{fold_num}/train.feather'
    test_path = f'{cfg.DIR_DATASET}validation/fold{fold_num}/valid.feather'
    print(f'Loading data from {train_path}')
    data_train = load_data(train_path)
    print(f'Loading data from {test_path}')
    data_test = load_data(test_path)
    print(f'Data loaded, Fold: {fold_num}. Train size: {len(data_train)}, Test size: {len(data_test)}')
    return data_train, data_test

# Function to calculate similarity between protein features
def get_top_protein_simi(x_feature, y_feature, y_uniprot_id, topk):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_euclidean = executor.submit(simitool.get_euclidean_distances, x_feature, y_feature)
        future_cosine = executor.submit(simitool.get_cosine_similarity, x_feature, y_feature)
        res_euclidean = np.round(future_euclidean.result()[0], 6)
        res_cosine = np.round(future_cosine.result()[0], 6)
    res_euclidean_dict = dict(zip(y_uniprot_id, res_euclidean))
    final_res_euclidean = sorted(res_euclidean_dict.items(), key=lambda item: item[1])[:topk]
    res_cosine_dict = dict(zip(y_uniprot_id, res_cosine))
    final_res_cosine = sorted(res_cosine_dict.items(), key=lambda item: item[1], reverse=True)[:topk]
    return final_res_euclidean, final_res_cosine

# Function to run evaluation on a specific fold
def run_fold(fold_num: int, embd_methd: str):
    print(f'Running fold {fold_num} with embedding method {embd_methd}...')

    # Load training and test datasets
    data_train, data_test = load_dataset(fold_num=fold_num)

    # Load feature vectors
    embd_featuers = load_feature(method=embd_methd)

    # Merge features with training and test data
    data_train, data_test = merge_features(embdfeatures=embd_featuers, data_train=data_train, data_test=data_test)

    # Prepare bank features for similarity calculations
    dim_features = len(data_train.features.values[0])
    bank_features = np.concatenate(data_train.features.values).reshape(-1, dim_features)

    # Testing mode: limit the number of test data points
    # data_test = data_test.head(10)

    print('Calculating similarity between test set and bank set...')
    tqdm.pandas()
    data_test[['euclidean', 'cosine']] = data_test.progress_apply(
        lambda x: get_top_protein_simi(
            x_feature=x.features[np.newaxis, :],
            y_feature=bank_features,
            y_uniprot_id=data_train.uniprot_id.to_list(),
            topk=100
        ),
        axis=1, result_type="expand"
    )

    data_test = data_test[['uniprot_id', 'reaction_id', 'ec_number', 'euclidean', 'cosine']]

    # Save the results to HDF5 file
    try:
        warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
        with pd.HDFStore(f'{cfg.RESULTS_DIR}simi/fold_{fold_num}_{embd_methd}_results.h5', 'w') as h5:
            h5['data'] = data_test
    except Exception as e:
        print(f"Error writing HDF5 file: {e}")
    
    print(f'Similarity calculation for fold {fold_num} and embedding method {embd_methd} is done, file saved to {cfg.RESULTS_DIR}simi/fold_{fold_num}_{embd_methd}_results.h5')

# Entry point of the script
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Blast and DeepEC Evaluation")
    parser.add_argument('--foldnum', type=int, default=1, help='Fold number to run')
    parser.add_argument('--embdmethod', type=str, default='esm', choices=['esm', 'unirep', 't5'], help='Embedding method to use')
    args = parser.parse_args()

    run_fold(fold_num=args.foldnum, embd_methd=args.embdmethod)
