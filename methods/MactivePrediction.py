from math import cos
import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../../')
from config import conf as cfg

import Mactive
import torch
import pandas as pd
from types import SimpleNamespace


def load_model():
    
    mcfg = SimpleNamespace(
    #模型参数
    batch_size=3,
    esm_out_dim=1280,
    gru_h_dim=512,
    att_dim=32,
    dropout_rate=0.2,
    freeze_esm_layers = 32, #冻结层数,
    output_dimensions=10479,
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    model_weight_path=cfg.FILE_WEIGHT_PRODUCTION_BEST_MODEL,
    dict_path = cfg.FILE_DS_DICT_ID2RXN
    )


    model = Mactive.BGRU(
        input_dimensions=mcfg.esm_out_dim,
        device = mcfg.device,
        gru_h_size=mcfg.gru_h_dim,
        attention_size=mcfg.att_dim,
        dropout=mcfg.dropout_rate,
        output_dimensions=mcfg.output_dimensions,
        freeze_esm_layers=mcfg.freeze_esm_layers,
    )
    
    
    
    return model, mcfg

if __name__ == '__main__':
    
    print('Start active learning prediction')
    
    # load datasets
    # ds_test = pd.read_feather(cfg.FILE_DS_TEST)
    
    # ds_train = pd.read_feather(cfg.FILE_DS_TRAIN)
    # ds_train['seq_len'] = ds_train['seq'].apply(lambda x: len(x))
    # ds_trian_10000aa = ds_train[(ds_train.seq_len>10000)&(ds_train.seq_len<120000)].reset_index(drop=True)
    
        
    print('loading model')
    model, mcfg = load_model() 
    print(f'Predction parameters: {mcfg}')


    for i in range(4,11):
        print(f'Loading dataset:{cfg.DIR_DATASET}validation/fold{i}/valid.feather')
        ds_test = pd.read_feather(f'{cfg.DIR_DATASET}validation/fold{i}/valid.feather').reset_index(drop=True)
        sequences = ds_test['seq']
        
        
        print('loading model')
        model, mcfg = load_model() 
        print(f'Predction parameters: {mcfg}')
        print('Predicting reactions')
        consine_rep = Mactive.predict_sequences(model=model, 
                                                sequences=ds_test.seq, 
                                                model_weight_path=mcfg.model_weight_path, 
                                                dict_path=mcfg.dict_path, 
                                                batch_size=mcfg.batch_size,
                                                device=mcfg.device)
        ds_test['alfp_active_prediction_full'] = consine_rep
        ds_test.to_feather(f'alfp_active_prediction_fold{i}_valid.feather')
        
        print(f'Done save to alfp_active_prediction_fold{i}_valid.feather')