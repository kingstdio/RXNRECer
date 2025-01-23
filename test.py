import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../../')
sys.path.insert(1,'./')
from config import conf as cfg
import rxnrecer as production
from modules import commonfunction as cmfunc
from tqdm import tqdm
from pandarallel import pandarallel # 导入pandaralle
pandarallel.initialize(progress_bar=False)
import json
import numpy as np
from modules.embedding import seqEmbedding as ebdseq
from modules.embedding import t5Embedding as ebdt5
import pandas as pd



if __name__ == '__main__':

    input_protein_df = pd.read_feather(cfg.FILE_DS_TEST)[['uniprot_id', 'seq']]
    res  = production.multi_batch_run_prediction(input_protein_df.iloc[0:20], getEquation=False, Ensemble=True, batch_size=10)
    res