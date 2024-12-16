# Standard Library Imports
import os
import sys

# Third-party Imports
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import plotly.graph_objects as go
from IPython.display import HTML
from pandarallel import pandarallel  # Importing pandarallel for parallel processing

# Setting up the path for the module
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1, '../')

# Local Imports
from config import conf as cfg
from tools import btools
import evTools
# Initialize parallel processing
pandarallel.initialize(progress_bar=False)




# 从 JSON 文件加载反应编码字典
with open(cfg.FILE_DS_DICT_RXN2ID, "r") as json_file:
    dict_rxn2id = json.load(json_file)
    print(f'加载反应编码字典完成，共有 {len(dict_rxn2id)} 个反应。')  # 打印加载的数据
    
print('Loading validation datasets feather path ...')


std_blast, metrics_blast, ec_no_rxn_blat =evTools.get_eval_results(baselineName= 'blast', dict_rxn2id=dict_rxn2id, method_type='ec')
evTools.display_html_results(metrics = metrics_blast, std_mean = std_blast, no_pred=ec_no_rxn_blat, eva_name ='Blast')