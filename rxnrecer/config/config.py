'''
Author: Zhenkun Shi
Date: 2022-04-08 15:57:16
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-01-28
FilePath: /preaction/config/conf.py
Description: RXNRECer configuration file - cleaned version

Copyright (c) 2022 by tibd, All Rights Reserved. 
'''
import os
from rxnrecer.utils import file_utils

CODE_MODE = 'PRODUCTION'
# CODE_MODE = 'DEBUG'
# Cache configuration
CACHE_ENABLED = True

# Project path
DIR_PROJECT_ROOT = str(file_utils.get_project_root())

# Data directories
DATA_ROOT = DIR_PROJECT_ROOT + '/data/'

# Cache directory
CACHE_DIR = DIR_PROJECT_ROOT + '/results/cache/'

# 临时目录
TEMP_DIR = DIR_PROJECT_ROOT + '/temp/'

# 示例文件路径
SAMPLE_DIR = DATA_ROOT + 'sample/'

# Model directory
CKPT_PROSTT5 = DIR_PROJECT_ROOT + '/ckpt/prostt5'

# 生产环境Model files
FILE_MOLEL_PRODUCTION_BEST_MODEL = DIR_PROJECT_ROOT + '/ckpt/rxnrecer/production_185846best.pth'

# Separator
SPLITER = ';'

# 核心数据文件
FILE_RHEA_REACTION = DATA_ROOT + 'rhea/rhea_reactions.feather'
FILE_DS_DICT_ID2RXN = DATA_ROOT + 'dict/dict_id2rxn.json'
FILE_DS_DICT_RXN2ID = DATA_ROOT + 'dict/dict_rxn2id.json'
FILE_DS_RHEA_REACTIONS = DATA_ROOT + 'rhea/ds_rhea_reactions.feather'
FILE_DS_TRAIN = DATA_ROOT + 'datasets/task240524/ds_train.feather'
FILE_DS_TEST = DATA_ROOT + 'datasets/task240524/ds_test.feather'
FILE_DS_TRAIN_FASTA = DATA_ROOT + 'datasets/task240524/ds_train.fasta'
FILE_DS_TEST_FASTA = DATA_ROOT + 'datasets/task240524/ds_test.fasta'
FILE_DS_DMND = DATA_ROOT + 'datasets/task240524/ds_train.dmnd'
FILE_DS_CHEBI_CPD = DATA_ROOT + 'chebi/ds_chebi_cpd.feather'

# 特征文件
FILE_PRODUCTION_FEATURES = DATA_ROOT + 'feature_bank/featureBank.feather'
FILE_PRODUCTION_FEATURES_T5 = DIR_PROJECT_ROOT + '/files/features/featureBank_t5.feather'

# 字典文件
DICT_UNIPROT_RHEA = DATA_ROOT + 'dict/dict_uniprot_rhea.json'
DICT_EC_RHEA = DATA_ROOT + 'dict/dict_ec_rhea.json'
FILE_DICT_RXNRECERS3_PROMPT = DATA_ROOT + 'dict/dict_rxnrecers3_prompt.json'

# 反应和化合物相关目录
DIR_RXN_JSON = DATA_ROOT + 'rxn_json/'
DIR_CPD_SVG = DATA_ROOT + 'cpd_svg/'

# 数据集目录
DIR_DATASET = DATA_ROOT + 'datasets/task240524/'

# Results directory
RESULTS_DIR = DIR_PROJECT_ROOT + '/results/'

# LLM API配置
LLM_API_URL = os.environ.get('LLM_API_URL', '')
LLM_API_KEY = os.environ.get('LLM_API_KEY', '')



# Uniprot相关
UNIPROT_DATA_DIR = DATA_ROOT + 'uniprot/'
FILE_UNIPROT_PROTEIN_REACTION_RELATION = UNIPROT_DATA_DIR + 'sprot_rhea_relation.feather'
URL_SPROT_SNAP201801 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2018_01/knowledgebase/uniprot_sprot-only2018_01.tar.gz'
URL_SPROT_SNAP202401 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2024_01/knowledgebase/uniprot_sprot-only2024_01.tar.gz'
FILE_SPROT_SNAP201801 = UNIPROT_DATA_DIR + 'uniprot_sprot-only2018_01.tar.gz'
FILE_SPROT_SNAP202401 = UNIPROT_DATA_DIR + 'uniprot_sprot-only2024_01.tar.gz'

# RHEA相关
RHEA_DATA_DIR = DATA_ROOT + 'rhea/'
URL_RHEA_REACTION_SMILES = 'https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv'
URL_RHEA_REACTION_EC = 'https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv'

# CHEBI相关
CHEBI_DATA_DIR = DATA_ROOT + 'chebi/'

# Web数据
WEB_DATA_DIR = DATA_ROOT + 'web/'
FILE_WEB_REACTIONS = WEB_DATA_DIR + 'web_reactions.feather'
FILE_WEB_REACTION_ENZYME_RELATION = WEB_DATA_DIR + 'web_reaction_enzyme_relation.feather'
FILE_WEB_PROTEIONS = WEB_DATA_DIR + 'web_proteins.feather'
FILE_WEB_EC = WEB_DATA_DIR + 'web_ec.feather'
FILE_SUP_SPROT = UNIPROT_DATA_DIR + 'uniprot_sprot_info.feather'




