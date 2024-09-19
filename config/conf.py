'''
Author: Zhenkun Shi
Date: 2022-04-08 15:57:16
LastEditors: Zhenkun Shi
LastEditTime: 2023-10-06 17:11:07
FilePath: /preaction/config/conf.py
Description: 

Copyright (c) 2022 by tibd, All Rights Reserved. 
'''
import os

# CODE_MODE = 'PRODUCTION'
CODE_MODE = 'DEBUG'

#工程路径
DIR_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

#数据目录
DATA_ROOT =  DIR_PROJECT_ROOT + '/data/'

# 临时目录
TEMP_DIR = DIR_PROJECT_ROOT + '/temp/'
# 日志目录
LOG_DIR = DATA_ROOT + 'log/'


# 图片文件目录
PIC_DIR = f'production/files/pic' 

# 图片文件目录
RXN_JSON_DIR = f'production/files/json/' 

# 生产环境结果目录
DIR_PRODUCTION_RES = f'production/files/results/'


# 示例文件路径
SAMPLE_DIR = DATA_ROOT + 'samples/'

# WEB 用数据目录
WEB_DATA_DIR = DATA_ROOT +'web/'
UNIPROT_DATA_DIR = DATA_ROOT +'uniprot/'

SPLITER=';'

#Uniprot
URL_SPROT_SNAP201801 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2018_01/knowledgebase/uniprot_sprot-only2018_01.tar.gz'
URL_SPROT_SNAP202301 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2023_01/knowledgebase/uniprot_sprot-only2023_01.tar.gz'
URL_SPROT_SNAP202401 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2024_01/knowledgebase/uniprot_sprot-only2024_01.tar.gz'

FILE_SPROT_SNAP201801 = UNIPROT_DATA_DIR+'uniprot_sprot-only2018_01.tar.gz'
FILE_SPROT_SNAP202301 = UNIPROT_DATA_DIR+'uniprot_sprot-only2023_01.tar.gz'
FILE_SPROT_SNAP202401 = UNIPROT_DATA_DIR+'uniprot_sprot-only2024_01.tar.gz'


FILE_UNIPROT_PROTEIN_REACTION_RELATION= UNIPROT_DATA_DIR + 'sprot_rhea_relation.feather'


#RHEA
URL_RHEA_REACTION_SMILES = 'https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv'
URL_RHEA_REACTION_EC ='https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv'
FILE_RHEA_REACTION = UNIPROT_DATA_DIR +'rhea_reactions.feather'


# 特征目录
FEATURE_DIR = DATA_ROOT + 'common/featurebank/'

#ESM
DIR_FEATURE_ESM2_L33_650M = FEATURE_DIR + 'esm2l33650m/'
FILE_EMBD_PROTEIN_ESM2_L33_650M = FEATURE_DIR + 'esm2_l33_650m.feather'
FILE_COMPUTED_PROTIEN_FASTA_ESM2_L33_650M = DATA_ROOT + 'common/featurebank/esm2l33650m.fasta'

#T5
DIR_FEATURE_T5_SEQ = FEATURE_DIR + 't5/seq/'
DIR_FEATURE_T5_3DI = FEATURE_DIR + 't5/3di/'
FILE_EMBD_PROTEIN_T5_SEQ = FEATURE_DIR + 't5seq.feather'

#Unirep
DIR_FEATURE_UNIREP = FEATURE_DIR + 'unirep/seq/'
FILE_EMBD_PROTEIN_UNIREP = FEATURE_DIR + 'unirep.feather'


# 数据集
DIR_DATASET = f'{DATA_ROOT}datasets/task240524/'

FILE_DS_TRAIN = f'{DIR_DATASET}ds_train.feather'
FILE_DS_TEST = f'{DIR_DATASET}ds_test.feather'
FILE_DS_TRAIN_FASTA = f'{DIR_DATASET}ds_train.fasta'
FILE_DS_TEST_FASTA = f'{DIR_DATASET}ds_test.fasta'

FILE_DS_DICT_RXN2ID = f'{DIR_DATASET}dict_rxn2id.json'
FILE_DS_DICT_ID2RXN = f'{DIR_DATASET}dict_id2rxn.json'


FILE_DS_CASE_ECOLI = f'{DIR_DATASET}case/ds_case_ecoli.feather'


# 反应相关
FILE_DS_RHEA_REACTIONS = f'{DIR_DATASET}ds_rhea_reactions.feather'
FILE_DS_CHEBI_CPD = f'{DIR_DATASET}ds_chebi_cpd.feather'
FILE_DS_RHEA = f'{DIR_DATASET}/ds_rhea.feather'



# 字典目录
DIR_DICT = DATA_ROOT +'dict/'
DICT_LABEL_RHEA = DIR_DICT+'dict_label_rhea.npy'


# 原始数据
FILE_WEB_REACTIONS = f'{WEB_DATA_DIR}web_reactions.feather' #反应数据
FILE_WEB_REACTION_ENZYME_RELATION = f'{WEB_DATA_DIR}web_reaction_enzyme_relation.feather' #反应关系数据
FILE_WEB_PROTEIONS = f'{WEB_DATA_DIR}web_proteins.feather' #蛋白数据
FILE_WEB_EC = f'{WEB_DATA_DIR}web_ec.feather' #反应数据
FILE_SUP_SPROT = f'{UNIPROT_DATA_DIR}uniprot_sprot_info.feather' #uniprot 蛋白补充信息






# 训练好的模型

DIR_MODEL = DIR_PROJECT_ROOT + '/model/'


# 结果目录
RESULTS_DIR = DIR_PROJECT_ROOT + '/results/'

#结果文件
FILE_RESULTS_MSA = f'{RESULTS_DIR}baselines/exp_test_msa.tsv'
FILE_RESULTS_SIMI = f'{RESULTS_DIR}exp_test_simi.h5'

FILE_RESULTS_SIMI_ESM = f'{RESULTS_DIR}simi/exp_test_esm.h5'
FILE_RESULTS_SIMI_UNIREP = f'{RESULTS_DIR}simi/exp_test_unirep.h5'
FILE_RESULTS_SIMI_T5 = f'{RESULTS_DIR}simi/exp_test_t5.h5'

FILE_RESULTS_SIMI_ESM_REACTION = f'{RESULTS_DIR}simi/exp_test_esm_reaction.h5'
FILE_RESULTS_SIMI_UNIREP_REACTION = f'{RESULTS_DIR}simi/exp_test_unirep_reaction.h5'
FILE_RESULTS_SIMI_T5_REACTION = f'{RESULTS_DIR}simi/exp_test_t5_reaction.h5'


FILE_RESULTS_APLF = f'{RESULTS_DIR}exp_test_aplf.xlsx'
FILE_RESULTS_APLF_AP = f'{RESULTS_DIR}exp_test_aplf_ap.tsv'


#EC-based Results
FILE_RESULTS_DEEPEC = f'{RESULTS_DIR}baselines/exp_test_deepec.tsv'
FILE_RESULTS_CLEAN = f'{RESULTS_DIR}baselines/exp_test_clean.tsv'
FILE_RESULTS_ECRECER = f'{RESULTS_DIR}baselines/exp_test_ecrecer.tsv'
FILE_RESULTS_ECPRED = f'{RESULTS_DIR}baselines/exp_test_ecpred.tsv'
FILE_RESULTS_BLAST_EC = f'{RESULTS_DIR}baselines/exp_test_blast_ec.tsv'
FILE_RESULTS_CATFAM = f'{RESULTS_DIR}baselines/exp_test_catfam.tsv'
FILE_RESULTS_PRIAM = f'{RESULTS_DIR}baselines/exp_test_priam.tsv'




FILE_RESULTS_BLAST_DIRECT = f'{RESULTS_DIR}baselines/exp_test_blast_direct.tsv'
