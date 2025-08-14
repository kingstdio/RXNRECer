'''
Author: Zhenkun Shi
Date: 2022-04-08 15:57:16
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-05-28 09:55:41
FilePath: /preaction/config/conf.py
Description: 

Copyright (c) 2022 by tibd, All Rights Reserved. 
'''
import os

# CODE_MODE = 'PRODUCTION'
CODE_MODE = 'DEBUG'

#工程路径
DIR_PROJECT_ROOT = '/hpcfs/fhome/shizhenkun/codebase/rxnrecer_production'

#数据目录
DATA_ROOT =  DIR_PROJECT_ROOT + '/data/'

# 临时目录
TEMP_DIR = DIR_PROJECT_ROOT + '/temp/'
# 日志目录
LOG_DIR = DATA_ROOT + 'log/'


# 图片文件目录
PIC_DIR = f'production/files/pic' 

# 图片文件目录
DIR_RXN_JSON = f'files/rxn_json/' 
DIR_CPD_SVG = f'files/cpd_svg/' 
# 生产环境结果目录
DIR_PRODUCTION_RES = f'files/results/'
FILE_PRODUCTION_FEATURES = f'{DIR_PROJECT_ROOT}/files/features/featureBank.feather'
FILE_PRODUCTION_FEATURES_T5 = f'{DIR_PROJECT_ROOT}/files/features/featureBank_t5.feather'


FILE_MOLEL_PRODUCTION_BEST_MODEL = f'{DATA_ROOT}model/production_185846best.pth'


# 示例文件路径
SAMPLE_DIR = DATA_ROOT + 'samples/'


#PDB 结构文件目录
DIR_PDB_BEST = DATA_ROOT + 'structure/pdb/rcsb2uniprot/'

# WEB 用数据目录
WEB_DATA_DIR = DATA_ROOT +'web/'
UNIPROT_DATA_DIR = DATA_ROOT +'uniprot/'
RHEA_DATA_DIR = DATA_ROOT +'rhea/'

SPLITER=';'

#Uniprot
URL_SPROT_SNAP201801 = 'https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_major_releases/release-2018_01/knowledgebase/uniprot_sprot-only2018_01.tar.gz'
URL_SPROT_SNAP202301 = 'https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_major_releases/release-2023_01/knowledgebase/uniprot_sprot-only2023_01.tar.gz'
URL_SPROT_SNAP202401 = 'https://ftp.ebi.ac.uk/pub/databases/uniprot/previous_major_releases/release-2024_01/knowledgebase/uniprot_sprot-only2024_01.tar.gz'

FILE_SPROT_SNAP201801 = UNIPROT_DATA_DIR+'uniprot_sprot-only2018_01.tar.gz'
FILE_SPROT_SNAP202301 = UNIPROT_DATA_DIR+'uniprot_sprot-only2023_01.tar.gz'
FILE_SPROT_SNAP202401 = UNIPROT_DATA_DIR+'uniprot_sprot-only2024_01.tar.gz'


FILE_UNIPROT_PROTEIN_REACTION_RELATION= UNIPROT_DATA_DIR + 'sprot_rhea_relation.feather'


#RHEA
URL_RHEA_REACTION_SMILES = 'https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv'
URL_RHEA_REACTION_EC ='https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv'
FILE_RHEA_REACTION = RHEA_DATA_DIR +'rhea_reactions.feather'


# 特征目录
FEATURE_DIR = DATA_ROOT + 'featurebank/'

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
DIR_DATASET_3DI = f'{DATA_ROOT}3di/'

FILE_DS_TRAIN = f'{DIR_DATASET}ds_train.feather'
FILE_DS_TEST = f'{DIR_DATASET}ds_test.feather'
FILE_DS_TRAIN_FASTA = f'{DIR_DATASET}ds_train.fasta'
FILE_DS_TEST_FASTA = f'{DIR_DATASET}ds_test.fasta'

FILE_DS_DICT_RXN2ID = f'{DATA_ROOT}dict/dict_rxn2id.json'
FILE_DS_DICT_ID2RXN = f'{DATA_ROOT}dict/dict_id2rxn.json'

FILE_DS_DMND = f'{DATA_ROOT}datasets/task240524/ds_train.dmnd'

FILE_DS_PDB_LIST = f'{DIR_DATASET}/ds_all_pdb_map.feather' # 所有pdb 结构列表
FILE_DS_3DI_LIST = f'{DIR_DATASET}/ds_all_3di.feather' # 所有pdb 结构列表

# 3DI
FILE_EMBD_PROTEIN_TDIT5 = f'{DIR_DATASET}ds_all_3di_embedding.feather' # 所有pdb 结构列表

FILE_DS_CASE_ECOLI = f'{DATA_ROOT}datasets/case/ds_case_ecoli.feather'


# 反应相关
FILE_DS_RHEA_REACTIONS = f'{DIR_DATASET}ds_rhea_reactions.feather'
FILE_DS_CHEBI_CPD = f'{DIR_DATASET}ds_chebi_cpd.feather'
FILE_DS_RHEA = f'{DIR_DATASET}ds_rhea.feather'



# 字典目录
DIR_DICT = DATA_ROOT +'dict/'
DICT_LABEL_RHEA = DIR_DICT+'dict_label_rhea.npy'
DICT_UNIPROT_RHEA = DIR_DICT +'dict_uniprot_rhea.json'
DICT_EC_RHEA = DIR_DICT +'dict_ec_rhea.json'
DICT_RHEA_EC = DIR_DICT +'dict_rhea_ec.json'


# 原始数据
FILE_WEB_REACTIONS = f'{WEB_DATA_DIR}web_reactions.feather' #反应数据
FILE_WEB_REACTION_ENZYME_RELATION = f'{WEB_DATA_DIR}web_reaction_enzyme_relation.feather' #反应关系数据
FILE_WEB_PROTEIONS = f'{WEB_DATA_DIR}web_proteins.feather' #蛋白数据
FILE_WEB_EC = f'{WEB_DATA_DIR}web_ec.feather' #反应数据
FILE_SUP_SPROT = f'{UNIPROT_DATA_DIR}uniprot_sprot_info.feather' #uniprot 蛋白补充信息


#FOLD SEEK 数据
DIR_FOLDSEEK_PDB = f'{DIR_PROJECT_ROOT}/results/intermediate/foldseek/'


# 训练好的模型
DIR_MODEL = DIR_PROJECT_ROOT + '/data/model/'
FILE_WEIGHT_PRODUCTION_BEST_MODEL = f'{DIR_MODEL}production_185846best.pth'


# 结果目录
RESULTS_DIR = DIR_PROJECT_ROOT + '/results/'

#结果文件
FILE_RESULTS_MSA = f'{RESULTS_DIR}baselines/exp_test_msa.tsv'
FILE_RESULTS_SIMI = f'{RESULTS_DIR}exp_test_simi.h5'
FILE_RESULTS_SIMI_ESM = f'{RESULTS_DIR}simi/exp_test_esm.h5'
FILE_RESULTS_SIMI_UNIREP = f'{RESULTS_DIR}simi/exp_test_unirep.h5'
FILE_RESULTS_SIMI_T5 = f'{RESULTS_DIR}simi/exp_test_t5.h5'
FILE_RESULTS_APLF = f'{RESULTS_DIR}exp_test_aplf.xlsx'
FILE_RESULTS_APLF_AP = f'{RESULTS_DIR}exp_test_aplf_ap.tsv'


#EC-based Results
DIR_RES_BASELINE = f'{DIR_PROJECT_ROOT}/baselines/'

# Case RESULTS 2018later 
CASE_DIR = f'{DIR_PROJECT_ROOT}/case/'
CASE_2018LATER = f'{CASE_DIR}2018later/'
# EC-based methods
FILE_CASE_RESULTS_BLAST_EC = f'{CASE_2018LATER}res/exp_test_blast_ec.tsv'
FILE_CASE_RESULTS_DEEPEC = f'{CASE_2018LATER}res/exp_test_deepec.tsv'
FILE_CASE_RESULTS_CLEAN = f'{CASE_2018LATER}res/exp_test_clean.tsv'
FILE_CASE_RESULTS_ECRECER = f'{CASE_2018LATER}res/exp_test_ecrecer.tsv'
FILE_CASE_RESULTS_ECPRED = f'{CASE_2018LATER}res/exp_test_ecpred.tsv'
FILE_CASE_RESULTS_CATFAM = f'{CASE_2018LATER}res/exp_test_catfam.tsv'
FILE_CASE_RESULTS_PRIAM = f'{CASE_2018LATER}res/exp_test_priam.tsv'

#Direct methods
FILE_CASE_RESULTS_BLAST_DIRECT = f'{CASE_2018LATER}res/exp_test_blast_direct.tsv'
FILE_CASE_RESULTS_SIMI_ESM_REACTION = f'{CASE_2018LATER}res/exp_test_esm_reaction.h5'
FILE_CASE_RESULTS_SIMI_UNIREP_REACTION = f'{CASE_2018LATER}res/exp_test_unirep_reaction.h5'
FILE_CASE_RESULTS_SIMI_T5_REACTION = f'{CASE_2018LATER}res/exp_test_t5_reaction.h5'
FILE_CASE_RESULTS_RXNRECER_REACTION = f'{CASE_2018LATER}res/exp_test_rxnrecer_reaction.pkl'
