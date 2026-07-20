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
from pathlib import Path
from rxnrecer.utils import file_utils

CODE_MODE = 'PRODUCTION'
# Cache configuration
CACHE_ENABLED = True

# Project path
DIR_PROJECT_ROOT = str(file_utils.get_project_root())
_PROJECT_ROOT_PATH = Path(DIR_PROJECT_ROOT)


def _first_existing_path(*relative_paths: str) -> str:
    """Return the first existing path, otherwise the first candidate."""
    candidates = [_PROJECT_ROOT_PATH / rel for rel in relative_paths]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return str(candidates[0])

# Data directories
DATA_ROOT = DIR_PROJECT_ROOT + '/data/'

# Cache directory
CACHE_DIR = DIR_PROJECT_ROOT + '/results/cache/'

# Temporary directory
TEMP_DIR = DIR_PROJECT_ROOT + '/temp/'

# Sample file path
SAMPLE_DIR = DATA_ROOT + 'sample/'

# Model directory
CKPT_PROSTT5 = DIR_PROJECT_ROOT + '/ckpt/prostt5'

# Production model files
FILE_MOLEL_PRODUCTION_BEST_MODEL = DIR_PROJECT_ROOT + '/ckpt/rxnrecer/production_185849best.pth'

# Separator
SPLITER = ';'

# Core data files
FILE_RHEA_REACTION = DATA_ROOT + 'rhea/rhea_reactions.feather'
FILE_DS_DICT_ID2RXN = DATA_ROOT + 'dict/dict_id2rxn.json'
FILE_DS_DICT_RXN2ID = DATA_ROOT + 'dict/dict_rxn2id.json'
FILE_DS_RHEA_REACTIONS = DATA_ROOT + 'rhea/ds_rhea_reactions.feather'
FILE_DS_TRAIN = DATA_ROOT + 'datasets/task240524/ds_train.feather'
FILE_DS_RCP = DATA_ROOT + 'datasets/task240524/ds_test.feather'
FILE_DS_TRAIN_FASTA = DATA_ROOT + 'datasets/task240524/ds_train.fasta'
FILE_DS_RCP_FASTA = DATA_ROOT + 'datasets/task240524/ds_test.fasta'
FILE_DS_DMND = DATA_ROOT + 'datasets/task240524/ds_train.dmnd'
FILE_DS_CHEBI_CPD = DATA_ROOT + 'chebi/ds_chebi_cpd.feather'

# Feature files
FILE_PRODUCTION_FEATURES = DATA_ROOT + 'feature_bank/featureBank.feather'
FILE_PRODUCTION_FEATURES_T5 = DATA_ROOT + 'legacy/features/featureBank_t5.feather'
FEATURE_DIR = _first_existing_path('data/featurebank', 'data/feature_bank') + '/'
DIR_FEATURE_ESM2_L33_650M = FEATURE_DIR + 'esm2l33650m/'
DIR_FEATURE_T5_SEQ = FEATURE_DIR + 't5/seq/'
DIR_FEATURE_T5_3DI = FEATURE_DIR + 't5/3di/'
DIR_FEATURE_UNIREP = FEATURE_DIR + 'unirep/seq/'
FILE_EMBD_PROTEIN_ESM2_L33_650M = _first_existing_path(
    'data/featurebank/esm2_l33_650m.feather',
    'data/feature_bank/esm2_l33_650m.feather',
)
FILE_EMBD_PROTEIN_T5_SEQ = _first_existing_path(
    'data/featurebank/t5seq.feather',
    'data/feature_bank/t5seq.feather',
)
FILE_EMBD_PROTEIN_UNIREP = _first_existing_path(
    'data/featurebank/unirep.feather',
    'data/feature_bank/unirep.feather',
)
FILE_COMPUTED_PROTIEN_FASTA_ESM2_L33_650M = DATA_ROOT + 'common/featurebank/esm2l33650m.fasta'

# Dictionary files
DICT_UNIPROT_RHEA = DATA_ROOT + 'dict/dict_uniprot_rhea.json'
DICT_EC_RHEA = DATA_ROOT + 'dict/dict_ec_rhea.json'
FILE_DICT_RXNRECERS3_PROMPT = DATA_ROOT + 'dict/dict_rxnrecers3_prompt.json'

# Reaction and compound related directories
DIR_RXN_JSON = DATA_ROOT + 'rxn_json/'
DIR_CPD_SVG = DATA_ROOT + 'cpd_svg/'

# Dataset directory
DIR_DATASET = DATA_ROOT + 'datasets/task240524/'
DIR_DATASET_3DI = DATA_ROOT + '3di/'

# Results directory
RESULTS_DIR = DIR_PROJECT_ROOT + '/results/'
DIR_RES_BASELINE = DIR_PROJECT_ROOT + '/baselines/'
CASE_DIR = DIR_PROJECT_ROOT + '/case/'
DIR_FOLDSEEK_PDB = DIR_PROJECT_ROOT + '/results/intermediate/foldseek/'
DIR_MODEL = DIR_PROJECT_ROOT + '/ckpt/rxnrecer/'

# Legacy compatibility aliases
FILE_DS_PDB_LIST = DIR_DATASET + 'ds_all_pdb_map.feather'
FILE_DS_3DI_LIST = DIR_DATASET + 'ds_all_3di.feather'
FILE_EMBD_PROTEIN_TDIT5 = DIR_DATASET + 'ds_all_3di_embedding.feather'
FILE_DS_RHEA = DIR_DATASET + 'ds_rhea.feather'
FILE_DS_RHEA_REACTIONS = DIR_DATASET + 'ds_rhea_reactions.feather'
FILE_DS_CASE_ECOLI = DATA_ROOT + 'datasets/case/ds_case_ecoli.feather'
DIR_DICT = DATA_ROOT + 'dict/'
DICT_LABEL_RHEA = DIR_DICT + 'dict_label_rhea.npy'
DICT_RHEA_EC = DIR_DICT + 'dict_rhea_ec.json'
FILE_RESULTS_SIMI = RESULTS_DIR + 'simi/exp_test_simi.pkl'
FILE_RESULTS_SIMI_ESM = RESULTS_DIR + 'simi/exp_test_esm.pkl'
FILE_RESULTS_SIMI_UNIREP = RESULTS_DIR + 'simi/exp_test_unirep.pkl'
FILE_RESULTS_SIMI_T5 = RESULTS_DIR + 'simi/exp_test_t5.pkl'

# LLM API configuration
LLM_API_URL = os.environ.get('LLM_API_URL', '')
LLM_API_KEY = os.environ.get('LLM_API_KEY', '')
LLM_MODEL = os.environ.get('LLM_MODEL', 'openai/gpt-4.1')


# UniProt related
UNIPROT_DATA_DIR = DATA_ROOT + 'uniprot/'
FILE_UNIPROT_PROTEIN_REACTION_RELATION = UNIPROT_DATA_DIR + 'sprot_rhea_relation.feather'
URL_SPROT_SNAP201801 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2018_01/knowledgebase/uniprot_sprot-only2018_01.tar.gz'
URL_SPROT_SNAP202401 = 'https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2024_01/knowledgebase/uniprot_sprot-only2024_01.tar.gz'
FILE_SPROT_SNAP201801 = UNIPROT_DATA_DIR + 'uniprot_sprot-only2018_01.tar.gz'
FILE_SPROT_SNAP202401 = UNIPROT_DATA_DIR + 'uniprot_sprot-only2024_01.tar.gz'

# RHEA related
RHEA_DATA_DIR = DATA_ROOT + 'rhea/'
URL_RHEA_REACTION_SMILES = 'https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv'
URL_RHEA_REACTION_EC = 'https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv'

# CHEBI related
CHEBI_DATA_DIR = DATA_ROOT + 'chebi/'

# Web data
WEB_DATA_DIR = DATA_ROOT + 'web/'
FILE_WEB_REACTIONS = WEB_DATA_DIR + 'web_reactions.feather'
FILE_WEB_REACTION_ENZYME_RELATION = WEB_DATA_DIR + 'web_reaction_enzyme_relation.feather'
FILE_WEB_PROTEIONS = WEB_DATA_DIR + 'web_proteins.feather'
FILE_WEB_EC = WEB_DATA_DIR + 'web_ec.feather'
FILE_SUP_SPROT = UNIPROT_DATA_DIR + 'uniprot_sprot_info.feather'
