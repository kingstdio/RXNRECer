import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/../')
from rxnrecer.config import config as cfg
import pandas as pd
import numpy as np
import json
from rxnrecer.utils import file_utils as ftool
import tempfile
import subprocess
from rxnrecer.lib.embedding import esm_embedding as ebdseq
from rxnrecer.lib.embedding import t5_embedding as ebdt5
from rxnrecer.lib.smi import  simi_caculator as simitool


def getmsa(df_test, k=1):
    
    print('Running MSA ...')
    with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.tsv') as res_blast:


        ftool.table2fasta(df_test, fasta_test.name)

        cmd = [f"{cfg.DIR_PROJECT_ROOT}/extools/msa/diamond", "blastp",  "-d", cfg.FILE_DS_DMND, "-q", fasta_test.name, "-o", res_blast.name, "-b5", "-c1", "-k", str(k), "--quiet"]


        subprocess.run(cmd, check=True)


        blast_res = pd.read_csv(res_blast.name, sep='\t', names=['id', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        

        with open(cfg.DICT_UNIPROT_RHEA, "r") as json_file:
            dict_uniprot2rhea = json.load(json_file)

        blast_res['rxn_msa'] = blast_res.sseqid.apply(lambda x: dict_uniprot2rhea.get(x, 'None'))
        blast_res = blast_res[['id', 'rxn_msa']].rename(columns={'id':'uniprot_id'})
        

        blast_res = df_test[['uniprot_id']].merge(blast_res, on='uniprot_id', how='left').fillna('NO-PREDICTION')
        
        
    return blast_res  

def getcatfam(df_test):
    
    print('Running CatFam ...')
    with tempfile.NamedTemporaryFile(delete=False, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=False, suffix='.tsv') as res_catfam:
        

        ftool.table2fasta(df_test, fasta_test.name)
        

        cmd = [
            "singularity", "exec", f"{cfg.DIR_PROJECT_ROOT}/extools/ec/catfam.sif",
            "/catfam/source/catsearch.pl",
            "-d", "/catfam/CatFamDB/CatFam_v2.0/CatFam4D99R",
            "-i", fasta_test.name,
            "-o", res_catfam.name
        ]


        subprocess.run(cmd, check=True)


        res_catfam = pd.read_csv(res_catfam.name, sep='\t', names=['uniprot_id', 'ec_catfam']).fillna('-')
        

        with open(cfg.DICT_EC_RHEA, "r") as json_file:
            dict_ec2rhea = json.load(json_file)

        res_catfam['rxn_catfam'] = res_catfam.ec_catfam.apply(lambda x: dict_ec2rhea.get(x, 'None'))
        

        res_catfam = df_test[['uniprot_id']].merge(res_catfam, on='uniprot_id', how='left').fillna('NO-PREDICTION')

    return res_catfam



def getecrecer(df_test):
    
    print('Running ECRECer ...')
    with tempfile.NamedTemporaryFile(delete=False, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=False, suffix='.tsv') as res_ecrecer:
        

        ftool.table2fasta(df_test, fasta_test.name)
        

        cmd = [
            "singularity", "exec", "--nv", f"{cfg.DIR_PROJECT_ROOT}/extools/ec/ecrecer.sif",
            "python", "/ecrecer/production.py",
            "-i", fasta_test.name,
            "-o", res_ecrecer.name,
            "-mode", "p",
            "-topk", "20"
        ]
        

        subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


        res_df = pd.read_csv(res_ecrecer.name, sep='\t', names=['uniprot_id', 'ec_ecrecer']).fillna('-')


        with open(cfg.DICT_EC_RHEA, "r") as json_file:
            dict_ec2rhea = json.load(json_file)


        res_df['rxn_ecrecer'] = res_df['ec_ecrecer'].apply(lambda x: dict_ec2rhea.get(x, 'None'))
        

        res_df = df_test[['uniprot_id']].merge(res_df, on='uniprot_id', how='left').fillna('NO-PREDICTION')

    return res_df



# Function to calculate similarity between protein features
def get_top_protein_simi(x_feature, y_feature, y_uniprot_id, dict_featureBank,dict_uniprot2rhea, topk):
    future_cosine = simitool.get_cosine_similarity(x_feature, y_feature)
    future_cosine = future_cosine.transpose()   
    topk_indices = np.argsort(future_cosine, axis=1)[:, -topk:][:, ::-1]
    topk_values = np.round(np.take_along_axis(future_cosine, topk_indices, axis=1), 6)
    

    simi_uniprot_id = list(map(dict_featureBank.get, topk_indices.flatten()))
    simi_rhea = list(map(dict_uniprot2rhea.get, simi_uniprot_id))
    result_matrix = np.empty(len(simi_uniprot_id), dtype=object)

    result_matrix[:] = list(zip(simi_rhea, topk_values.flatten()))
    # Reshape to the original shape
    result_matrix = result_matrix.reshape(topk_indices.shape)
    # Create the DataFrame
    res = pd.DataFrame({
        'uniprot_id': y_uniprot_id,
        'simi': list(result_matrix)
    })

    return res


def getT5(df_test, topk=3):
    
    print('Running T5 ...')
    featureBank = pd.read_feather(cfg.FILE_PRODUCTION_FEATURES)
    dict_featureBank = pd.Series( featureBank['uniprot_id'], featureBank.index.values).to_dict()
    

    with open(cfg.DICT_UNIPROT_RHEA, "r") as json_file:
        dict_uniprot2rhea = json.load(json_file)

    # ESM embedding   
    embd_esm = ebdseq.getEsm(df_test.rename(columns={'uniprot_id':'id'}))
    # ESM similarity
    esm_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.esm), 
                                    y_feature=np.vstack(embd_esm.esm), 
                                    y_uniprot_id=embd_esm.id, 
                                    dict_featureBank=dict_featureBank, 
                                    dict_uniprot2rhea = dict_uniprot2rhea,
                                    topk=topk).rename(columns={'simi':'esm'})
        
    # T5 Embedding    
    embd_t5 = ebdt5.get_embd_seq(seqdfwithid=df_test.rename(columns={'uniprot_id':'id'}), batch_size=20)
        
    # T5 similarity
    t5_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.t5), 
                                y_feature=np.vstack(embd_t5.t5), 
                                y_uniprot_id=embd_t5.id, 
                                dict_featureBank=dict_featureBank, 
                                dict_uniprot2rhea = dict_uniprot2rhea,
                                topk=topk).rename(columns={'simi':'t5'})
    
    res = esm_cos.merge(t5_cos, on='uniprot_id', how='left')
    return res



if __name__ =='__main__':
    pass