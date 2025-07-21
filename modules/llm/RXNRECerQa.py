import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../../')


import json
import hashlib
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from config import conf as cfg
from modules.llm import chat as llmchat
import modules.llm.prompt as prompt

APIKEY = os.environ.get("OPENROUTER_API_KEY")


def get_rxn_detail(rxn_id, rxn_bank):
    if not rxn_id or rxn_id == '-':
        return {
            'reaction id': '-',
            'reaction equation': '-'
        }

    match = rxn_bank[rxn_bank.reaction_id == rxn_id]
    if match.empty:
        return {}

    row = match.iloc[0]
    return {
        'reaction id': row.reaction_id,
        'reaction equation': row.equation,
        'reaction equation in ChEBI format': row.equation_chebi,
        'reaction equation in SMILES format': row.equation_smiles,
        'reaction associated Enzyme Commission Number': row.ec_number
    }


def make_user_query(uniprot_id, seq, rxn_ids, rxn_bank):
    rxns = rxn_ids.split(cfg.SPLITER) if rxn_ids and rxn_ids != '-' and cfg.SPLITER in rxn_ids else []
    reaction_info = {
        f"predicted reaction {i+1}": get_rxn_detail(rid.strip(), rxn_bank)
        for i, rid in enumerate(rxns)
    } if rxns else {}

    query_info = {
        "protein information": {
            "uniprot id": uniprot_id,
            "protein amino acid sequence": seq
        },
        "reaction information": reaction_info
    }

    if uniprot_id != '-' and uniprot_id.strip():
        system_prompt = prompt.p_sys_reranking_with_uniprotid
    elif reaction_info:
        system_prompt = prompt.p_sys_reranking_without_uniprotid
    else:
        system_prompt = prompt.p_sys_evidence_seeking

    return system_prompt, json.dumps(query_info, indent=2)


def session_chat(chat_bot, 
                 rxn_bank, 
                 seq, 
                 rxn_ids, 
                 model_name, 
                 uniprot_id='-',
                 debug=True, debug_level=1, cache_dir=None):
    cache_path = None
    
    model_name_short = model_name.split('/')[-1]

    if cache_dir:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = uniprot_id.strip() if uniprot_id and uniprot_id != '-' else f"SEQ_{hashlib.md5(seq.encode()).hexdigest()[:10]}"

        cache_path = cache_dir / f"{key}_{model_name_short}.json"

        if cache_path.exists():
            try:
                with open(cache_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"[WARN] Cache {cache_path} corrupted. Recomputing... {e}")
    

    system_prompt, input_str = make_user_query(uniprot_id, seq, rxn_ids, rxn_bank)

    if debug and debug_level >= 2:
        print("=== SYSTEM PROMPT ===\n", system_prompt)
        print("=== INPUT STRING ===\n", input_str)

    try:
        response = chat_bot.chat(
            message=input_str,
            system_prompt=system_prompt,
            response_format={"type": "json_object"}
        )
        content_str = response.choices[0].message.content

        if debug and debug_level >= 1:
            print("=== RAW RESPONSE ===\n", content_str)

        records = json.loads(content_str)

        if cache_path:
            with open(cache_path, 'w') as f:
                json.dump(records, f, indent=2)

        return records

    except Exception as e:
        print(f"[ERROR] LLM query/parse failed: {e}")
        if debug and 'content_str' in locals():
            print("[RAW CONTENT]:", content_str)
        return {"error": str(e)}
    
    
def single_chat(rxn_bank, seq, rxn_ids, model_name, uniprot_id='-', debug=True, debug_level=1):


    chatcli = llmchat.Chat(
        name=model_name,
        url='https://openrouter.ai/api/v1',
        api_key=APIKEY,
    )
    rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
        
    system_prompt, input_str = make_user_query(uniprot_id, seq, rxn_ids, rxn_bank)

    if debug and debug_level >= 2:
        print("=== SYSTEM PROMPT ===\n", system_prompt)
        print("=== INPUT STRING ===\n", input_str)

    try:
        response = chat_bot.chat(
            message=input_str,
            system_prompt=system_prompt,
            response_format={"type": "json_object"}
        )
        content_str = response.choices[0].message.content

        if debug and debug_level >= 1:
            print("=== RAW RESPONSE ===\n", content_str)

        records = json.loads(content_str)

        return records

    except Exception as e:
        print(f"[ERROR] LLM query/parse failed: {e}")
        if debug and 'content_str' in locals():
            print("[RAW CONTENT]:", content_str)
        return {"error": str(e)}


def chat_case(model_name):
    print(f'Use model {model_name} for case study')
    # data = pd.read_pickle('/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/llm/res/crop_newfinding_20250624.pkl')
    data = pd.read_pickle('/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/llm/res/crop_nonenzyme2enzyme_20250630.pkl')
    rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)


    cache_dir = Path(cfg.TEMP_DIR) / 'rxnrecer-s3'
    cache_dir.mkdir(parents=True, exist_ok=True)

    chatcli = llmchat.Chat(
        name=model_name,
        url='https://openrouter.ai/api/v1',
        api_key=APIKEY,
    )

    results = []
    for idx, row in tqdm(data.iterrows(), total=len(data), desc="LLM RXNRECer-S3 Inference"):
        try:
            result = session_chat(
                chat_bot=chatcli,
                rxn_bank=rxn_bank,
                seq=row['seq'],
                rxn_ids=row['RXNRECer'],
                uniprot_id=row.get('uniprot_id', '-'),
                debug=False,
                model_name=model_name,
                cache_dir=cache_dir
            )
        except Exception as e:
            print(f"[ERROR] Index {idx} ({row.get('uniprot_id', '-')}) failed: {e}")
            result = {"error": str(e)}
        results.append(result)

    data['RXNRECer-S3'] = results
    return data



def format_s3(s3_str):
    rxn_s3 = '-'
    json_s3 ={}
    if 'results' in s3_str:
        s3_str = s3_str['results']
        df = pd.DataFrame(s3_str)
    else:
        if isinstance(s3_str, list):
            df = pd.DataFrame(s3_str)
        elif isinstance(s3_str, dict):
            df = pd.DataFrame.from_dict(s3_str,orient='index').T
        else:
            df = pd.DataFrame.from_dict([s3_str])
    # df = pd.DataFrame(s3_str)
    if len(df)>0:
        try:
            df =df[df.selected=='yes']
            rxn_s3 = cfg.SPLITER.join(df['reaction_id'].values)
            json_s3 = df.to_json(orient='records')
        except Exception as e:
            rxn_s3 = '-'
            json_s3 ={}
            print(f"[ERROR] Format S3 failed: {e}")
    return rxn_s3, json_s3

if __name__ == "__main__":
    
    model1 = 'openai/gpt-4.1'
    model2 = 'anthropic/claude-sonnet-4'
    
    model = model2
    
    # df = chat_case(model_name = model1)
    df = chat_case(model_name = model)
    
    df.to_pickle(f"{cfg.TEMP_DIR}/crop_gpt_rxnrecers3-{model.split('/')[-1]}.pkl")
    print(df)
