import json

import pandas as pd

from rxnrecer.config import config as cfg
from rxnrecer.lib.llm import chat as llmchat
from rxnrecer.utils import bio_utils as butils
from rxnrecer.utils import file_utils as futils


def dataframe_apply(df, fn, axis=1):
    apply_fn = getattr(df, "parallel_apply", None)
    if apply_fn is None:
        apply_fn = df.apply
    return apply_fn(fn, axis=axis)


def prompt_selector(rxnrecer_s1, rxnrecer_s2):
    
    set_s1 = set(rxnrecer_s1.split(cfg.SPLITER))
    set_s2 = set(rxnrecer_s2.split(cfg.SPLITER))
    
    if set_s1 == set_s2:
        return 'prompt1'
    else:
        return 'prompt2'

def make_query_string(query_row, sys_prompt_dict):
    uniprot_id = query_row['input_id']
    seq = query_row['seq']
    system_prompt = sys_prompt_dict[query_row['prompt_sys']]['prompt_text']
    rxn_details = query_row['RXN_details']
    

    if rxn_details and len(rxn_details) > 0:

        if hasattr(rxn_details[0], 'reaction_id'):
            reaction_ids = [rxn.reaction_id for rxn in rxn_details if rxn is not None]
        else:

            reaction_ids = rxn_details
    else:
        reaction_ids = []
    
    query_info = {
        "protein information": {
            "uniprot id": uniprot_id,
            "protein amino acid sequence": seq
        },
        "reaction information": {f"predicted reaction {i+1}": rid for i, rid in enumerate(reaction_ids)}
    }
        
    res_dict = {
    "PROMPT_SYS": system_prompt,
    "PROMPT_USER": f"INPUT:\n{json.dumps(query_info, indent=2)}"
    }
    return res_dict


def make_query_batch(res_rxnrecer):
    sys_prompt_dict =futils.read_json_file(cfg.FILE_DICT_RXNRECERS3_PROMPT)
    rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
    res_rxnrecer = butils.get_rxn_details_batch(df_rxns=res_rxnrecer, rxn_bank=rxn_bank, rxn_id_column='RXNRECer-S2')
    res_rxnrecer['prompt_sys'] = dataframe_apply(res_rxnrecer, lambda x: prompt_selector(x['RXNRECer-S1'], x['RXNRECer-S2']), axis=1)
    res_rxnrecer['s3_query'] = dataframe_apply(res_rxnrecer, lambda x: make_query_string(x, sys_prompt_dict), axis=1)
    return res_rxnrecer

def batch_chat(res_rxnrecer, api_key=None, api_url=None, llm_model=None, debug=False):
    llm_model = llm_model or cfg.LLM_MODEL
    res = make_query_batch(res_rxnrecer)
    res['RXNRECER-S3'] = res.apply(lambda x: single_chat(x.s3_query, llm_model, api_key, api_url, debug)['results'], axis=1)
    res = res[['input_id', 'RXNRECer-S1', 'RXNRECer-S2', 'RXNRECer_with_prob', 'RXN_details', 'RXNRECER-S3']]
    return res

def single_chat(dict_query, llm_model=None, api_key=None, api_url=None, debug=False):

    chatcli = llmchat.Chat(name=llm_model, url=api_url, api_key=api_key)
    
    if debug:
        print(dict_query)

    try:
        response = chatcli.chat(
            message=dict_query['PROMPT_SYS'],
            system_prompt=dict_query['PROMPT_USER'],
            response_format={"type": "json_object"}
        )
        content_str = response.choices[0].message.content

        if debug:
            print("=== RAW RESPONSE ===\n", content_str)
        records = json.loads(content_str)
        return records

    except Exception as e:
        print(f"[ERROR] LLM query/parse failed: {e}")
        if debug and 'content_str' in locals():
            print("[RAW CONTENT]:", content_str)

        return {"results": [], "error": str(e)}


if __name__ == "__main__":
    
    model1 = 'openai/gpt-4.1'
    model2 = 'anthropic/claude-sonnet-4'    
    print('rxnrecer/lib/llm/qa.py')
    
