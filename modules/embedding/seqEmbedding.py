import sys, os
import pandas as pd
import numpy as np
import re
from transformers import T5Tokenizer, T5EncoderModel
from tqdm import tqdm
import torch
import esm
from jax_unirep import get_reps

# Adding the project directory to system path for imports
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1, '../')



def getUnirep(sequences, batch_size=10):
    """
    This function computes the unirep (h_avg) for each sequence in the input dataframe.
    The sequences are processed in batches to improve efficiency.
    
    Args:
        sequences (pd.DataFrame): DataFrame with columns 'id' and 'seq'.
        batch_size (int): The batch size for processing sequences.

    Returns:
        pd.DataFrame: DataFrame with 'id' and corresponding 'unirep' (h_avg).
    """
    # Initialize an empty list to store the batch results
    result_batches = []

    # Process sequences in batches
    for start in tqdm(range(0, len(sequences), batch_size)):
        end = min(start + batch_size, len(sequences))
        batch = sequences.iloc[start:end]
        
        # Zip id and seq together for the current batch
        batch_ids = batch['id'].tolist()
        batch_seqs = batch['seq'].tolist()

        # Call get_reps for the current batch
        h_avg, _, _ = get_reps(batch_seqs)
        
        # Combine the results (id, h_avg) for the current batch
        batch_result = pd.DataFrame({
            'id': batch_ids,
            'unirep': [h_avg[i] for i in range(len(h_avg))]
        })

        # Append the current batch result to the list
        result_batches.append(batch_result)

    # Concatenate all batch results into a single DataFrame
    final_result = pd.concat(result_batches, ignore_index=True)
    
    return final_result


# region 将字符串拆分成固定长度
def cut_text(text,lenth):
    """[将字符串拆分成固定长度]

    Args:
        text ([string]): [input string]
        lenth ([int]): [sub_sequence length]

    Returns:
        [string list]: [string results list]
    """
    textArr = re.findall('.{'+str(lenth)+'}', text) 
    textArr.append(text[(len(textArr)*lenth):]) 
    return textArr 
#endregion

#region 对单个序列进行embedding
def get_rep_single_seq_inter(seqid, sequence, model, batch_converter,rep_layers=[33], seqthres=1022, device=0):
    """[对单个序列进行embedding]

    Args:
        seqid ([string]): [sequence name]]
        sequence ([sting]): [sequence]
        model ([model]): [ embedding model]]
        batch_converter ([object]): [description]
        seqthres (int, optional): [max sequence length]. Defaults to 1022.

    Returns:
        [type]: [description]
    """
    
    if len(sequence) < seqthres:
        data =[(seqid, sequence)]
    else:
        seqArray = cut_text(sequence, seqthres)
        data=[]
        for item in seqArray:
            data.append((seqid, item))
            
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.to(device=f'cuda:{device}', non_blocking=True)
        
    MINI_SIZE = len(batch_labels)
    
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=rep_layers, return_contacts=False)    
    

    representations = {layer: t.to(device="cpu") for layer, t in results["representations"].items()}
    result ={}
    result["label"] = batch_labels[0]

    for i in range(MINI_SIZE):
        if i ==0:
            result["mean_representations"] = {layer: t[i, 1 : len(batch_strs[0]) + 1].mean(0).clone() for layer, t in representations.items()}
        else:
            for index, layer in enumerate(rep_layers):
                result["mean_representations"][layer] += {layer: t[i, 1 : len(batch_strs[0]) + 1].mean(0).clone() for layer, t in representations.items()}[layer]

    for index, layer in enumerate(rep_layers):
        result["mean_representations"][layer] = result["mean_representations"][layer] /MINI_SIZE
    
    return result['mean_representations'][33]
#endregion


#region 对多个序列进行embedding
def getEsm(sequences, seqthres=1022, device=0):
    """[对多个序列进行embedding]
    Args:
        sequences ([DataFrame]): [ sequence info]]
        seqthres (int, optional): [description]. Defaults to 1022.

    Returns:
        [DataFrame]: [final_rep0, final_rep32, final_rep33]
    """

    final_rep33 =[]
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    if torch.cuda.is_available():
            model = model.cuda(device)
            print("Transferred model to GPU")

    for i in tqdm(range(len(sequences))):
        rep33 = get_rep_single_seq_inter(seqid = sequences.iloc[i].id, 
                                         sequence=sequences.iloc[i].seq, 
                                         model=model, 
                                         batch_converter=batch_converter, 
                                         seqthres=seqthres, 
                                         rep_layers=[33], device=device)
        final_rep33.append(np.array(rep33)) 

    res = sequences[['id']].copy()
    res['esm'] = final_rep33
    return res
#endregion

