import json
from collections.abc import Iterable

import numpy as np
import pandas as pd

from rxnrecer.config import config as cfg
from rxnrecer.lib.evaluation import eva


DEFAULT_ZERO_SENTINELS = frozenset({"EC-WITHOUT-REACTION", "NO-PREDICTION"})


def make_label(
    reaction_id,
    rxn_label_dict,
    zero_sentinels: Iterable[str] = DEFAULT_ZERO_SENTINELS,
):
    """Encode reaction ids into a multi-hot label vector.

    This is the single canonical reaction-label encoder for the codebase.
    `reaction_id` may contain multiple RHEA ids separated by `cfg.SPLITER`.
    Values listed in `zero_sentinels` are treated as process sentinels and
    encoded as an all-zero vector.
    """
    res_array = np.zeros(len(rxn_label_dict), dtype=int)

    if reaction_id in set(zero_sentinels):
        return res_array

    for item in reaction_id.split(cfg.SPLITER):
        array_index = rxn_label_dict.get(item)
        if array_index is not None:
            res_array[array_index] = 1
    return res_array


def rxn_eva_metric(eva_df, eva_name, methods, average_type="weighted"):
    results = []
    for method in methods:
        metric = eva.caculateMetrix(
            groundtruth=np.stack(eva_df.lb_rxn_groundtruth),
            predict=np.stack(eva_df[f"lb_rxn_{method}"]),
            baselineName=eva_name,
            type="multi",
            print_flag=False,
            averege_type=average_type,
        )
        results.append(metric)
    return pd.DataFrame(results, columns=["mAccuracy", "mPrecision", "mRecall", "mF1", "avgType"])


def rxn_eva_metric_with_colName(eva_df, col_groundtruth, col_pred, eva_name="", average_type="weighted"):
    groundtruth = np.stack(eva_df[col_groundtruth])
    pred = np.stack(eva_df[col_pred])
    metric = eva.caculateMetrix(
        groundtruth=groundtruth,
        predict=pred,
        baselineName=eva_name,
        type="multi",
        print_flag=False,
        averege_type=average_type,
    )
    return pd.DataFrame([metric], columns=["mAccuracy", "mPrecision", "mRecall", "mF1", "avgType"])


def retrival_reaction_from_ec(ec_pred, ec_reaction_map):
    if ec_pred == "-":
        return "-"
    if ec_pred == "NO-PREDICTION":
        return "NO-PREDICTION"

    if ";" not in ec_pred:
        reaction = ec_reaction_map[ec_reaction_map.ec == ec_pred].reaction_id.to_list()
        if reaction:
            reaction = reaction[0].split(cfg.SPLITER)
    else:
        ecs = [item.strip() for item in ec_pred.split(cfg.SPLITER)]
        reaction = []
        for item_ec in ecs:
            reaction += ec_reaction_map[ec_reaction_map.ec == item_ec].reaction_id.to_list()
        reaction = list({item for group in reaction for item in group.split(cfg.SPLITER)})

    result = cfg.SPLITER.join(sorted(set(reaction)))
    return result if result else "EC-WITHOUT-REACTION"


def load_clean_resluts(res_file):
    data = pd.read_csv(res_file, sep="\t").rename(columns={"Entry": "uniprot_id"})

    def format_ec(eclist):
        ecs = [item.split("/")[0].replace("EC:", "") for item in eclist.split(";") if item is not None]
        return cfg.SPLITER.join(ecs)

    data["ec_clean"] = data.apply(lambda x: format_ec(x.clean_pred_ec_maxsep), axis=1)
    return data[["uniprot_id", "ec_clean"]]


def load_deepec_resluts(filepath):
    result = pd.read_csv(filepath, sep="\t", names=["id", "ec_number"], header=0)
    result.ec_number = result.ec_number.apply(lambda x: x.replace("EC:", ""))
    result.columns = ["id", "ec_deepec"]
    rows = []
    for _, group in result.groupby("id"):
        if len(group) == 1:
            rows.append([group.id.values[0], group.ec_deepec.values[0]])
        else:
            rows.append([group.id.values[0], ";".join(group.ec_deepec.values)])
    return pd.DataFrame(rows, columns=["id", "ec_deepec"])


def load_praim_res(resfile):
    with open(resfile) as handle:
        line = handle.readline()
        counter = 0
        rows = []
        current_id = ""
        subec = []
        while line:
            if ">" in line:
                if counter != 0:
                    rows.append([current_id, ";".join(subec)])
                    subec = []
                current_id = line.replace(">", "").strip()
            elif line.strip():
                ecarray = line.split("\t")
                subec.append(ecarray[0].replace("#", "").strip())
            line = handle.readline()
            counter += 1
    return pd.DataFrame(rows, columns=["id", "ec_priam"])


def load_catfam_res(resfile):
    result = pd.read_csv(resfile, sep="\t", names=["id", "ec_catfam"])
    return result.groupby("id", as_index=False)["ec_catfam"].agg(lambda x: ";".join(x.dropna())).replace("", "-")


def load_ecpred_res(resfile):
    result = pd.read_csv(resfile, sep="\t", header=0)
    result = result.rename(columns={"Protein ID": "id", "EC Number": "ec_ecpred", "Confidence Score(max 1.0)": "pident_ecpred"})
    result["ec_ecpred"] = result.ec_ecpred.apply(lambda x: "-" if x == "non Enzyme" else x)
    return result


def read_h5_file(file_path):
    if file_path.endswith(".feather"):
        return pd.read_feather(file_path)
    if file_path.endswith(".parquet"):
        return pd.read_parquet(file_path)
    if file_path.endswith(".pkl"):
        return pd.read_pickle(file_path)
    with pd.HDFStore(file_path, "r") as store:
        return store["data"]


def get_simi_Pred(pred_list, uniprot_rxn_dict, topk=3):
    uniprot_ids = [item[0] for item in pred_list][:topk]
    rxn_ids = [uniprot_rxn_dict.get(uniprot_id) for uniprot_id in uniprot_ids]
    return cfg.SPLITER.join(set(rxn_ids))


def load_dict_rxn2ec():
    with open(cfg.DICT_RHEA_EC, "r") as handle:
        return json.load(handle)
    
def load_dict_rxn2id():
    with open(cfg.FILE_DS_DICT_RXN2ID, "r") as handel:
        return json.load(handel)
    
def load_dict_id2rxn():
    with open(cfg.FILE_DS_DICT_ID2RXN, "r") as handle:
        return json.load(handle)


def load_dict_ec2rxn():
    with open(cfg.DICT_EC_RHEA, "r") as handle:
        return json.load(handle)


def transRXN2EC(rxns, dict_rxn2ec):
    if rxns == "-":
        return "-"
    ec_list = [dict_rxn2ec.get(rxn, "-") for rxn in rxns.split(";")]
    return f"({' )('.join(ec_list)})".replace(" )(", ")(")


def json_to_dataframe(json_data):
    uniprot_id = list(json_data.keys())[0]
    for entry in json_data[uniprot_id]:
        entry["uniprot_id"] = uniprot_id
    df = pd.DataFrame(json_data[uniprot_id])
    cols = ["uniprot_id"] + [col for col in df.columns if col != "uniprot_id"]
    return df[cols]
