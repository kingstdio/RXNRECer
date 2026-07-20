"""
Dataset splitting and evaluation utilities migrated from legacy rxnrecer.

This module provides functions for:
- Dataset splitting strategies (ta, tb, tc, td)
- Dataset statistics display
- Train/valid/test format conversion
- Enzyme classification evaluation
"""

from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np
import pandas as pd
from pandarallel import pandarallel


def ta_ds_spliter(dataset: pd.DataFrame, test_ratio: float = 0.2) -> pd.DataFrame:
    """
    Split dataset by random sampling (both reactions and proteins can appear in train/test).

    Args:
        dataset: DataFrame with columns including 'reaction_id', 'uniprot_id'
        test_ratio: Ratio of test samples

    Returns:
        DataFrame with 'ds_type' column ('train' or 'test')
    """
    test_size = math.ceil(len(dataset) * test_ratio)

    if test_size * 3 < len(dataset):
        sample_size = test_size * 3
    else:
        sample_size = len(dataset)

    test_df = dataset.sample(n=sample_size)
    train_df = dataset[~dataset.index.isin(test_df.index.values)]

    train_reaction_ids = set(train_df.reaction_id)
    train_uniprot_ids = set(train_df.uniprot_id)

    pandarallel.initialize()
    test_df["ds_type"] = test_df.parallel_apply(
        lambda x: "test"
        if (x.reaction_id in train_reaction_ids) and (x.uniprot_id in train_uniprot_ids)
        else "train",
        axis=1,
    )

    test_df = test_df.sample(n=test_size)

    train = dataset[~dataset.index.isin(test_df.index.values)]
    test = dataset[dataset.index.isin(test_df.index.values)]

    with pd.option_context("mode.chained_assignment", None):
        train["ds_type"] = "train"
        test["ds_type"] = "test"

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df


def tb_ds_spliter(dataset: pd.DataFrame, test_ratio: float = 0.2) -> pd.DataFrame:
    """
    Split dataset by protein (proteins in test set are not in training set).

    Args:
        dataset: DataFrame with columns including 'reaction_id', 'uniprot_id'
        test_ratio: Ratio of test proteins

    Returns:
        DataFrame with 'ds_type' column ('train' or 'test')
    """
    test_size = math.ceil(len(set(dataset.uniprot_id)) * test_ratio)

    if test_size * 3 < len(dataset):
        sample_size = test_size * 3
    else:
        sample_size = len(dataset)

    proteins = dataset["uniprot_id"].drop_duplicates().reset_index(drop=True)
    proteins = proteins.sample(n=sample_size)

    test_df = dataset[dataset.uniprot_id.isin(proteins.values)]
    train_df = dataset[~dataset.uniprot_id.isin(test_df.uniprot_id)]

    train_reaction_ids = list(set(train_df.reaction_id))

    with pd.option_context("mode.chained_assignment", None):
        pandarallel.initialize()
        test_df["ds_type"] = test_df.parallel_apply(
            lambda x: "test" if x.reaction_id.strip() in train_reaction_ids else "train", axis=1
        )

    test_df = test_df[~test_df.uniprot_id.isin(test_df[test_df["ds_type"] == "train"].uniprot_id)]
    proteins_2 = test_df["uniprot_id"].drop_duplicates().reset_index(drop=True)
    proteins_2 = proteins_2.sample(n=test_size)

    test = dataset[dataset.uniprot_id.isin(list(set(proteins_2.values)))]
    train = dataset[~dataset.uniprot_id.isin(list(set(proteins_2.values)))]

    with pd.option_context("mode.chained_assignment", None):
        train["ds_type"] = "train"
        test["ds_type"] = "test"

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df


def tc_ds_spliter(dataset: pd.DataFrame, test_ratio: float = 0.2) -> pd.DataFrame:
    """
    Split dataset by reaction (reactions in test set are not in training set).

    Args:
        dataset: DataFrame with columns including 'reaction_id', 'uniprot_id'
        test_ratio: Ratio of test reactions

    Returns:
        DataFrame with 'ds_type' column ('train' or 'test')
    """
    test_size = math.ceil(len(set(dataset.reaction_id)) * test_ratio)

    if test_size * 3 < len(dataset):
        sample_size = test_size * 3
    else:
        sample_size = len(dataset)

    reactions = dataset["reaction_id"].drop_duplicates().reset_index(drop=True)
    reactions = reactions.sample(n=sample_size)

    test_df = dataset[dataset.reaction_id.isin(reactions.values)]
    train_df = dataset[~dataset.reaction_id.isin(test_df.reaction_id)]

    train_protein_ids = list(set(train_df.uniprot_id))

    with pd.option_context("mode.chained_assignment", None):
        pandarallel.initialize()
        test_df["ds_type"] = test_df.parallel_apply(
            lambda x: "test" if x.uniprot_id in train_protein_ids else "train", axis=1
        )

    test_df = test_df[~test_df.reaction_id.isin(test_df[test_df["ds_type"] == "train"].reaction_id)]
    reactions_2 = test_df["reaction_id"].drop_duplicates().reset_index(drop=True)
    reactions_2 = reactions_2.sample(n=test_size)

    test = dataset[dataset.reaction_id.isin(list(set(reactions_2.values)))]
    train = dataset[~dataset.reaction_id.isin(list(set(reactions_2.values)))]

    with pd.option_context("mode.chained_assignment", None):
        train["ds_type"] = "train"
        test["ds_type"] = "test"

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df


def td_ds_spliter(dataset: pd.DataFrame, test_ratio: float = 0.2, gamma: float = 0.85) -> pd.DataFrame:
    """
    Split dataset with both reaction and protein separation (TDA-like strategy).

    Args:
        dataset: DataFrame with columns including 'reaction_id', 'uniprot_id'
        test_ratio: Ratio of test reactions
        gamma: Ratio of reactions to use for initial test set

    Returns:
        DataFrame with 'ds_type' column ('train' or 'test')
    """
    test_size = math.ceil(len(set(dataset.reaction_id)) * test_ratio)

    reactions = dataset["reaction_id"].drop_duplicates().reset_index(drop=True)
    reactions = reactions.sample(n=math.ceil(test_size * gamma))

    test_df = dataset[dataset.reaction_id.isin(reactions.values)]
    train_df = dataset[~dataset.reaction_id.isin(test_df.reaction_id)]

    train_reaction_ids = set(train_df.reaction_id)
    train_uniprot_ids = set(train_df.uniprot_id)

    test = test_df[
        (~test_df.uniprot_id.isin(list(train_uniprot_ids)))
        & (~test_df.reaction_id.isin(list(train_reaction_ids)))
    ]
    train = dataset[~dataset.index.isin(test.index)]

    with pd.option_context("mode.chained_assignment", None):
        train["ds_type"] = "train"
        test["ds_type"] = "test"

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df


def show_stat_datasets(datasets: pd.DataFrame, setname: str) -> None:
    """
    Display statistics for train/test datasets.

    Args:
        datasets: DataFrame with 'ds_type' column
        setname: Name of the dataset for display
    """
    train = datasets[datasets.ds_type == "train"]
    test = datasets[datasets.ds_type == "test"]

    protein_in_train = len(set(train.uniprot_id))
    reaction_in_train = len(set(train.reaction_id))
    protein_in_test = len(set(test.uniprot_id))
    reaction_in_test = len(set(test.reaction_id))

    print(60 * "-")
    print(setname)
    print(60 * "-")
    print(f"Trainning records:             {len(train)}")
    print(f"Proteins in trainning set:     {protein_in_train}")
    print(f"Reactions in trainning set:    {reaction_in_train}")
    print(f"Testing records:               {len(test)}")
    print(f"Proteins in testing set:       {protein_in_test}")
    print(f"Reactions in testing set:      {reaction_in_test}")
    print("")


def format_train_valid_test(
    dataset: pd.DataFrame,
    feature_reactions: pd.DataFrame,
    feature_proteins: pd.DataFrame,
    vali_ratio: float = 0.2,
) -> Tuple[np.ndarray, List, np.ndarray, List, np.ndarray, List]:
    """
    Format dataset into train/valid/test arrays for model training.

    Args:
        dataset: DataFrame with 'reaction_id', 'uniprot_id', 'label', 'ds_type'
        feature_reactions: DataFrame with reaction features
        feature_proteins: DataFrame with protein features
        vali_ratio: Ratio of training data to use for validation

    Returns:
        Tuple of (train_X, train_Y, vali_X, vali_Y, test_X, test_Y)
    """
    dataset = dataset[["reaction_id", "uniprot_id", "label", "ds_type"]].reset_index(drop=True)
    dataset = dataset.merge(
        feature_reactions, on="reaction_id", how="left"
    ).merge(feature_proteins.rename(columns={"id": "uniprot_id"}), on="uniprot_id", how="left")

    train = dataset[dataset.ds_type == "train"].reset_index(drop=True)
    test = dataset[dataset.ds_type == "test"].reset_index(drop=True)

    vali = train.sample(frac=vali_ratio)
    train = train[~train.index.isin(vali.index)]

    vali = vali.reset_index(drop=True)
    train = train.reset_index(drop=True)

    train_X = train.iloc[:, 4:].values
    train_Y = train.label.tolist()

    vali_X = vali.iloc[:, 4:].values
    vali_Y = vali.label.tolist()

    test_X = test.iloc[:, 4:].values
    test_Y = test.label.tolist()

    return train_X, train_Y, vali_X, vali_Y, test_X, test_Y


def eva_isenzyme(
    baseline_name: str, res_df: pd.DataFrame, category: str, print_flag: bool = False
) -> List:
    """
    Evaluate enzyme/non-enzyme classification.

    Args:
        baseline_name: Name of the baseline method
        res_df: DataFrame with prediction and ground truth columns
        category: 'ec' or 'rxn' for different evaluation modes
        print_flag: Whether to print results

    Returns:
        List containing [baselineName, accuracy, precision, recall, ppv, npv, f1, tp, fp, fn, tn, up, un]
    """
    import re

    if category == "ec":
        tp = res_df[
            (res_df[f"{category}_{baseline_name}"].str.contains(r"\d", na=False))
            & (res_df.isenzyme_groundtruth == True)
        ].shape[0]
        tn = res_df[
            (res_df[f"{category}_{baseline_name}"] == "-") & (res_df.isenzyme_groundtruth == False)
        ].shape[0]
        fp = res_df[
            (res_df[f"{category}_{baseline_name}"].str.contains(r"\d", na=False))
            & (res_df.isenzyme_groundtruth == False)
        ].shape[0]
        fn = res_df[
            (res_df[f"{category}_{baseline_name}"] == "-") & (res_df.isenzyme_groundtruth == True)
        ].shape[0]
    elif category == "rxn":
        tp = res_df[
            (res_df[f"{category}_{baseline_name}"].str.contains("RHEA", na=False))
            & (res_df.rxn_groundtruth != "-")
        ].shape[0]
        tn = res_df[
            (res_df[f"{category}_{baseline_name}"] == "-") & (res_df.rxn_groundtruth == "-")
        ].shape[0]
        fp = res_df[
            (res_df[f"{category}_{baseline_name}"].str.contains("RHEA", na=False))
            & (res_df.rxn_groundtruth == "-")
        ].shape[0]
        fn = res_df[
            (res_df[f"{category}_{baseline_name}"] == "-") & (res_df.rxn_groundtruth.str.contains("RHEA"))
        ].shape[0]
    else:
        raise ValueError("Invalid category. Please choose 'ec' or 'rxn'.")

    up = res_df[
        (res_df[f"{category}_{baseline_name}"] == "NO-PREDICTION")
        & (
            res_df.isenzyme_groundtruth == True
            if category == "ec"
            else res_df.rxn_groundtruth.str.contains("RHEA")
        )
    ].shape[0]
    un = res_df[
        (res_df[f"{category}_{baseline_name}"] == "NO-PREDICTION")
        & (
            res_df.isenzyme_groundtruth == False
            if category == "ec"
            else res_df.rxn_groundtruth == "-"
        )
    ].shape[0]

    acc = (tp + tn) / (tp + tn + fp + fn + up + un + 1.4e-45)
    precision = tp / (tp + fp + up + un + 1.4e-45)
    recall = tp / (tp + fn + up + un + 1.4e-45)
    f1 = 2 * precision * recall / (precision + recall + 1.4e-45)
    ppv = tp / (tp + fp + up + un + 1.4e-45)
    npv = tn / (tn + fn + un + un + 1.4e-45)

    if print_flag:
        print(
            "{:<20} {:<14} {:<14} {:<15} {:<15} {:<20} {:<15} {:<10} {:<10} {:<10}{:<10}{:<10}{}".format(
                "baselineName",
                "accuracy",
                "precision",
                "recall",
                "PPV(Sensitivity)",
                "NPV(Specificity)",
                "F1",
                "TP",
                "FP",
                "FN",
                "TN",
                "UP",
                "UN",
            )
        )
        print(
            "{:<20} {:<14.6f} {:<14.6f} {:<15.6f} {:<15.6f} {:<20.6f} {:<15.6f} {:<10} {:<10} {:<10}{:<10}{:<10}{}".format(
                baseline_name, acc, precision, recall, ppv, npv, f1, tp, fp, fn, tn, up, un
            )
        )

    return [baseline_name, acc, precision, recall, ppv, npv, f1, tp, fp, fn, tn, up, un]
