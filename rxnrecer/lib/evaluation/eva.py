"""
Unified Evaluation Module for RXNRECer

This module consolidates evaluation functions from:
- rxnrecer/lib/evaluation/eva.py
- tools/evaluation.py
- rxnrecer/utils/evaluation_utils.py
"""

from __future__ import annotations

import itertools
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn import metrics
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    classification_report,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from tqdm import tqdm

try:
    from IPython.display import HTML
except ImportError:  # pragma: no cover - optional notebook dependency
    HTML = None

try:
    import plotly.graph_objects as go
except ImportError:  # pragma: no cover - optional notebook dependency
    go = None


def _effective_cpu_count() -> int:
    """返回当前进程真正可用的 CPU 数。"""
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus:
        try:
            return max(1, int(slurm_cpus))
        except ValueError:
            pass

    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        pass

    return max(1, os.cpu_count() or 1)


# =============================================================================
# Core Evaluation Functions (from rxnrecer/lib/evaluation/eva.py)
# =============================================================================


def read_10fold_res_csv_files(file_paths):
    """Read 10-fold cross-validation result CSV files."""
    return [pd.read_csv(file, sep="\t") for file in file_paths]


def process_no_res(res_list, eckey, rxnkey):
    """Process results without predictions."""
    summary = [
        [
            len(pred_detail),
            len(pred_detail[pred_detail[eckey].str.contains("NO-PREDICTION")]),
            len(pred_detail[pred_detail[rxnkey].str.contains("EC-WITHOUT-REACTION")]),
        ]
        for pred_detail in res_list
    ]
    return pd.DataFrame(summary, columns=["test_size", "no_prediction", "ec_without_rxn"])


def calculate_metrics_multi_joblib(
    groundtruth, predict, average_type, print_flag=False, n_jobs=4
):
    """Calculate multiple metrics in parallel using joblib."""
    is_multilabel = len(groundtruth.shape) > 1 and groundtruth.shape[1] > 1
    if average_type == "samples" and not is_multilabel:
        raise ValueError(
            "Samplewise metrics are not available outside of multilabel classification."
        )

    acc = metrics.accuracy_score(groundtruth, predict)
    p, r, f, _ = metrics.precision_recall_fscore_support(
        groundtruth, predict, average=average_type, zero_division=0
    )
    
    results = [acc, p, r, f]

    if print_flag:
        print(
            f"{results[0]:.6f}\t{results[1]:.6f}\t{results[2]:.6f}\t{results[3]:.6f}\t{average_type:>12s}"
        )

    return results + [average_type]


def caculateMetrix(
    groundtruth, predict, baselineName="", type="binary", averege_type="macro", print_flag=True
):
    """Legacy metric entrypoint retained for compatibility."""
    if type == "multi":
        return calculate_metrics_multi_joblib(
            groundtruth=groundtruth,
            predict=predict,
            average_type=averege_type,
            print_flag=print_flag,
            n_jobs=4,
        )

    if type == "binary":
        acc = metrics.accuracy_score(groundtruth, predict)
        precision = metrics.precision_score(groundtruth, predict, zero_division=0)
        recall = metrics.recall_score(groundtruth, predict, zero_division=0)
        f1 = metrics.f1_score(groundtruth, predict, zero_division=0)
        return [baselineName, acc, precision, recall, f1]

    raise ValueError(f"Unsupported type: {type}")


def eva_one_fold(
    eva_df,
    lb_groundtruth,
    lb_predict,
    fold_num=None,
    n_jobs=4,
    avg_types=None,
):
    """Evaluate one fold of cross-validation."""
    groundtruth = np.stack(eva_df[lb_groundtruth])
    predict = np.stack(eva_df[lb_predict])
    is_multilabel = len(groundtruth.shape) > 1 and groundtruth.shape[1] > 1

    if avg_types is None:
        average_types = ["weighted", "micro", "macro"]
        if is_multilabel:
            average_types.append("samples")
    else:
        average_types = list(avg_types)
        if "samples" in average_types and not is_multilabel:
            raise ValueError(
                "Samplewise metrics are not available outside of multilabel classification."
            )

    # Accuracy only needs to be calculated once
    acc = metrics.accuracy_score(groundtruth, predict)

    def calc_metrics_for_avg_type(avg_type):
        p, r, f, _ = metrics.precision_recall_fscore_support(
            groundtruth, predict, average=avg_type, zero_division=0
        )
        return [acc, p, r, f, avg_type]

    if int(n_jobs) == 1:
        results = [calc_metrics_for_avg_type(avg_type) for avg_type in average_types]
    else:
        results = Parallel(n_jobs=n_jobs)(
            delayed(calc_metrics_for_avg_type)(avg_type) for avg_type in average_types
        )

    res = pd.DataFrame(results, columns=["mAccuracy", "mPrecision", "mRecall", "mF1", "avgType"])
    if fold_num is not None:
        res.insert(0, "runFold", fold_num)
    return res


def eva_cross_validation(
    res_df,
    lb_groundtruth,
    lb_predict,
    num_folds=10,
    max_workers=None,
    metric_n_jobs=None,
    avg_types=None,
):
    """Evaluate cross-validation results.

    ``max_workers`` parallelizes across folds. By default, all folds are
    evaluated concurrently, capped by CPU count. Set ``max_workers=1`` to force
    serial execution. When fold-level parallelism is enabled, per-fold metric
    jobs default to 1 to avoid nested oversubscription.
    """
    if max_workers is None:
        max_workers = min(num_folds, _effective_cpu_count())
    if metric_n_jobs is None:
        metric_n_jobs = 1 if max_workers and max_workers > 1 else 4

    def evaluate_fold(runfold):
        return eva_one_fold(
            eva_df=res_df[res_df.run_fold == runfold].reset_index(drop=True),
            lb_groundtruth=lb_groundtruth,
            lb_predict=lb_predict,
            fold_num=runfold,
            n_jobs=metric_n_jobs,
            avg_types=avg_types,
        )

    runfolds = list(range(1, num_folds + 1))
    if not max_workers or max_workers <= 1:
        eva_metrics = [evaluate_fold(runfold) for runfold in tqdm(runfolds)]
        return pd.concat(eva_metrics, axis=0).reset_index(drop=True)

    eva_metrics = [None] * len(runfolds)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(evaluate_fold, runfold): index
            for index, runfold in enumerate(runfolds)
        }
        for future in tqdm(as_completed(futures), total=len(futures)):
            eva_metrics[futures[future]] = future.result()
    return pd.concat(eva_metrics, axis=0).reset_index(drop=True)


def calculate_classification_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: Optional[np.ndarray] = None,
    average: str = "weighted",
) -> Dict[str, float]:
    """Calculate generic classification metrics."""
    metrics_result = {
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, average=average, zero_division=0),
        "recall": recall_score(y_true, y_pred, average=average, zero_division=0),
        "f1": f1_score(y_true, y_pred, average=average, zero_division=0),
    }

    if y_prob is not None:
        try:
            if len(y_prob.shape) == 1:
                metrics_result["roc_auc"] = roc_auc_score(y_true, y_prob)
                metrics_result["pr_auc"] = average_precision_score(y_true, y_prob)
            else:
                metrics_result["roc_auc"] = roc_auc_score(
                    y_true, y_prob, multi_class="ovr", average=average
                )
                metrics_result["pr_auc"] = average_precision_score(y_true, y_prob, average=average)
        except Exception as e:
            print(f"Warning: Could not calculate probability-based metrics: {e}")

    return metrics_result


def calculate_eval_dataframe_metrics(
    eva_df: pd.DataFrame,
    ground_truth_col: str,
    pred_col: str,
    metric_fn: Optional[Callable] = None,
    eva_name: str = "",
    avg_method: str = "weighted",
) -> pd.DataFrame:
    """Calculate legacy evaluation metrics for label-array DataFrame columns."""
    if metric_fn is not None:
        res = metric_fn(
            eva_df=eva_df,
            col_groundtruth=ground_truth_col,
            col_pred=pred_col,
            eva_name=eva_name,
            average_type=avg_method,
        )
        if "evaName" not in res.columns:
            res.insert(0, "evaName", eva_name)
        return res

    groundtruth = np.stack(eva_df[ground_truth_col])
    predict = np.stack(eva_df[pred_col])
    metric = caculateMetrix(
        groundtruth=groundtruth,
        predict=predict,
        baselineName=eva_name,
        type="multi",
        averege_type=avg_method,
        print_flag=False,
    )
    res = pd.DataFrame(
        [metric],
        columns=["mAccuracy", "mPrecision", "mRecall", "mF1", "avgType"],
    )
    res.insert(0, "evaName", eva_name)
    return res


def calculate_metrics(*args, **kwargs):
    """Compatibility dispatcher for legacy DataFrame metrics and modern array metrics.

    Legacy usage:
        calculate_metrics(eva_df, ground_truth_col, pred_col, metric_fn, eva_name, avg_method)
        calculate_metrics(eva_df=..., ground_truth_col=..., pred_col=...)

    Modern usage:
        calculate_metrics(y_true, y_pred, y_prob=None, average="weighted")
        calculate_metrics(y_true=..., y_pred=...)
    """
    if args and isinstance(args[0], pd.DataFrame):
        eva_df = args[0]
        ground_truth_col = args[1] if len(args) > 1 else kwargs.pop("ground_truth_col")
        pred_col = args[2] if len(args) > 2 else kwargs.pop("pred_col")
        metric_fn = args[3] if len(args) > 3 else kwargs.pop("metric_fn", None)
        eva_name = args[4] if len(args) > 4 else kwargs.pop("eva_name", "")
        avg_method = args[5] if len(args) > 5 else kwargs.pop("avg_method", "weighted")
        return calculate_eval_dataframe_metrics(
            eva_df=eva_df,
            ground_truth_col=ground_truth_col,
            pred_col=pred_col,
            metric_fn=metric_fn,
            eva_name=eva_name,
            avg_method=avg_method,
        )

    if "eva_df" in kwargs:
        return calculate_eval_dataframe_metrics(
            eva_df=kwargs.pop("eva_df"),
            ground_truth_col=kwargs.pop("ground_truth_col"),
            pred_col=kwargs.pop("pred_col"),
            metric_fn=kwargs.pop("metric_fn", None),
            eva_name=kwargs.pop("eva_name", ""),
            avg_method=kwargs.pop("avg_method", "weighted"),
        )

    y_true = args[0] if args else kwargs.pop("y_true")
    y_pred = args[1] if len(args) > 1 else kwargs.pop("y_pred")
    y_prob = args[2] if len(args) > 2 else kwargs.pop("y_prob", None)
    average = kwargs.pop("average", kwargs.pop("avg_method", "weighted"))
    return calculate_classification_metrics(y_true, y_pred, y_prob=y_prob, average=average)


def calculate_metrics_parallel(
    res_df,
    ground_truth_col,
    pred_col,
    metric_fn: Optional[Callable] = None,
    avg_method="weighted",
    max_workers=None,
):
    """Calculate metrics in parallel across folds."""

    def run_metric_evaluation(index):
        return calculate_metrics(
            eva_df=res_df[index],
            ground_truth_col=ground_truth_col,
            pred_col=pred_col,
            metric_fn=metric_fn,
            eva_name=f"fold{index + 1}",
            avg_method=avg_method,
        )

    results = [None] * len(res_df)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_metric_evaluation, i): i for i in range(len(res_df))}
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    return pd.concat(results, axis=0).reset_index(drop=True)


def get_fold_mean_std_metrics(input_df):
    """Calculate mean and std across folds."""
    res_fold_std = input_df[
        ["baselineName", "avgType", "mAccuracy", "mPrecision", "mRecall", "mF1"]
    ]
    res_fold_std = res_fold_std.groupby(["baselineName", "avgType"]).agg(["mean", "std"]).reset_index()
    res_fold_std.columns = ["_".join(filter(None, col)).strip() for col in res_fold_std.columns]
    melted = res_fold_std.melt(
        id_vars=["baselineName", "avgType"], var_name="Metric_Statistic", value_name="Value"
    )
    melted[["Metric", "Statistic"]] = melted["Metric_Statistic"].str.rsplit("_", n=1, expand=True)
    melted = melted.sort_values(by=["baselineName", "Metric", "avgType"]).reset_index(drop=True)
    melted = melted.drop(columns=["Metric_Statistic"])
    pivot = melted.pivot_table(
        index=["baselineName", "avgType", "Metric"], columns="Statistic", values="Value"
    ).reset_index()
    pivot.columns.name = None
    return pivot


def statistic_no_res(res_df, name_col_ec, name_col_rxn, type="ec"):
    """Statistic no prediction results."""
    if type == "ec":
        return res_df.groupby("run_fold").agg(
            test_size=("run_fold", "count"),
            no_prediction_count=(name_col_ec, lambda x: (x == "NO-PREDICTION").sum()),
            ec_without_reaction_count=(
                name_col_rxn,
                lambda x: (x == "EC-WITHOUT-REACTION").sum(),
            ),
        ).reset_index()

    return res_df.groupby("run_fold").agg(
        test_size=("run_fold", "count"),
        no_prediction_count=(name_col_rxn, lambda x: (x == "NO-PREDICTION").sum()),
    ).reset_index()


def display_html_results(metrics_df=None, std_mean=None, no_pred=None, eva_name="", **kwargs):
    """Display evaluation results as HTML."""
    if metrics_df is None:
        metrics_df = kwargs.pop("metrics", None)
    if kwargs:
        raise TypeError(f"Unexpected keyword arguments: {', '.join(kwargs)}")
    if metrics_df is None:
        raise ValueError("metrics_df is required")
    if std_mean is None:
        raise ValueError("std_mean is required")
    if HTML is None:
        raise RuntimeError("IPython is required for HTML rendering")
    return HTML(
        f"""
         <div style="float:left; width:auto; margin-right:30px;">
              <h2 style='color:blue'>{eva_name} Evaluation 10 Fold Details</h2>
              {metrics_df.to_html()}
         </div>
         <div  style="float:left; width:auto; margin-right:30px;" >
              <h2 style='color:blue' >{eva_name} Evaluation 10 Fold std and mean Overview</h2>
                   {std_mean.to_html()}
         </div>
         <div style="float:left; width:auto;">
         <h2 style='color:blue' >{eva_name} Evaluation 10 Fold No Prediction Overview</h2>
              {no_pred.to_html() if no_pred is not None else ''}
         </div>
         """
    )


# =============================================================================
# Extended Evaluation Functions (from tools/evaluation.py)
# =============================================================================


def calculate_single_metric(method, avg_type, df, ground_truth_col):
    """Calculate single method metric."""
    try:
        res = calculate_metrics(
            eva_df=df,
            ground_truth_col=ground_truth_col,
            pred_col=method,
            eva_name=method,
            avg_method=avg_type,
        )
        return res
    except Exception as e:
        print(f"Error calculating {method} {avg_type}: {e}")
        return None


def get_metrics(pred_cols, avg_types, df, ground_truth_col="rxn_groundtruth", n_jobs=20):
    """Get metrics for multiple prediction columns in parallel."""
    tasks = []
    for method in pred_cols:
        for avg_type in avg_types:
            tasks.append((method, avg_type, df, ground_truth_col))

    results = Parallel(n_jobs=n_jobs, backend="threading")(
        delayed(calculate_single_metric)(method, avg_type, df, ground_truth_col)
        for method, avg_type, df, ground_truth_col in tasks
    )

    valid_results = [r for r in results if r is not None]
    if valid_results:
        return pd.concat(valid_results, ignore_index=True)
    return pd.DataFrame()


# =============================================================================
# Ensemble Methods
# =============================================================================


def recall_boosted_ensemble(res_array):
    """Recall-boosted ensemble (union of predictions).

    Rule: if any position is 1, result is 1; only all 0 gives 0.
    """
    return np.maximum.reduce(res_array).astype(int)


def majority_vote_with_priority(arrays, priority_array, s1first=False):
    """Majority vote ensemble with priority array.

    Rule: priority_array=1 gives 1; otherwise majority vote (tie prefers 1).
    """
    result = priority_array.copy().astype(int)

    if s1first:
        zero_mask = priority_array == 0
    else:
        zero_mask = np.ones_like(priority_array, dtype=bool)

    if np.any(zero_mask):
        stacked = np.stack(arrays, axis=0)
        ones_count = np.sum(stacked, axis=0)
        majority = ones_count >= len(arrays) / 2
        result[zero_mask] = majority[zero_mask].astype(int)

    return result


def get_combinations_without_order(items, r=None):
    """Generate unordered combinations of items."""
    n = len(items)
    if n == 0:
        return []

    if r is None:
        sizes = range(2, n + 1)
    elif isinstance(r, int):
        sizes = [r]
    else:
        sizes = []
        for k in r:
            try:
                k = int(k)
                if 1 <= k <= n:
                    sizes.append(k)
            except Exception:
                continue
        sizes = list(dict.fromkeys(sizes))
        if not sizes:
            return []

    res = []
    for k in sizes:
        if 1 <= k <= n:
            res.extend(itertools.combinations(items, k))

    return [list(item) for item in res]


def calculate_ensemble_metrics(
    cb_methods,
    baseline_method="ESMwithCLF",
    df=None,
    ground_truth_col="rxn_groundtruth",
    avg_types=["weighted"],
    combination_sizes=[2, 3, 4],
    ensemble_types=["recall_boosted", "majority"],
):
    """Calculate ensemble method metrics."""
    all_metrics = []

    for size in combination_sizes:
        combinations = get_combinations_without_order(cb_methods, r=size)
        print(f"Computing {size + 1} method combinations, total {len(combinations)}")

        for ensemble_type in ensemble_types:
            print(f"  Computing {ensemble_type} ensemble...")

            col_names = []
            for combo in combinations:
                col_name = "_".join(combo) + f"_{baseline_method}_{ensemble_type}"
                col_names.append(col_name)

            for i, combo in enumerate(combinations):
                if ensemble_type == "recall_boosted":
                    df[col_names[i]] = df.apply(
                        lambda x: recall_boosted_ensemble(
                            [x[method] for method in combo] + [x[baseline_method]]
                        ),
                        axis=1,
                    )
                elif ensemble_type == "majority":
                    df[col_names[i]] = df.apply(
                        lambda x: majority_vote_with_priority(
                            [x[method] for method in combo], x[baseline_method], s1first=True
                        ),
                        axis=1,
                    )

            metrics = get_metrics(
                pred_cols=col_names,
                avg_types=avg_types,
                df=df,
                ground_truth_col=ground_truth_col,
            )
            all_metrics.append(metrics)

    if all_metrics:
        return pd.concat(all_metrics, ignore_index=True)
    return pd.DataFrame()


# =============================================================================
# Visualization Functions
# =============================================================================


def show_ec_methods_10_eva_fig(res_metrics_data):
    """Display EC methods 10-fold evaluation figure."""
    if go is None:
        raise RuntimeError("plotly is required for visualization")

    colors = [
        "#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC",
        "#E7DAD2", "#999999", "#A1D3B2", "#F5C98A", "#F9988C",
        "#6FB7AA", "#FFA351", "#E45756", "#5A9BD4", "#9D81BA", "#C9C1B6",
    ]

    methods = (
        res_metrics_data[res_metrics_data.Metric == "mF1"]
        .sort_values(by=["mean"])
        .baselineName.tolist()
    )

    if "blast_via_ec" in methods and "blast_via_rxn" in methods:
        methods = [
            "priam", "deepec", "clean", "catfam", "blast_via_ec",
            "ecrecer", "blast_via_rxn", "unirep_euclidean", "unirep_cosine",
            "esm_euclidean", "esm_cosine", "t5_euclidean", "t5_cosine", "RXNRECer",
        ]
        bar_width = 0.04
        title_text = "Performance Comparison of All Reaction Prediction Methods"
    elif "RXNRECer" in methods and "blast_via_ec" not in methods:
        methods = [
            "blast", "unirep_euclidean", "unirep_cosine", "esm_euclidean",
            "esm_cosine", "t5_euclidean", "t5_cosine", "RXNRECer",
        ]
        bar_width = 0.08
        title_text = "Performance Comparison of Direct Reaction Prediction Methods"
    else:
        bar_width = 0.11
        title_text = "Performance Comparison of EC-based Reaction Prediction Methods"

    fig = go.Figure()
    showMetics = ["mAccuracy", "mPrecision", "mRecall", "mF1"]

    for idx, method in enumerate(methods):
        df_method = res_metrics_data[res_metrics_data["baselineName"] == method].copy()
        df_method["Metric"] = pd.Categorical(
            df_method["Metric"], categories=showMetics, ordered=True
        )
        df_method = df_method.sort_values(by="Metric").reset_index(drop=True)

        if not df_method.empty:
            fig.add_trace(
                go.Bar(
                    name=method,
                    x=showMetics,
                    y=df_method["mean"],
                    width=bar_width,
                    error_y=dict(type="data", array=df_method["std"], visible=True),
                    marker_color=colors[idx % len(colors)],
                    marker_line=dict(color="black", width=0.5),
                )
            )

    fig.update_layout(
        yaxis=dict(
            showline=True,
            linecolor="black",
            linewidth=1,
            showgrid=True,
            gridcolor="gray",
            gridwidth=1,
            minor=dict(showgrid=True, griddash="dash", gridcolor="gray", gridwidth=0.5, dtick=0.05),
            dtick=0.1,
            range=[0, 1.2],
        ),
        xaxis=dict(showline=True, linecolor="#000000", linewidth=1),
        title=dict(
            text=title_text,
            x=0.5,
            y=0.05,
            xanchor="center",
            font=dict(size=20, weight=1000),
        ),
        yaxis_title="Mean Value",
        width=1600,
        height=600,
        barmode="group",
        bargap=0.29,
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=0.99,
            xanchor="center",
            x=0.5,
            tracegroupgap=2,
            borderwidth=1,
        ),
    )

    return fig


def show_methods_eva_fig(res_metrics_data):
    """Display methods evaluation figure without error bars."""
    if go is None:
        raise RuntimeError("plotly is required for visualization")

    colors = [
        "#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2",
        "#BEB8DC", "#E7DAD2", "#999999", "#A1D3B2",
        "#F5C98A", "#F9988C",
    ]

    methods_sorted = (
        res_metrics_data.sort_values(by="mF1", ascending=False)["baselineName"].tolist()
    )

    unique_groups = res_metrics_data["methodGroup"].unique()
    bar_width = 0.08
    title_text = ""

    if len(unique_groups) > 1:
        direct_order = [
            "MSA-via-RXN", "Unirep-cosine", "Unirep-euclidean", "ESM-cosine",
            "ESM-euclidean", "T5-cosine", "T5-euclidean", "RXNRECer-S1",
        ]
        ec_order = ["PRIAM", "DeepEC", "CLEAN", "CatFam", "MSA-via-EC", "ECRECer"]
        methods_order = []
        for m in direct_order:
            if m in res_metrics_data["baselineName"].values:
                methods_order.append(m)
        for m in ec_order:
            if m in res_metrics_data["baselineName"].values:
                methods_order.append(m)
        methods = methods_order
        bar_width = 0.04
        title_text = "Performance Comparison of All Reaction Prediction Methods"
    elif unique_groups[0] == "direct":
        methods = methods_sorted
        bar_width = 0.08
        title_text = "Performance Comparison of Direct Reaction Prediction Methods"
    else:
        methods = methods_sorted
        bar_width = 0.11
        title_text = "Performance Comparison of EC-based Reaction Prediction Methods"

    metrics = ["mAccuracy", "mPrecision", "mRecall", "mF1"]
    fig = go.Figure()

    for idx, method in enumerate(methods):
        df_method = res_metrics_data[res_metrics_data["baselineName"] == method]
        if not df_method.empty:
            row = df_method.iloc[0]
            y_values = [row[metric] for metric in metrics]
            fig.add_trace(
                go.Bar(
                    name=method,
                    x=metrics,
                    y=y_values,
                    width=bar_width,
                    marker_color=colors[idx % len(colors)],
                    marker_line=dict(color="black", width=0.5),
                )
            )

    fig.update_layout(
        yaxis=dict(
            showline=True,
            linecolor="black",
            linewidth=1,
            showgrid=True,
            gridcolor="gray",
            gridwidth=1,
            minor=dict(showgrid=True, griddash="dash", gridcolor="gray", gridwidth=0.5, dtick=0.05),
            dtick=0.1,
            range=[0, 1.2],
        ),
        xaxis=dict(showline=True, linecolor="#000000", linewidth=1),
        title=dict(
            text=title_text,
            x=0.5,
            y=0.05,
            xanchor="center",
            font=dict(size=20, weight=1000),
        ),
        yaxis_title="Mean Value",
        width=1600,
        height=600,
        barmode="group",
        bargap=0.29,
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            tracegroupgap=2,
            borderwidth=1,
        ),
    )

    return fig


# =============================================================================
# Modern Evaluation Utilities (from rxnrecer/utils/evaluation_utils.py)
# =============================================================================


def calculate_metrics_new(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: Optional[np.ndarray] = None,
    average: str = "weighted",
) -> Dict[str, float]:
    """Calculate comprehensive evaluation metrics.

    Args:
        y_true: True labels
        y_pred: Predicted labels
        y_prob: Predicted probabilities (optional)
        average: Averaging method for multi-class metrics

    Returns:
        Dictionary with evaluation metrics
    """
    metrics_result = {}

    metrics_result["accuracy"] = accuracy_score(y_true, y_pred)
    metrics_result["precision"] = precision_score(y_true, y_pred, average=average, zero_division=0)
    metrics_result["recall"] = recall_score(y_true, y_pred, average=average, zero_division=0)
    metrics_result["f1"] = metrics.f1_score(y_true, y_pred, average=average, zero_division=0)

    if y_prob is not None:
        try:
            if len(y_prob.shape) == 1:
                metrics_result["roc_auc"] = roc_auc_score(y_true, y_prob)
                metrics_result["pr_auc"] = average_precision_score(y_true, y_prob)
            else:
                metrics_result["roc_auc"] = roc_auc_score(
                    y_true, y_prob, multi_class="ovr", average=average
                )
                metrics_result["pr_auc"] = average_precision_score(y_true, y_prob, average=average)
        except Exception as e:
            print(f"Warning: Could not calculate probability-based metrics: {e}")

    return metrics_result


def calculate_reaction_metrics(
    ground_truth: pd.DataFrame,
    predictions: pd.DataFrame,
    reaction_col: str = "reaction_id",
    top_k: int = 5,
) -> Dict[str, float]:
    """Calculate reaction prediction specific metrics."""
    metrics_result = {}

    common_idx = ground_truth.index.intersection(predictions.index)
    gt = ground_truth.loc[common_idx]
    pred = predictions.loc[common_idx]

    if f"top_{top_k}_predictions" in pred.columns:
        top_k_correct = 0
        for idx in common_idx:
            true_rxn = gt.loc[idx, reaction_col]
            pred_rxns = pred.loc[idx, f"top_{top_k}_predictions"]
            if isinstance(pred_rxns, str):
                pred_rxns = pred_rxns.split(";")
            if true_rxn in pred_rxns:
                top_k_correct += 1

        metrics_result[f"top_{top_k}_accuracy"] = top_k_correct / len(common_idx)

    if "prediction_rank" in pred.columns:
        reciprocal_ranks = []
        for idx in common_idx:
            rank = pred.loc[idx, "prediction_rank"]
            if rank > 0:
                reciprocal_ranks.append(1.0 / rank)
            else:
                reciprocal_ranks.append(0.0)

        metrics_result["mrr"] = np.mean(reciprocal_ranks)

    return metrics_result


def split_dataset(
    dataset: pd.DataFrame,
    test_ratio: float = 0.2,
    split_method: str = "random",
    stratify_col: Optional[str] = None,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split dataset into train and test sets."""
    if split_method == "random":
        if stratify_col:
            train_df, test_df = train_test_split(
                dataset,
                test_size=test_ratio,
                stratify=dataset[stratify_col],
                random_state=random_state,
            )
        else:
            train_df, test_df = train_test_split(
                dataset, test_size=test_ratio, random_state=random_state
            )
    elif split_method == "protein":
        proteins = dataset["uniprot_id"].unique()
        test_proteins = np.random.choice(
            proteins, size=int(len(proteins) * test_ratio), replace=False
        )
        test_df = dataset[dataset["uniprot_id"].isin(test_proteins)]
        train_df = dataset[~dataset["uniprot_id"].isin(test_proteins)]
    elif split_method == "reaction":
        reactions = dataset["reaction_id"].unique()
        test_reactions = np.random.choice(
            reactions, size=int(len(reactions) * test_ratio), replace=False
        )
        test_df = dataset[dataset["reaction_id"].isin(test_reactions)]
        train_df = dataset[~dataset["reaction_id"].isin(test_reactions)]
    else:
        raise ValueError(f"Unknown split method: {split_method}")

    train_df = train_df.copy()
    test_df = test_df.copy()
    train_df["split"] = "train"
    test_df["split"] = "test"

    return train_df, test_df


def calculate_confidence_intervals(
    metrics_dict: Dict[str, List[float]], confidence_level: float = 0.95
) -> Dict[str, Dict[str, float]]:
    """Calculate confidence intervals for metrics."""
    intervals = {}

    for metric_name, values in metrics_dict.items():
        if len(values) < 2:
            continue

        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)
        n = len(values)

        t_critical = 1.96

        margin_of_error = t_critical * (std_val / np.sqrt(n))

        intervals[metric_name] = {
            "mean": mean_val,
            "std": std_val,
            "lower": mean_val - margin_of_error,
            "upper": mean_val + margin_of_error,
            "n": n,
        }

    return intervals


def format_classification_report(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    target_names: Optional[List[str]] = None,
) -> str:
    """Format classification report."""
    return classification_report(y_true, y_pred, target_names=target_names)


def calculate_hit_rate(
    predictions: pd.DataFrame,
    ground_truth: pd.DataFrame,
    top_k: int = 1,
    reaction_col: str = "reaction_id",
) -> float:
    """Calculate hit rate (percentage of correct predictions in top-k)."""
    hits = 0
    total = 0

    for idx in ground_truth.index:
        if idx not in predictions.index:
            continue

        true_rxn = ground_truth.loc[idx, reaction_col]

        if f"top_{top_k}_predictions" in predictions.columns:
            pred_rxns = predictions.loc[idx, f"top_{top_k}_predictions"]
            if isinstance(pred_rxns, str):
                pred_rxns = pred_rxns.split(";")

            if true_rxn in pred_rxns[:top_k]:
                hits += 1
            total += 1

    return hits / total if total > 0 else 0.0


def calculate_coverage(
    predictions: pd.DataFrame, unique_reactions: Optional[List[str]] = None
) -> float:
    """Calculate prediction coverage."""
    if unique_reactions is None:
        all_pred_rxns = []
        for col in predictions.columns:
            if col.startswith("top_") and col.endswith("_predictions"):
                pred_rxns = predictions[col].dropna()
                for pred_str in pred_rxns:
                    if isinstance(pred_str, str):
                        all_pred_rxns.extend(pred_str.split(";"))
        unique_reactions = list(set(all_pred_rxns))

    predicted_reactions = set()
    for col in predictions.columns:
        if col.startswith("top_") and col.endswith("_predictions"):
            pred_rxns = predictions[col].dropna()
            for pred_str in pred_rxns:
                if isinstance(pred_str, str):
                    predicted_reactions.update(pred_str.split(";"))

    return len(predicted_reactions) / len(unique_reactions) if unique_reactions else 0.0


def evaluate_ensemble_predictions(
    predictions_list: List[pd.DataFrame],
    ground_truth: pd.DataFrame,
    ensemble_method: str = "voting",
    weights: Optional[List[float]] = None,
) -> Dict[str, float]:
    """Evaluate ensemble predictions."""
    if ensemble_method == "voting":
        ensemble_pred = voting_ensemble(predictions_list)
    elif ensemble_method == "weighted_voting":
        if weights is None:
            weights = [1.0] * len(predictions_list)
        ensemble_pred = weighted_voting_ensemble(predictions_list, weights)
    else:
        raise ValueError(f"Unknown ensemble method: {ensemble_method}")

    metrics_result = calculate_reaction_metrics(ground_truth, ensemble_pred)
    return metrics_result


def voting_ensemble(predictions_list: List[pd.DataFrame]) -> pd.DataFrame:
    """Simple voting ensemble."""
    # Placeholder - requires implementation specific to data format
    return predictions_list[0] if predictions_list else pd.DataFrame()


def weighted_voting_ensemble(
    predictions_list: List[pd.DataFrame], weights: List[float]
) -> pd.DataFrame:
    """Weighted voting ensemble."""
    # Placeholder - requires implementation specific to data format
    return predictions_list[0] if predictions_list else pd.DataFrame()
