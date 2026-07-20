"""Cross-validation evaluation workflows for RXNRECer."""

from __future__ import annotations

import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable
from rxnrecer.lib.ml import mlcommon as btools
import numpy as np
import pandas as pd
from tqdm import tqdm

from rxnrecer.config import config as cfg
from rxnrecer.lib.evaluation import eva

display_html_results = eva.display_html_results
show_ec_methods_10_eva_fig = eva.show_ec_methods_10_eva_fig


def read_10fold_res_csv_files(file_paths: Iterable[str | Path]) -> list[pd.DataFrame]:
    """Read tab-separated fold result files."""
    return eva.read_10fold_res_csv_files([str(file_path) for file_path in file_paths])


def make_10folds_labels(
    resdf: list[pd.DataFrame],
    columns_dict: dict[str, str],
    rxn_label_dict: dict[str, int],
    fold_num: int | None = None,
) -> pd.DataFrame:
    """Convert reaction-id columns in per-fold results to label-vector columns."""


    if fold_num is None:
        fold_num = len(resdf)
    if fold_num > len(resdf):
        raise ValueError(
            "fold_num cannot be larger than the number of result DataFrames"
        )

    labeled_folds = []
    for index, frame in enumerate(resdf[:fold_num], start=1):
        fold_frame = frame.copy()
        for src_col, target_col in columns_dict.items():
            fold_frame[target_col] = fold_frame[src_col].apply(
                lambda reaction_id: _encode_reaction_label(
                    reaction_id=str(reaction_id),
                    rxn_label_dict=rxn_label_dict,
                )
            )
        fold_frame["run_fold"] = index
        labeled_folds.append(fold_frame)

    return pd.concat(labeled_folds, axis=0).reset_index(drop=True)


def _encode_reaction_label(
    reaction_id: str,
    rxn_label_dict: dict[str, int],
) -> np.ndarray:
    """Compatibility wrapper around the canonical reaction-label encoder.

    Current evaluation semantics keep `"-"` as a real label when it exists in
    `dict_rxn2id`. Only process sentinels are forced to the all-zero vector.
    """
    return btools.make_label(reaction_id, rxn_label_dict)


def calculate_metrics_parallel(
    res_df: list[pd.DataFrame],
    ground_truth_col: str,
    pred_col: str,
    avg_method: str = "weighted",
    max_workers: int | None = None,
) -> pd.DataFrame:
    """Calculate one metric row per fold in parallel."""
    return eva.calculate_metrics_parallel(
        res_df=res_df,
        ground_truth_col=ground_truth_col,
        pred_col=pred_col,
        avg_method=avg_method,
        max_workers=max_workers,
    )


def statistic_no_res(
    res_df: pd.DataFrame,
    name_col_ec: str,
    name_col_rxn: str,
    type: str = "ec",
) -> pd.DataFrame:
    """Count no-prediction and EC-without-reaction rows per fold."""
    return eva.statistic_no_res(
        res_df=res_df,
        name_col_ec=name_col_ec,
        name_col_rxn=name_col_rxn,
        type=type,
    )


def _method_output_dir(method_type: str) -> Path:
    if method_type == "ec":
        return Path(cfg.RESULTS_DIR) / "intermediate" / "ecmethods"
    if method_type == "direct":
        return Path(cfg.RESULTS_DIR) / "intermediate" / "direct"
    if method_type == "structural":
        return Path(cfg.RESULTS_DIR) / "intermediate" / "structural"
    raise ValueError("method_type must be one of: ec, direct, structural")


def _labels_need_refresh(
    labeled: pd.DataFrame,
    source_col: str,
    label_col: str,
) -> bool:
    """Detect stale cached labels for empty-prediction sentinels."""
    if source_col not in labeled.columns or label_col not in labeled.columns:
        return True

    sentinel_mask = labeled[source_col].isin(["-", "NO-PREDICTION", "EC-WITHOUT-REACTION"])
    if not sentinel_mask.any():
        return False

    return any(int(pd.Series(label).sum()) != 0 for label in labeled.loc[sentinel_mask, label_col])


def _try_read_feather(path: Path) -> pd.DataFrame | None:
    """Return a feather DataFrame, or None if the cache is missing/corrupted."""
    if not path.exists():
        return None
    try:
        return pd.read_feather(path)
    except Exception:
        return None


def _log_progress(verbose: bool, message: str) -> None:
    if verbose:
        print(f"[cross_validation] {message}", flush=True)


def _log_step_elapsed(
    verbose: bool,
    baseline_name: str,
    method_type: str,
    step_name: str,
    started_at: float,
) -> None:
    elapsed = time.time() - started_at
    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: {step_name} finished in {elapsed:.2f}s",
    )


def _normalize_avg_types(avg_types: Iterable[str] | None) -> list[str] | None:
    if avg_types is None:
        return None
    return list(avg_types)


def _avg_types_cache_suffix(avg_types: Iterable[str] | None) -> str:
    normalized = _normalize_avg_types(avg_types)
    if normalized is None:
        return ""
    safe_names = [str(avg_type).replace("/", "_") for avg_type in normalized]
    return "_" + "_".join(safe_names)


def _load_dict_rxn2id(verbose: bool = True) -> dict[str, int]:
    """加载反应编码字典。

    高层评估接口不再要求调用方显式传 `dict_rxn2id`。真正需要重建 label cache
    时，再从统一配置文件中加载，避免 notebook/脚本重复写样板代码。
    """
    with open(cfg.FILE_DS_DICT_RXN2ID, "r") as json_file:
        dict_rxn2id = json.load(json_file)
    _log_progress(verbose, f"loaded dict_rxn2id with {len(dict_rxn2id)} reactions")
    return dict_rxn2id


def _effective_cpu_count() -> int:
    """返回当前进程真正可用的 CPU 数，而不是整机总核数。

    优先级：
    1. SLURM_CPUS_PER_TASK：slurm 任务最明确的 CPU 配额；
    2. os.sched_getaffinity(0)：当前进程实际可绑定的 CPU 集合；
    3. os.cpu_count()：最后兜底。
    """
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


def _method_fold_paths(
    baseline_name: str,
    method_type: str,
    num_folds: int,
) -> list[Path]:
    if method_type == "ec":
        fallback_dir = Path(cfg.DIR_RES_BASELINE) / "results" / "ec_methods" / baseline_name
        base_dir = Path(cfg.RESULTS_DIR) / "intermediate" / "ecmethods" / baseline_name
        if base_dir.exists():
            return [base_dir / f"fold{fold}.tsv" for fold in range(1, num_folds + 1)]
        return [fallback_dir / f"fold{fold}.tsv" for fold in range(1, num_folds + 1)]
    if method_type == "structural":
        base_dir = Path(cfg.RESULTS_DIR) / "intermediate" / "structural"
    elif method_type == "direct":
        base_dir = Path(cfg.RESULTS_DIR) / "intermediate" / "direct"
    else:
        raise ValueError("method_type must be one of: ec, direct, structural")

    # Some direct baselines are stored as one TSV family per fold, containing
    # multiple prediction columns (e.g. `unirep_fold*.tsv` has both
    # `rxn_unirep_cosine` and `rxn_unirep_euclidean`). In that case, the fold TSV
    # is named by the prefix before the last underscore.
    fallback_prefix = baseline_name.rsplit("_", 1)[0] if "_" in baseline_name else None
    paths: list[Path] = []
    for fold in range(1, num_folds + 1):
        primary = base_dir / f"{baseline_name}_fold{fold}.tsv"
        if primary.exists() or not fallback_prefix:
            paths.append(primary)
            continue
        fallback = base_dir / f"{fallback_prefix}_fold{fold}.tsv"
        paths.append(fallback)
    return paths


def _build_label_cache(
    baseline_name: str,
    dict_rxn2id: dict[str, int],
    method_type: str,
    num_folds: int,
    label_file: Path,
    verbose: bool,
) -> Path:
    """从原始 fold TSV 直接重建 label cache。"""
    step_started_at = time.time()
    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: reading {num_folds} fold TSV files and rebuilding labels",
    )
    res = read_10fold_res_csv_files(_method_fold_paths(baseline_name, method_type, num_folds))
    columns_dict = {
        "rxn_groundtruth": "lb_rxn_groundtruth",
        f"rxn_{baseline_name}": f"lb_rxn_{baseline_name}",
    }
    labeled = make_10folds_labels(
        res,
        columns_dict,
        dict_rxn2id,
        fold_num=num_folds,
    )
    labeled.to_feather(label_file)
    _log_progress(verbose, f"{baseline_name}/{method_type}: saved label cache -> {label_file}")
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "label rebuild",
        step_started_at,
    )
    return label_file


def _ensure_label_cache(
    baseline_name: str,
    dict_rxn2id: dict[str, int],
    method_type: str,
    num_folds: int,
    use_label_cache: bool,
    verbose: bool,
) -> Path:
    """确保 label cache 可用；必要时重建，但不触发 overall 评估。"""
    dir_path = _method_output_dir(method_type)
    label_file = dir_path / f"{baseline_name}_10folds_labels_res.feather"

    if not label_file.exists():
        return _build_label_cache(
            baseline_name=baseline_name,
            dict_rxn2id=dict_rxn2id,
            method_type=method_type,
            num_folds=num_folds,
            label_file=label_file,
            verbose=verbose,
        )

    labeled = _try_read_feather(label_file)
    if labeled is None:
        return _build_label_cache(
            baseline_name=baseline_name,
            dict_rxn2id=dict_rxn2id,
            method_type=method_type,
            num_folds=num_folds,
            label_file=label_file,
            verbose=verbose,
        )

    if use_label_cache:
        return label_file

    if _labels_need_refresh(
        labeled=labeled,
        source_col=f"rxn_{baseline_name}",
        label_col=f"lb_rxn_{baseline_name}",
    ):
        return _build_label_cache(
            baseline_name=baseline_name,
            dict_rxn2id=dict_rxn2id,
            method_type=method_type,
            num_folds=num_folds,
            label_file=label_file,
            verbose=verbose,
        )

    return _build_label_cache(
        baseline_name=baseline_name,
        dict_rxn2id=dict_rxn2id,
        method_type=method_type,
        num_folds=num_folds,
        label_file=label_file,
        verbose=verbose,
    )


def _build_overall_metrics_cache(
    baseline_name: str,
    method_type: str,
    label_file: Path,
    metrics_file: Path,
    num_folds: int,
    max_workers: int | None,
    metric_n_jobs: int | None,
    avg_types: Iterable[str] | None,
    verbose: bool,
) -> Path:
    """基于 label cache 计算一次 overall 10-fold 指标并写盘。"""
    step_started_at = time.time()
    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: calculating overall 10-fold metrics "
        f"(fold_workers={max_workers or 'auto'}, avg_types={avg_types or 'all'})",
    )
    labeled = _try_read_feather(label_file)
    if labeled is None:
        raise ValueError(f"Label cache is missing or corrupted: {label_file}")
    metrics = eva.eva_cross_validation(
        res_df=labeled,
        lb_groundtruth="lb_rxn_groundtruth",
        lb_predict=f"lb_rxn_{baseline_name}",
        num_folds=num_folds,
        max_workers=max_workers,
        metric_n_jobs=metric_n_jobs,
        avg_types=avg_types,
    )
    metrics.insert(0, "baselineName", baseline_name)
    metrics.to_feather(metrics_file)
    _log_progress(verbose, f"{baseline_name}/{method_type}: saved metrics cache -> {metrics_file}")
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "overall 10-fold metrics",
        step_started_at,
    )
    return metrics_file


def _ensure_overall_metrics_cache(
    baseline_name: str,
    method_type: str,
    label_file: Path,
    metrics_file: Path,
    num_folds: int,
    use_metrics_cache: bool,
    max_workers: int | None,
    metric_n_jobs: int | None,
    avg_types: Iterable[str] | None,
    verbose: bool,
) -> Path:
    """确保 overall metrics cache 可用。"""
    if use_metrics_cache and metrics_file.exists() and _try_read_feather(metrics_file) is not None:
        _log_progress(verbose, f"{baseline_name}/{method_type}: loading metrics cache -> {metrics_file}")
        return metrics_file

    return _build_overall_metrics_cache(
        baseline_name=baseline_name,
        method_type=method_type,
        label_file=label_file,
        metrics_file=metrics_file,
        num_folds=num_folds,
        max_workers=max_workers,
        metric_n_jobs=metric_n_jobs,
        avg_types=avg_types,
        verbose=verbose,
    )


def get_eval_results(
    baseline_name: str,
    method_type: str,
    num_folds: int = 10,
    force: bool = False,
    use_label_cache: bool = True,
    use_metrics_cache: bool = True,
    max_workers: int | None = None,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """Build or load cached cross-validation metrics for one baseline method."""

    dir_path = _method_output_dir(method_type)
    os.makedirs(dir_path, exist_ok=True)
    avg_types = None
    avg_suffix = ""

    label_file = dir_path / f"{baseline_name}_10folds_labels_res.feather"
    metrics_file = dir_path / f"{baseline_name}_10_folds_metrics_res{avg_suffix}.feather"

    if use_metrics_cache and not use_label_cache:
        use_metrics_cache = False

    if force:
        use_label_cache = False
        use_metrics_cache = False

    dict_rxn2id = _load_dict_rxn2id(verbose=verbose)

    label_file = _ensure_label_cache(
        baseline_name=baseline_name,
        dict_rxn2id=dict_rxn2id,
        method_type=method_type,
        num_folds=num_folds,
        use_label_cache=use_label_cache,
        verbose=verbose,
    )
    metrics_file = _ensure_overall_metrics_cache(
        baseline_name=baseline_name,
        method_type=method_type,
        label_file=label_file,
        metrics_file=metrics_file,
        num_folds=num_folds,
        use_metrics_cache=use_metrics_cache,
        max_workers=max_workers,
        metric_n_jobs=None,
        avg_types=None,
        verbose=verbose,
    )

    labeled = _try_read_feather(label_file)
    metrics = _try_read_feather(metrics_file)
    if labeled is None:
        raise ValueError(f"Label cache is missing or corrupted: {label_file}")
    if metrics is None:
        raise ValueError(f"Metrics cache is missing or corrupted: {metrics_file}")
    std_mean = eva.get_fold_mean_std_metrics(metrics)

    method_rxn_col = f"rxn_{baseline_name}"
    method_ec_col = f"ec_{baseline_name}"
    no_pred = None
    if method_rxn_col in labeled.columns:
        no_pred_type = (
            "ec" if method_type == "ec" and method_ec_col in labeled.columns else "rxn"
        )
        no_pred = statistic_no_res(
            labeled,
            name_col_ec=method_ec_col,
            name_col_rxn=method_rxn_col,
            type=no_pred_type,
        )

    return std_mean, metrics, no_pred


def _identity_cache_paths(
    dir_path: Path,
    baseline_name: str,
    split_enzyme: bool,
    avg_suffix: str,
) -> tuple[Path, Path, Path]:
    """构造 identity-bin 评价结果的 3 个缓存文件路径。

    这个函数只负责路径命名，不读写文件。缓存分成三类：
    1. std_mean：每个 bin 的 10-fold 均值和标准差，通常用于画图和表格展示。
    2. metrics：每个 bin、每个 fold 的原始指标明细，便于后续排查。
    3. no_pred：每个 bin 中无预测样本的统计信息。

    当前目录只保留一份“最新”缓存，不再把版本号写进文件名。
    """
    suffix = "_enzyme" if split_enzyme else ""
    cache_std_mean_file = (
        dir_path / f"{baseline_name}_identity_bins{suffix}{avg_suffix}_std_mean.feather"
    )
    cache_metrics_file = (
        dir_path / f"{baseline_name}_identity_bins{suffix}{avg_suffix}_metrics.feather"
    )
    cache_no_pred_file = (
        dir_path / f"{baseline_name}_identity_bins{suffix}{avg_suffix}_no_pred.feather"
    )
    return cache_std_mean_file, cache_metrics_file, cache_no_pred_file


def _combined_identity_cache_paths(
    dir_path: Path,
    baseline_name: str,
    avg_suffix: str,
) -> tuple[Path, Path, Path]:
    """构造“完整 identity 大表”的缓存路径。

    这套缓存对应论文最终使用的大表，内部已经同时包含：
    - Overall(all)
    - identity bins(all)
    - Overall_enzyme / Overall_non_enzyme
    - identity bins(enzyme / non_enzyme)
    """
    cache_std_mean_file = (
        dir_path / f"{baseline_name}_identity_bins_combined{avg_suffix}_std_mean.feather"
    )
    cache_metrics_file = (
        dir_path / f"{baseline_name}_identity_bins_combined{avg_suffix}_metrics.feather"
    )
    cache_no_pred_file = (
        dir_path / f"{baseline_name}_identity_bins_combined{avg_suffix}_no_pred.feather"
    )
    return cache_std_mean_file, cache_metrics_file, cache_no_pred_file


def _load_identity_cache(
    baseline_name: str,
    method_type: str,
    cache_std_mean_file: Path,
    cache_metrics_file: Path,
    cache_no_pred_file: Path,
    verbose: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame] | None:
    """尝试读取 identity-bin 缓存。

    只有 std_mean 和 metrics 两个核心缓存都存在且能正常读取时，才认为缓存可用。
    no_pred 是辅助统计，缺失时返回空 DataFrame，不阻断主流程。

    返回 None 表示缓存不可用，调用方应该重新计算。
    """
    if not cache_std_mean_file.exists() or not cache_metrics_file.exists():
        return None

    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: loading identity-bin cache -> {cache_std_mean_file}",
    )
    final_std_mean = _try_read_feather(cache_std_mean_file)
    final_metrics = _try_read_feather(cache_metrics_file)
    final_no_pred = _try_read_feather(cache_no_pred_file)
    if final_std_mean is None or final_metrics is None:
        return None
    if final_no_pred is None:
        final_no_pred = pd.DataFrame()
    return final_std_mean, final_metrics, final_no_pred
def _ensure_identity_label_cache(
    baseline_name: str,
    dict_rxn2id: dict[str, int],
    method_type: str,
    num_folds: int,
    force: bool,
    use_label_cache: bool,
    use_metrics_cache: bool,
    max_workers: int | None,
    metric_n_jobs: int | None,
    avg_types: Iterable[str] | None,
    verbose: bool,
) -> Path:
    """identity 评价只依赖 label cache，本函数不再顺带触发 overall metrics。"""
    del force, use_metrics_cache, max_workers, metric_n_jobs, avg_types
    return _ensure_label_cache(
        baseline_name=baseline_name,
        dict_rxn2id=dict_rxn2id,
        method_type=method_type,
        num_folds=num_folds,
        use_label_cache=use_label_cache,
        verbose=verbose,
    )


def _load_validation_identity_features(num_folds: int, verbose: bool, baseline_name: str, method_type: str) -> pd.DataFrame:
    """读取每个 validation fold 的 sequence identity 特征。

    评价分箱依赖 validation/fold*/valid.feather 中的 `pident` 字段：
    - pident == 0.0：该验证样本相对训练集没有 BLAST hit，归入 `No BLAST Hit`；
    - pident > 0.0：按 identity_bins 切到 `(0,30]`、`(30,40]` 等区间。

    如果 valid.feather 里存在 `isenzyme`，也会保留下来。EC 方法在不拆 enzyme/non-enzyme
    时会只评价 enzyme 样本，避免大量非酶 true negative 抬高准确率。
    """
    _log_progress(verbose, f"{baseline_name}/{method_type}: loading validation pident features")
    valid_list = []
    for fold in range(1, num_folds + 1):
        valid_path = Path(cfg.DIR_DATASET) / "validation" / f"fold{fold}" / "valid.feather"
        if valid_path.exists():
            temp_df = pd.read_feather(valid_path)
            cols = ["uniprot_id", "pident"]
            if "isenzyme" in temp_df.columns:
                cols.append("isenzyme")
            df_valid = temp_df[cols].copy()
            df_valid["run_fold"] = fold
            valid_list.append(df_valid)
        else:
            print(f"Warning: Validation file not found: {valid_path}")

    if not valid_list:
        raise FileNotFoundError("No validation feather files found.")

    return pd.concat(valid_list, axis=0).reset_index(drop=True)


def _merge_labels_with_identity_bins(
    labeled: pd.DataFrame,
    valid_features: pd.DataFrame,
    identity_bins: list[float],
    baseline_name: str,
    method_type: str,
    verbose: bool,
) -> pd.DataFrame:
    """把预测 label 表和 validation identity 特征合并，生成基础 identity bin。

    分箱逻辑是整个函数最关键的地方：
    - `No BLAST Hit` 只由 validation 特征里的 `pident == 0.0` 决定；
    - 它不是某个方法是否有预测的标记，所以 RXNRECer、embedding 方法、EC 方法等
      在这个 bin 内仍然应该正常计算指标；
    - `pident > 0.0` 的样本通过 pandas.cut 进入 identity 区间。

    注意这里不使用 `include_lowest=True`。因为 pident==0 已经单独归入
    `No BLAST Hit`，最低普通区间应显示为 `(0.0, 30.0]`，不需要 pandas 为了
    包含左边界而生成 `(-0.001, 30.0]` 这种显示标签。
    """
    _log_progress(verbose, f"{baseline_name}/{method_type}: merging labels with validation features")
    merged = labeled.merge(valid_features, on=["uniprot_id", "run_fold"], how="inner")

    no_blast_hit_mask = merged["pident"] == 0.0
    merged["identity_bin_base"] = ""
    merged.loc[no_blast_hit_mask, "identity_bin_base"] = "No BLAST Hit"
    merged.loc[~no_blast_hit_mask, "identity_bin_base"] = pd.cut(
        merged.loc[~no_blast_hit_mask, "pident"],
        bins=identity_bins,
        include_lowest=False,
    ).astype(str)

    return merged


def _prepare_identity_eval_frame(
    merged_base: pd.DataFrame,
    baseline_name: str,
    method_type: str,
    split_enzyme: bool,
    verbose: bool,
) -> pd.DataFrame:
    """基于同一份 merge 结果，派生 split / non-split 两种 identity 视图。"""
    merged = merged_base.copy()
    merged["identity_bin"] = merged["identity_bin_base"]

    if split_enzyme and "isenzyme" in merged.columns:
        merged["identity_bin"] = merged.apply(
            lambda row: f"{row['identity_bin_base']}_enzyme"
            if row["isenzyme"]
            else f"{row['identity_bin_base']}_non_enzyme",
            axis=1,
        )

    if method_type == "ec" and "isenzyme" in merged.columns and not split_enzyme:
        merged = merged[merged["isenzyme"] == True].copy()
        _log_progress(
            verbose,
            f"{baseline_name}/{method_type}: restricted EC identity-bin evaluation to enzyme rows ({len(merged):,})",
        )

    return merged


def _prepare_enzyme_binary_eval_frame(
    merged_base: pd.DataFrame,
    baseline_name: str,
) -> pd.DataFrame:
    """基于同一份 merge 结果，派生 enzyme-vs-non-enzyme 二分类视图。

    这条任务与主反应多标签任务并行存在：
    - 主任务：预测具体反应（含 `-` 类）；
    - 二分类任务：判断该样本是否为 enzyme。

    二分类标签定义：
    - truth：`rxn_groundtruth != "-"` 记为 1，否则 0；
    - pred：
      - 原始预测字符串中包含 `RHEA:` 记为 1；
      - 原始预测字符串恰好为 `"-"` 记为 0；
      - `NO-PREDICTION` / `EC-WITHOUT-REACTION` 视为“无效预测”，无论真值是什么都记错。

    最后一条不能直接用普通 0/1 映射表达，因为：
    - 若真值为 non-enzyme，把无效预测记成 0 会被误算为正确；
    - 若真值为 enzyme，把无效预测记成 0 才是错误。
    因此这里把无效预测统一映射成 `1 - truth`，确保它始终记错。
    """
    merged = merged_base.copy()
    merged["identity_bin"] = merged["identity_bin_base"]
    merged["truth_enzyme_binary"] = merged["rxn_groundtruth"].astype(str).ne("-").astype(int)

    pred_col = f"rxn_{baseline_name}"
    if pred_col not in merged.columns:
        raise KeyError(f"Missing prediction column for binary enzyme task: {pred_col}")

    pred_str = merged[pred_col].astype(str)
    rhea_mask = pred_str.str.contains("RHEA:", na=False)
    dash_mask = pred_str.eq("-")
    invalid_mask = pred_str.isin(["NO-PREDICTION", "EC-WITHOUT-REACTION"])

    merged["pred_enzyme_binary"] = 0
    merged.loc[rhea_mask, "pred_enzyme_binary"] = 1
    merged.loc[dash_mask, "pred_enzyme_binary"] = 0
    merged.loc[invalid_mask, "pred_enzyme_binary"] = 1 - merged.loc[invalid_mask, "truth_enzyme_binary"]
    return merged


def _build_identity_groups_to_eval(
    merged: pd.DataFrame,
    split_enzyme: bool,
    verbose: bool,
    baseline_name: str,
    method_type: str,
) -> tuple[list[tuple[str, pd.DataFrame]], int, pd.Series, pd.Series]:
    """根据 `identity_bin` 拆出需要评价的分组。

    返回值包括：
    - groups_to_eval：每个元素是 `(bin_name, bin_dataframe)`；
    - total_samples：当前评价范围内的总样本数，用来计算 bin_percentage；
    - bin_sizes：每个 bin 的样本数。

    当 `split_enzyme=True` 时，除普通 identity bin 外，还会额外加入
    `Overall_enzyme` 和 `Overall_non_enzyme` 两个整体分组。
    """
    total_samples = len(merged)
    total_samples_by_fold = merged.groupby("run_fold").size()
    bin_sizes = merged["identity_bin"].value_counts(dropna=False)
    groups_to_eval = list(merged.groupby("identity_bin"))

    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: prepared {len(groups_to_eval)} identity bins from {total_samples:,} rows",
    )

    if split_enzyme and "isenzyme" in merged.columns:
        enzyme_overall_group = merged[merged["isenzyme"] == True].copy()
        enzyme_overall_group["identity_bin"] = "Overall_enzyme"

        non_enzyme_overall_group = merged[merged["isenzyme"] == False].copy()
        non_enzyme_overall_group["identity_bin"] = "Overall_non_enzyme"

        groups_to_eval.append(("Overall_enzyme", enzyme_overall_group))
        groups_to_eval.append(("Overall_non_enzyme", non_enzyme_overall_group))

    return groups_to_eval, total_samples, bin_sizes, total_samples_by_fold


def _base_bin_name(identity_bin: str) -> str:
    """把带后缀的 identity_bin 还原成基础 similarity bin 名称。"""
    bin_str = str(identity_bin)
    for suffix in (
        "_enzyme_vs_non_enzyme",
        "_non_enzyme",
        "_enzyme",
    ):
        if bin_str.endswith(suffix):
            return bin_str[: -len(suffix)]
    return bin_str


def _compute_shared_bin_percentage_stats(merged_base: pd.DataFrame) -> dict[str, dict[str, object]]:
    """基于全体 validation 样本，计算 fold-wise similarity bin 占比。

    这套统计只和：
    - fold
    - identity_bin_base
    有关，不再和方法或 enzyme/non-enzyme 子集绑定。

    返回值中每个 bin 会包含：
    - fold_pct_by_fold：每个 fold 的 bin 占比
    - mean / std：10-fold 均值和标准差
    - overall：10 折合并后的总体占比
    """
    total_samples = len(merged_base)
    total_samples_by_fold = merged_base.groupby("run_fold").size()

    stats: dict[str, dict[str, object]] = {}

    overall_fold_pct = pd.Series(1.0, index=total_samples_by_fold.index, dtype=float)
    stats["Overall"] = {
        "fold_pct_by_fold": overall_fold_pct,
        "mean": 1.0,
        "std": 0.0,
        "overall": 1.0,
    }

    for bin_name, group in merged_base.groupby("identity_bin_base"):
        group_counts_by_fold = group.groupby("run_fold").size()
        fold_pct_by_fold = total_samples_by_fold.index.to_series().map(
            lambda fold: group_counts_by_fold.get(fold, 0) / total_samples_by_fold.get(fold, 1)
        ).astype(float)
        stats[str(bin_name)] = {
            "fold_pct_by_fold": fold_pct_by_fold,
            "mean": float(fold_pct_by_fold.mean()) if not fold_pct_by_fold.empty else 0.0,
            "std": float(fold_pct_by_fold.std()) if len(fold_pct_by_fold) > 1 else 0.0,
            "overall": (len(group) / total_samples) if total_samples > 0 else 0.0,
        }

    return stats


def _finalize_identity_group_result(
    bin_val: str,
    group: pd.DataFrame,
    metrics: pd.DataFrame,
    baseline_name: str,
    method_type: str,
    total_samples: int,
    bin_sizes: pd.Series,
    total_samples_by_fold: pd.Series,
    shared_bin_pct_stats: dict[str, dict[str, object]],
    skip_blast_zero_rule: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """给一个 identity bin 的 fold 指标补齐元信息并计算均值/标准差。

    `metrics` 是该 bin 内每个 fold 的原始评价结果。这里会补三类信息：
    - baselineName：方法名称；
    - identity_bin：当前 bin 名称；
    - bin_percentage：当前 bin 样本数占总样本数的比例。

    特殊规则：只有 BLAST 方法本身在 `No BLAST Hit` bin 中才强制置 0。
    其他方法并不依赖 BLAST hit，必须在这个 bin 里正常评价。
    """
    _ = (total_samples, bin_sizes, total_samples_by_fold, skip_blast_zero_rule)
    shared_pct = shared_bin_pct_stats[_base_bin_name(str(bin_val))]
    fold_pct_by_fold = shared_pct["fold_pct_by_fold"]
    fold_bin_pct = metrics["runFold"].map(lambda fold: fold_pct_by_fold.get(fold, 0.0)).astype(float)
    bin_pct_mean = float(shared_pct["mean"])
    bin_pct_std = float(shared_pct["std"])
    bin_pct_display = f"{bin_pct_mean:.4f} ± {bin_pct_std:.4f}"
    bin_pct_overall = float(shared_pct["overall"])

    if baseline_name == "blast" and str(bin_val).startswith("No BLAST Hit"):
        metrics = metrics.copy()
        metrics[["mAccuracy", "mPrecision", "mRecall", "mF1"]] = 0.0

    metrics.insert(0, "baselineName", baseline_name)
    metrics.insert(1, "identity_bin", str(bin_val))
    metrics.insert(2, "bin_percentage", fold_bin_pct.values)

    std_mean = eva.get_fold_mean_std_metrics(metrics)
    std_mean.insert(0, "identity_bin", str(bin_val))
    std_mean.insert(1, "bin_percentage", bin_pct_display)
    std_mean.insert(2, "bin_percentage_mean", bin_pct_mean)
    std_mean.insert(3, "bin_percentage_std", bin_pct_std)
    std_mean.insert(4, "bin_percentage_overall", bin_pct_overall)

    method_rxn_col = f"rxn_{baseline_name}"
    method_ec_col = f"ec_{baseline_name}"
    if method_rxn_col not in group.columns:
        return metrics, std_mean, None

    no_pred_type = "ec" if method_type == "ec" and method_ec_col in group.columns else "rxn"
    no_pred = statistic_no_res(
        group,
        name_col_ec=method_ec_col,
        name_col_rxn=method_rxn_col,
        type=no_pred_type,
    )
    no_pred.insert(0, "identity_bin", str(bin_val))
    no_pred.insert(
        1,
        "bin_percentage",
        no_pred["run_fold"].map(
            lambda fold: fold_pct_by_fold.get(fold, 0.0)
        ).astype(float),
    )

    return metrics, std_mean, no_pred


def _evaluate_identity_group_serial(
    bin_val: str,
    group: pd.DataFrame,
    baseline_name: str,
    method_type: str,
    num_folds: int,
    max_workers: int | None,
    metric_n_jobs: int | None,
    avg_types: Iterable[str] | None,
    total_samples: int,
    bin_sizes: pd.Series,
    total_samples_by_fold: pd.Series,
    shared_bin_pct_stats: dict[str, dict[str, object]],
    verbose: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None] | None:
    """评价单个 identity bin。

    这个函数用于 `bin_workers <= 1` 的路径。它仍然可以通过 `max_workers`
    在 bin 内部并行计算 fold，但外层 bin 是串行的。优点是逻辑简单；
    缺点是在 bin 很多、每个 bin 样本量较大时，CPU 利用率不如 flat 并行路径。
    """
    if group.empty:
        return None

    try:
        _log_progress(
            verbose,
            f"{baseline_name}/{method_type}: calculating bin {bin_val} ({len(group):,} rows)",
        )
        metrics = eva.eva_cross_validation(
            res_df=group,
            lb_groundtruth="lb_rxn_groundtruth",
            lb_predict=f"lb_rxn_{baseline_name}",
            num_folds=num_folds,
            max_workers=max_workers,
            metric_n_jobs=metric_n_jobs,
            avg_types=avg_types,
        )
        return _finalize_identity_group_result(
            bin_val=bin_val,
            group=group,
            metrics=metrics,
            baseline_name=baseline_name,
            method_type=method_type,
            total_samples=total_samples,
            bin_sizes=bin_sizes,
            total_samples_by_fold=total_samples_by_fold,
            shared_bin_pct_stats=shared_bin_pct_stats,
        )
    except Exception as e:
        print(f"Warning: Failed to calculate metrics for bin {bin_val}: {e}")
        return None


def _evaluate_identity_groups_flat_parallel(
    groups_to_eval: list[tuple[str, pd.DataFrame]],
    baseline_name: str,
    method_type: str,
    num_folds: int,
    max_workers: int | None,
    metric_n_jobs: int | None,
    bin_workers: int,
    avg_types: Iterable[str] | None,
    total_samples: int,
    bin_sizes: pd.Series,
    total_samples_by_fold: pd.Series,
    shared_bin_pct_stats: dict[str, dict[str, object]],
    verbose: bool,
) -> list[tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None] | None]:
    """把 `(identity bin, fold)` 展平成一个全局线程池并行评价。

    这里不用“bin 外层线程池 + fold 内层线程池”的嵌套并行，因为嵌套线程池容易：
    - 线程数失控，导致 CPU 过度切换；
    - 外层任务等待内层任务，实际吞吐下降；
    - 在某些环境中出现看起来像卡住的同步等待。

    flat 并行的做法是把每个 bin 的每个 fold 都作为独立任务提交到同一个线程池。
    所有任务完成后，再按 bin 拼回 10-fold metrics，并进入统一的 finalize 流程。
    """
    non_empty_groups = [
        (index, bin_val, group)
        for index, (bin_val, group) in enumerate(groups_to_eval)
        if not group.empty
    ]
    if not non_empty_groups:
        return []

    effective_cpus = _effective_cpu_count()
    flat_workers = min(
        len(non_empty_groups) * num_folds,
        max(1, bin_workers),
        effective_cpus,
    )
    metric_jobs = metric_n_jobs if metric_n_jobs is not None else 1
    total_tasks = len(non_empty_groups) * num_folds

    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: calculating identity-bin metrics as {total_tasks} "
        f"(bin, fold) tasks with {flat_workers} workers, avg_types={avg_types or 'all'}",
    )

    group_fold_metrics = {index: [None] * num_folds for index, _, _ in non_empty_groups}
    failed_groups = set()
    with ThreadPoolExecutor(max_workers=flat_workers) as executor:
        futures = {}
        for index, bin_val, group in non_empty_groups:
            for fold_idx, runfold in enumerate(range(1, num_folds + 1)):
                fold_df = group[group.run_fold == runfold].reset_index(drop=True)
                futures[
                    executor.submit(
                        eva.eva_one_fold,
                        fold_df,
                        "lb_rxn_groundtruth",
                        f"lb_rxn_{baseline_name}",
                        runfold,
                        metric_jobs,
                        avg_types,
                    )
                ] = (index, bin_val, fold_idx)

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc=f"{baseline_name}/{method_type} bins",
            disable=not verbose,
        ):
            index, bin_val, fold_idx = futures[future]
            try:
                group_fold_metrics[index][fold_idx] = future.result()
            except Exception as e:
                failed_groups.add(index)
                print(f"Warning: Failed to calculate metrics for bin {bin_val}: {e}")

    group_results = [None] * len(groups_to_eval)
    group_lookup = {index: (bin_val, group) for index, bin_val, group in non_empty_groups}
    for index, fold_metrics in group_fold_metrics.items():
        if index in failed_groups or any(metric is None for metric in fold_metrics):
            continue
        bin_val, group = group_lookup[index]
        metrics = pd.concat(fold_metrics, axis=0).reset_index(drop=True)
        group_results[index] = _finalize_identity_group_result(
            bin_val=bin_val,
            group=group,
            metrics=metrics,
            baseline_name=baseline_name,
            method_type=method_type,
            total_samples=total_samples,
            bin_sizes=bin_sizes,
            total_samples_by_fold=total_samples_by_fold,
            shared_bin_pct_stats=shared_bin_pct_stats,
        )

    return group_results


def _evaluate_identity_groups(
    groups_to_eval: list[tuple[str, pd.DataFrame]],
    baseline_name: str,
    method_type: str,
    num_folds: int,
    max_workers: int | None,
    metric_n_jobs: int | None,
    bin_workers: int | None,
    avg_types: Iterable[str] | None,
    total_samples: int,
    bin_sizes: pd.Series,
    total_samples_by_fold: pd.Series,
    shared_bin_pct_stats: dict[str, dict[str, object]],
    verbose: bool,
) -> list[tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None] | None]:
    """根据 bin_workers 选择串行分 bin 或 flat 并行分 bin 的评价路径。"""
    if bin_workers and bin_workers > 1:
        return _evaluate_identity_groups_flat_parallel(
            groups_to_eval=groups_to_eval,
            baseline_name=baseline_name,
            method_type=method_type,
            num_folds=num_folds,
            max_workers=max_workers,
            metric_n_jobs=metric_n_jobs,
            bin_workers=bin_workers,
            avg_types=avg_types,
            total_samples=total_samples,
            bin_sizes=bin_sizes,
            total_samples_by_fold=total_samples_by_fold,
            shared_bin_pct_stats=shared_bin_pct_stats,
            verbose=verbose,
        )

    return [
        _evaluate_identity_group_serial(
            bin_val=bin_val,
            group=group,
            baseline_name=baseline_name,
            method_type=method_type,
            num_folds=num_folds,
            max_workers=max_workers,
            metric_n_jobs=metric_n_jobs,
            avg_types=avg_types,
            total_samples=total_samples,
            bin_sizes=bin_sizes,
            total_samples_by_fold=total_samples_by_fold,
            shared_bin_pct_stats=shared_bin_pct_stats,
            verbose=verbose,
        )
        for bin_val, group in groups_to_eval
    ]


def _combine_identity_group_results(
    group_results: list[tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None] | None],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """合并所有 identity bin 的结果。

    输入是每个 bin 的 `(metrics, std_mean, no_pred)`。这里会分别纵向拼接：
    - 所有 bin 的 fold 明细指标；
    - 所有 bin 的均值/标准差；
    - 所有 bin 的 no-prediction 统计。

    如果所有 bin 都计算失败或为空，返回三个空 DataFrame，让调用方保持兼容。
    """
    all_metrics = []
    all_std_mean = []
    all_no_pred = []

    for group_result in group_results:
        if group_result is None:
            continue
        metrics, std_mean, no_pred = group_result
        all_metrics.append(metrics)
        all_std_mean.append(std_mean)
        if no_pred is not None:
            all_no_pred.append(no_pred)

    if not all_std_mean:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    final_std_mean = pd.concat(all_std_mean, axis=0).reset_index(drop=True)
    final_metrics = pd.concat(all_metrics, axis=0).reset_index(drop=True)
    final_no_pred = (
        pd.concat(all_no_pred, axis=0).reset_index(drop=True)
        if all_no_pred
        else pd.DataFrame()
    )
    return final_std_mean, final_metrics, final_no_pred


def _save_identity_cache(
    final_std_mean: pd.DataFrame,
    final_metrics: pd.DataFrame,
    final_no_pred: pd.DataFrame,
    cache_std_mean_file: Path,
    cache_metrics_file: Path,
    cache_no_pred_file: Path,
    baseline_name: str,
    method_type: str,
    verbose: bool,
) -> None:
    """保存 identity-bin 评价缓存。

    std_mean 和 metrics 是主结果，始终保存。no_pred 是辅助统计，只有非空时保存，
    避免生成没有信息量的空 feather 文件。
    """
    _log_progress(verbose, f"{baseline_name}/{method_type}: saving identity-bin caches")
    final_std_mean.to_feather(cache_std_mean_file)
    final_metrics.to_feather(cache_metrics_file)
    if not final_no_pred.empty:
        final_no_pred.to_feather(cache_no_pred_file)
    _log_progress(verbose, f"{baseline_name}/{method_type}: identity-bin evaluation finished")


def _evaluate_identity_from_merged(
    merged_base: pd.DataFrame,
    baseline_name: str,
    method_type: str,
    split_enzyme: bool,
    num_folds: int,
    max_workers: int | None,
    metric_n_jobs: int | None,
    bin_workers: int | None,
    avg_types: Iterable[str] | None,
    cache_std_mean_file: Path,
    cache_metrics_file: Path,
    cache_no_pred_file: Path,
    force: bool,
    verbose: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """对同一份 merge 后数据执行一次 identity 评估。"""
    if not force:
        cached = _load_identity_cache(
            baseline_name=baseline_name,
            method_type=method_type,
            cache_std_mean_file=cache_std_mean_file,
            cache_metrics_file=cache_metrics_file,
            cache_no_pred_file=cache_no_pred_file,
            verbose=verbose,
        )
        if cached is not None:
            return cached

    step_started_at = time.time()
    merged = _prepare_identity_eval_frame(
        merged_base=merged_base,
        baseline_name=baseline_name,
        method_type=method_type,
        split_enzyme=split_enzyme,
        verbose=verbose,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "identity view preparation",
        step_started_at,
    )

    step_started_at = time.time()
    shared_bin_pct_stats = _compute_shared_bin_percentage_stats(merged_base)
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "shared bin percentage preparation",
        step_started_at,
    )

    step_started_at = time.time()
    groups_to_eval, total_samples, bin_sizes, total_samples_by_fold = _build_identity_groups_to_eval(
        merged=merged,
        split_enzyme=split_enzyme,
        verbose=verbose,
        baseline_name=baseline_name,
        method_type=method_type,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "identity group preparation",
        step_started_at,
    )

    step_started_at = time.time()
    group_results = _evaluate_identity_groups(
        groups_to_eval=groups_to_eval,
        baseline_name=baseline_name,
        method_type=method_type,
        num_folds=num_folds,
        max_workers=max_workers,
        metric_n_jobs=metric_n_jobs,
        bin_workers=bin_workers,
        avg_types=avg_types,
        total_samples=total_samples,
        bin_sizes=bin_sizes,
        total_samples_by_fold=total_samples_by_fold,
        shared_bin_pct_stats=shared_bin_pct_stats,
        verbose=verbose,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "identity-bin metrics",
        step_started_at,
    )

    step_started_at = time.time()
    final_std_mean, final_metrics, final_no_pred = _combine_identity_group_results(group_results)
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "identity result combine",
        step_started_at,
    )

    if final_std_mean.empty:
        return final_std_mean, final_metrics, final_no_pred

    step_started_at = time.time()
    _save_identity_cache(
        final_std_mean=final_std_mean,
        final_metrics=final_metrics,
        final_no_pred=final_no_pred,
        cache_std_mean_file=cache_std_mean_file,
        cache_metrics_file=cache_metrics_file,
        cache_no_pred_file=cache_no_pred_file,
        baseline_name=baseline_name,
        method_type=method_type,
        verbose=verbose,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "identity cache save",
        step_started_at,
    )

    return final_std_mean, final_metrics, final_no_pred


def _evaluate_enzyme_binary_from_merged(
    merged_base: pd.DataFrame,
    baseline_name: str,
    method_type: str,
    num_folds: int,
    max_workers: int | None,
    cache_std_mean_file: Path,
    cache_metrics_file: Path,
    force: bool,
    verbose: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """基于 merged_base 计算 enzyme-vs-non-enzyme 二分类指标。"""
    if not force:
        cached = _load_identity_cache(
            baseline_name=baseline_name,
            method_type=method_type,
            cache_std_mean_file=cache_std_mean_file,
            cache_metrics_file=cache_metrics_file,
            cache_no_pred_file=cache_metrics_file.parent / "__unused__.feather",
            verbose=verbose,
        )
        if cached is not None:
            cached_std, cached_metrics, _ = cached
            return cached_std, cached_metrics

    step_started_at = time.time()
    merged = _prepare_enzyme_binary_eval_frame(
        merged_base=merged_base,
        baseline_name=baseline_name,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "enzyme-vs-non-enzyme view preparation",
        step_started_at,
    )

    total_samples = len(merged)
    shared_bin_pct_stats = _compute_shared_bin_percentage_stats(merged_base)
    groups_to_eval = [("Overall", merged.copy())] + list(merged.groupby("identity_bin"))

    _log_progress(
        verbose,
        f"{baseline_name}/{method_type}: prepared {len(groups_to_eval)} enzyme-vs-non-enzyme bins "
        f"from {total_samples:,} rows",
    )

    all_metrics = []
    all_std_mean = []
    for bin_val, group in groups_to_eval:
        metrics = eva.eva_cross_validation(
            res_df=group,
            lb_groundtruth="truth_enzyme_binary",
            lb_predict="pred_enzyme_binary",
            num_folds=num_folds,
            max_workers=max_workers,
            metric_n_jobs=None,
            avg_types=None,
        )

        shared_pct = shared_bin_pct_stats[_base_bin_name(str(bin_val))]
        fold_pct_by_fold = shared_pct["fold_pct_by_fold"]
        fold_bin_pct = metrics["runFold"].map(lambda fold: fold_pct_by_fold.get(fold, 0.0)).astype(float)
        bin_pct_mean = float(shared_pct["mean"])
        bin_pct_std = float(shared_pct["std"])
        bin_pct_display = f"{bin_pct_mean:.4f} ± {bin_pct_std:.4f}"
        bin_pct_overall = float(shared_pct["overall"])

        metrics = metrics.copy()
        if baseline_name == "blast" and str(bin_val) == "No BLAST Hit":
            metrics[["mAccuracy", "mPrecision", "mRecall", "mF1"]] = 0.0
        metrics.insert(0, "baselineName", baseline_name)
        metrics.insert(1, "identity_bin", f"{bin_val}_enzyme_vs_non_enzyme")
        metrics.insert(2, "bin_percentage", fold_bin_pct.values)
        all_metrics.append(metrics)

        std_mean = eva.get_fold_mean_std_metrics(metrics)
        std_mean.insert(0, "identity_bin", f"{bin_val}_enzyme_vs_non_enzyme")
        std_mean.insert(1, "bin_percentage", bin_pct_display)
        std_mean.insert(2, "bin_percentage_mean", bin_pct_mean)
        std_mean.insert(3, "bin_percentage_std", bin_pct_std)
        std_mean.insert(4, "bin_percentage_overall", bin_pct_overall)
        all_std_mean.append(std_mean)

    final_std_mean = pd.concat(all_std_mean, axis=0).reset_index(drop=True)
    final_metrics = pd.concat(all_metrics, axis=0).reset_index(drop=True)
    final_std_mean.to_feather(cache_std_mean_file)
    final_metrics.to_feather(cache_metrics_file)
    return final_std_mean, final_metrics


def get_eval_results_by_identity(
    baseline_name: str,
    method_type: str,
    num_folds: int = 10,
    identity_bins: list[float] | None = None,
    force: bool = False,
    use_label_cache: bool = True,
    use_metrics_cache: bool = True,
    split_enzyme: bool = False,
    combine_views: bool = True,
    bin_workers: int = 10 ,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """按 sequence identity bin 计算一套评估结果。

    默认返回一张“完整 identity 大表”：
    1. Overall(all)
    2. 每个 identity bin 的 all
    3. Overall_enzyme / Overall_non_enzyme
    4. 每个 identity bin 的 enzyme / non_enzyme

    这样 notebook 只需要调一次，不再分别跑 split / non-split 再手工合并。
    """
    if identity_bins is None:
        identity_bins = [0.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

    # 统一并发策略：
    # - bin_workers: 唯一对外暴露的并发控制，主要控制 identity-bin 评价吞吐；
    # - overall 阶段 fold 并发直接复用 bin_workers；
    # - 单 fold 内部 metric joblib 不再单独开放控制，默认交给 eva 模块自动处理。
    fold_workers = bin_workers
    metric_n_jobs = None

    dir_path = _method_output_dir(method_type)
    os.makedirs(dir_path, exist_ok=True)
    avg_types = None
    avg_suffix = ""

    if combine_views:
        combined_paths = _combined_identity_cache_paths(
            dir_path=dir_path,
            baseline_name=baseline_name,
            avg_suffix=avg_suffix,
        )
        if not force:
            cached = _load_identity_cache(
                baseline_name=baseline_name,
                method_type=method_type,
                cache_std_mean_file=combined_paths[0],
                cache_metrics_file=combined_paths[1],
                cache_no_pred_file=combined_paths[2],
                verbose=verbose,
            )
            if cached is not None:
                cached_std, cached_metrics, cached_no_pred = cached
                has_binary_rows = (
                    "identity_bin" in cached_std.columns
                    and cached_std["identity_bin"].astype(str).str.endswith("_enzyme_vs_non_enzyme").any()
                )
                if has_binary_rows:
                    return cached_std, cached_metrics, cached_no_pred

        # Step 1. overall 只算一次
        overall_std, overall_metrics, overall_no_pred = get_eval_results(
            baseline_name=baseline_name,
            method_type=method_type,
            num_folds=num_folds,
            force=force,
            use_label_cache=use_label_cache,
            use_metrics_cache=use_metrics_cache,
            max_workers=fold_workers,
            verbose=verbose,
        )

        # Step 2. label / validation / merge 只做一次
        #
        # `get_eval_results()` 上一步已经保证 label cache 存在：
        # - 如果 force=True，它刚刚重建完 label cache；
        # - 如果 force=False，它已经确认现有 label cache 可用。
        #
        # 因此这里优先直接复用同一个 label cache，避免 stage=all 中
        # “overall 算一次后，identity 前又 rebuild 一次 labels”的重复开销。
        label_file = dir_path / f"{baseline_name}_10folds_labels_res.feather"
        if _try_read_feather(label_file) is None:
            label_file = _ensure_identity_label_cache(
                baseline_name=baseline_name,
                dict_rxn2id=_load_dict_rxn2id(verbose=verbose),
                method_type=method_type,
                num_folds=num_folds,
                force=force,
                use_label_cache=True,
                use_metrics_cache=use_metrics_cache,
                max_workers=fold_workers,
                metric_n_jobs=metric_n_jobs,
                avg_types=None,
                verbose=verbose,
            )
        _log_progress(verbose, f"{baseline_name}/{method_type}: loading label cache -> {label_file}")
        step_started_at = time.time()
        labeled = pd.read_feather(label_file)
        _log_step_elapsed(
            verbose,
            baseline_name,
            method_type,
            "label cache load",
            step_started_at,
        )

        step_started_at = time.time()
        valid_features = _load_validation_identity_features(
            num_folds=num_folds,
            verbose=verbose,
            baseline_name=baseline_name,
            method_type=method_type,
        )
        _log_step_elapsed(
            verbose,
            baseline_name,
            method_type,
            "validation identity feature load",
            step_started_at,
        )

        step_started_at = time.time()
        merged_base = _merge_labels_with_identity_bins(
            labeled=labeled,
            valid_features=valid_features,
            identity_bins=identity_bins,
            baseline_name=baseline_name,
            method_type=method_type,
            verbose=verbose,
        )
        _log_step_elapsed(
            verbose,
            baseline_name,
            method_type,
            "label/identity merge",
            step_started_at,
        )

        # Step 3. 从同一份 merged_base 分别得到 all 和 enzyme/non-enzyme
        identity_paths = _identity_cache_paths(dir_path, baseline_name, False, avg_suffix)
        identity_enzyme_paths = _identity_cache_paths(dir_path, baseline_name, True, avg_suffix)
        identity_binary_std_file = (
            dir_path / f"{baseline_name}_identity_bins_enzyme_vs_non_enzyme{avg_suffix}_std_mean.feather"
        )
        identity_binary_metrics_file = (
            dir_path / f"{baseline_name}_identity_bins_enzyme_vs_non_enzyme{avg_suffix}_metrics.feather"
        )

        identity_std, identity_metrics, identity_no_pred = _evaluate_identity_from_merged(
            merged_base=merged_base,
            baseline_name=baseline_name,
            method_type=method_type,
            split_enzyme=False,
            num_folds=num_folds,
            max_workers=fold_workers,
            metric_n_jobs=metric_n_jobs,
            bin_workers=bin_workers,
            avg_types=None,
            cache_std_mean_file=identity_paths[0],
            cache_metrics_file=identity_paths[1],
            cache_no_pred_file=identity_paths[2],
            force=force,
            verbose=verbose,
        )
        identity_enzyme_std, identity_enzyme_metrics, identity_enzyme_no_pred = _evaluate_identity_from_merged(
            merged_base=merged_base,
            baseline_name=baseline_name,
            method_type=method_type,
            split_enzyme=True,
            num_folds=num_folds,
            max_workers=fold_workers,
            metric_n_jobs=metric_n_jobs,
            bin_workers=bin_workers,
            avg_types=None,
            cache_std_mean_file=identity_enzyme_paths[0],
            cache_metrics_file=identity_enzyme_paths[1],
            cache_no_pred_file=identity_enzyme_paths[2],
            force=force,
            verbose=verbose,
        )
        identity_binary_std, identity_binary_metrics = _evaluate_enzyme_binary_from_merged(
            merged_base=merged_base,
            baseline_name=baseline_name,
            method_type=method_type,
            num_folds=num_folds,
            max_workers=fold_workers,
            cache_std_mean_file=identity_binary_std_file,
            cache_metrics_file=identity_binary_metrics_file,
            force=force,
            verbose=verbose,
        )

        # Step 4. 组装完整大表，避免 notebook 侧再做拼接。
        overall_std = overall_std.copy()
        overall_std["identity_bin"] = "Overall"
        overall_std["bin_percentage"] = "1.0000 ± 0.0000"
        overall_std["bin_percentage_mean"] = 1.0
        overall_std["bin_percentage_std"] = 0.0
        overall_std["bin_percentage_overall"] = 1.0

        overall_metrics = overall_metrics.copy()
        overall_metrics["identity_bin"] = "Overall"
        overall_metrics["bin_percentage"] = 1.0

        if overall_no_pred is not None and not overall_no_pred.empty:
            overall_no_pred = overall_no_pred.copy()
            overall_no_pred["identity_bin"] = "Overall"
            overall_no_pred["bin_percentage"] = 1.0

        final_std_mean = pd.concat(
            [overall_std, identity_std, identity_enzyme_std, identity_binary_std],
            axis=0,
            ignore_index=True,
        )
        final_metrics = pd.concat(
            [overall_metrics, identity_metrics, identity_enzyme_metrics, identity_binary_metrics],
            axis=0,
            ignore_index=True,
        )

        no_pred_parts = [df for df in [overall_no_pred, identity_no_pred, identity_enzyme_no_pred] if df is not None and not df.empty]
        final_no_pred = (
            pd.concat(no_pred_parts, axis=0, ignore_index=True)
            if no_pred_parts
            else pd.DataFrame()
        )

        _save_identity_cache(
            final_std_mean=final_std_mean,
            final_metrics=final_metrics,
            final_no_pred=final_no_pred,
            cache_std_mean_file=combined_paths[0],
            cache_metrics_file=combined_paths[1],
            cache_no_pred_file=combined_paths[2],
            baseline_name=baseline_name,
            method_type=method_type,
            verbose=verbose,
        )
        return final_std_mean, final_metrics, final_no_pred if not final_no_pred.empty else pd.DataFrame()

    cache_std_mean_file, cache_metrics_file, cache_no_pred_file = _identity_cache_paths(
        dir_path=dir_path,
        baseline_name=baseline_name,
        split_enzyme=split_enzyme,
        avg_suffix=avg_suffix,
    )

    # 先检查 identity-bin 缓存。命中时直接返回，避免白做 label / pident / merge。
    if not force:
        cached = _load_identity_cache(
            baseline_name=baseline_name,
            method_type=method_type,
            cache_std_mean_file=cache_std_mean_file,
            cache_metrics_file=cache_metrics_file,
            cache_no_pred_file=cache_no_pred_file,
            verbose=verbose,
        )
        if cached is not None:
            return cached

    label_file = _ensure_identity_label_cache(
        baseline_name=baseline_name,
        dict_rxn2id=_load_dict_rxn2id(verbose=verbose),
        method_type=method_type,
        num_folds=num_folds,
        force=force,
        use_label_cache=use_label_cache,
        use_metrics_cache=use_metrics_cache,
        max_workers=fold_workers,
        metric_n_jobs=metric_n_jobs,
        avg_types=None,
        verbose=verbose,
    )

    _log_progress(verbose, f"{baseline_name}/{method_type}: loading label cache -> {label_file}")
    step_started_at = time.time()
    labeled = pd.read_feather(label_file)
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "label cache load",
        step_started_at,
    )

    step_started_at = time.time()
    valid_features = _load_validation_identity_features(
        num_folds=num_folds,
        verbose=verbose,
        baseline_name=baseline_name,
        method_type=method_type,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "validation identity feature load",
        step_started_at,
    )

    step_started_at = time.time()
    merged_base = _merge_labels_with_identity_bins(
        labeled=labeled,
        valid_features=valid_features,
        identity_bins=identity_bins,
        baseline_name=baseline_name,
        method_type=method_type,
        verbose=verbose,
    )
    _log_step_elapsed(
        verbose,
        baseline_name,
        method_type,
        "label/identity merge",
        step_started_at,
    )

    return _evaluate_identity_from_merged(
        merged_base=merged_base,
        baseline_name=baseline_name,
        method_type=method_type,
        split_enzyme=split_enzyme,
        num_folds=num_folds,
        max_workers=fold_workers,
        metric_n_jobs=metric_n_jobs,
        bin_workers=bin_workers,
        avg_types=None,
        cache_std_mean_file=cache_std_mean_file,
        cache_metrics_file=cache_metrics_file,
        cache_no_pred_file=cache_no_pred_file,
        force=force,
        verbose=verbose,
    )


def get_eval_results_full(
    baseline_name: str,
    method_type: str,
    num_folds: int = 10,
    identity_bins: list[float] | None = None,
    force: bool = False,
    use_label_cache: bool = True,
    use_metrics_cache: bool = True,
    max_workers: int = 10,
    bin_workers: int = 10,
    verbose: bool = True,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]]:
    """一次性产出论文需要的三类结果：overall / identity / identity_enzyme。"""

    overall = get_eval_results(
        baseline_name=baseline_name,
        method_type=method_type,
        num_folds=num_folds,
        force=force,
        use_label_cache=use_label_cache,
        use_metrics_cache=use_metrics_cache,
        max_workers=max_workers,
        verbose=verbose,
    )

    if identity_bins is None:
        identity_bins = [0.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

    label_file = _ensure_label_cache(
        baseline_name=baseline_name,
        dict_rxn2id=_load_dict_rxn2id(verbose=verbose),
        method_type=method_type,
        num_folds=num_folds,
        use_label_cache=True,
        verbose=verbose,
    )
    _log_progress(verbose, f"{baseline_name}/{method_type}: loading label cache -> {label_file}")
    labeled = pd.read_feather(label_file)
    valid_features = _load_validation_identity_features(
        num_folds=num_folds,
        verbose=verbose,
        baseline_name=baseline_name,
        method_type=method_type,
    )
    merged_base = _merge_labels_with_identity_bins(
        labeled=labeled,
        valid_features=valid_features,
        identity_bins=identity_bins,
        baseline_name=baseline_name,
        method_type=method_type,
        verbose=verbose,
    )

    dir_path = _method_output_dir(method_type)
    avg_types = None
    avg_suffix = ""
    identity_paths = _identity_cache_paths(dir_path, baseline_name, False, avg_suffix)
    identity_enzyme_paths = _identity_cache_paths(dir_path, baseline_name, True, avg_suffix)

    identity = _evaluate_identity_from_merged(
        merged_base=merged_base,
        baseline_name=baseline_name,
        method_type=method_type,
        split_enzyme=False,
        num_folds=num_folds,
        max_workers=max_workers,
        bin_workers=bin_workers,
        avg_types=None,
        cache_std_mean_file=identity_paths[0],
        cache_metrics_file=identity_paths[1],
        cache_no_pred_file=identity_paths[2],
        force=force,
        verbose=verbose,
    )
    identity_enzyme = _evaluate_identity_from_merged(
        merged_base=merged_base,
        baseline_name=baseline_name,
        method_type=method_type,
        split_enzyme=True,
        num_folds=num_folds,
        max_workers=max_workers,
        bin_workers=bin_workers,
        avg_types=None,
        cache_std_mean_file=identity_enzyme_paths[0],
        cache_metrics_file=identity_enzyme_paths[1],
        cache_no_pred_file=identity_enzyme_paths[2],
        force=force,
        verbose=verbose,
    )

    return {
        "overall": overall,
        "identity": identity,
        "identity_enzyme": identity_enzyme,
    }


def get_eval_results_full_expanded(
    baseline_name: str,
    method_type: str,
    avg_type: str = "micro",
    num_folds: int = 10,
    identity_bins: list[float] | None = None,
    force: bool = False,
    use_label_cache: bool = True,
    use_metrics_cache: bool = True,
    max_workers: int = 10,
    bin_workers: int = 10,
    verbose: bool = True,
) -> tuple[pd.DataFrame, dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]]]:
    """一次性返回一个方法在论文里需要的完整大表。

    返回值包含两部分：
    1. expanded_table：已经合并好的展示表，包含
       - Overall(all)
       - 每个 identity bin 的 all
       - Overall_enzyme / Overall_non_enzyme
       - 每个 identity bin 的 enzyme / non_enzyme
    2. raw_results：底层原始结果，结构与 `get_eval_results_full()` 一致，便于继续排查。

    这个函数的目的就是把 notebook 里常见的两张表：
    - 不分 enzyme 的 identity 表
    - 分 enzyme/non-enzyme 的 identity 表
    直接合并成一个更直观的大表。
    """
    raw_results = get_eval_results_full(
        baseline_name=baseline_name,
        method_type=method_type,
        num_folds=num_folds,
        identity_bins=identity_bins,
        force=force,
        use_label_cache=use_label_cache,
        use_metrics_cache=use_metrics_cache,
        max_workers=max_workers,
        bin_workers=bin_workers,
        verbose=verbose,
    )

    overall_std, _, _ = raw_results["overall"]
    identity_std, _, _ = raw_results["identity"]
    identity_enzyme_std, _, _ = raw_results["identity_enzyme"]

    combined_std = pd.concat(
        [
            identity_std[identity_std["avgType"] == avg_type],
            identity_enzyme_std[identity_enzyme_std["avgType"] == avg_type],
        ],
        axis=0,
        ignore_index=True,
    )

    expanded_table = expand_binned_table(
        combined_std,
        std_overall=overall_std[overall_std["avgType"] == avg_type],
        method_type=method_type or "direct",
    )
    return expanded_table, raw_results


def load_method_table(
    baseline_name: str,
    method_type: str,
    avg_type: str = "micro",
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None, pd.DataFrame]:
    """加载单个方法的完整指标与展示表。

    返回：
    1. std_df：按 bin 聚合后的统计表；
    2. metrics_df：逐 fold 的原始指标；
    3. no_pred_df：no prediction / no reaction 统计；
    4. display_df：展开后的展示表，适合 notebook 直接查看。
    """
    std_df, metrics_df, no_pred_df = get_eval_results_by_identity(
        baseline_name=baseline_name,
        method_type=method_type,
        **kwargs,
    )
    display_df = expand_binned_table(
        std_df[std_df["avgType"] == avg_type],
        method_type=method_type,
    )
    return std_df, metrics_df, no_pred_df, display_df


def expand_binned_table(aa: pd.DataFrame, std_overall: pd.DataFrame | None = None, method_type: str = 'direct') -> pd.DataFrame:
    """
    Expand binned cross-validation results horizontally and append an 'Overall' row.
    
    Parameters:
    -----------
    aa : pd.DataFrame
        The binned results DataFrame (e.g., std_blast[std_blast.avgType == 'macro'])
    std_overall : pd.DataFrame, optional
        The overall (unbinned) metrics DataFrame (e.g. from get_eval_results).
        If None, the function will try to automatically load and compute it from disk.
    method_type : str, default 'direct'
        The method type ('direct' or 'ec') to locate overall metrics files on disk.
    """
    aa = aa.copy()
    
    # 1. Obtain or load overall metrics
    overall_list = []
    
    has_overall_rows = (
        'identity_bin' in aa.columns
        and aa['identity_bin'].astype(str).str.startswith('Overall').any()
    )

    if not has_overall_rows:
        if std_overall is not None:
            avg_types = aa['avgType'].unique()
            df_overall = std_overall[std_overall['avgType'].isin(avg_types)].copy()
            if not df_overall.empty:
                df_overall['identity_bin'] = 'Overall'
                df_overall['bin_percentage'] = 1.0
                overall_list.append(df_overall)
        else:
            for baseline in aa['baselineName'].unique():
                dir_path = _method_output_dir(method_type)
                metrics_candidates = [
                    dir_path / f"{baseline}_10_folds_metrics_res.feather",
                ]

                loaded_overall = False
                for metrics_file in metrics_candidates:
                    if not metrics_file.exists():
                        continue
                    try:
                        df_metrics = pd.read_feather(metrics_file)
                        avg_types = aa['avgType'].unique()
                        df_metrics_filtered = df_metrics[df_metrics['avgType'].isin(avg_types)]
                        
                        if not df_metrics_filtered.empty:
                            df_overall_std = eva.get_fold_mean_std_metrics(df_metrics_filtered)
                            df_overall_std['identity_bin'] = 'Overall'
                            df_overall_std['bin_percentage'] = 1.0
                            overall_list.append(df_overall_std)
                            loaded_overall = True
                            break
                    except Exception as e:
                        print(f"Warning: Failed to load overall metrics for {baseline}: {e}")
                if not loaded_overall:
                    print(f"Warning: Overall metrics file not found for {baseline} in {dir_path}")
                
    # Append overall rows to aa if any were found
    if overall_list:
        df_overall_all = pd.concat(overall_list, axis=0, ignore_index=True)
        # Align columns
        for col in aa.columns:
            if col not in df_overall_all.columns:
                df_overall_all[col] = method_type if col == 'methodGroup' else None
        aa = pd.concat([aa, df_overall_all], axis=0, ignore_index=True)

    # 2. Identify Index Columns
    index_cols = [c for c in aa.columns if c not in ['Metric', 'mean', 'std']]
    
    if 'methodGroup' not in aa.columns:
        aa['methodGroup'] = method_type
        if 'methodGroup' not in index_cols:
            index_cols.append('methodGroup')

    # 3. Pivot mean and std
    df_mean = aa.pivot(index=index_cols, columns='Metric', values='mean')
    df_mean.columns = [f"{c}_mean" for c in df_mean.columns]

    df_std = aa.pivot(index=index_cols, columns='Metric', values='std')
    df_std.columns = [f"{c}_std" for c in df_std.columns]

    # 4. Concatenate and format
    res_pivoted = pd.concat([df_mean, df_std], axis=1)
    
    metrics = ['mAccuracy', 'mPrecision', 'mRecall', 'mF1']
    for metric in metrics:
        mean_col = f"{metric}_mean"
        std_col = f"{metric}_std"
        if mean_col in res_pivoted.columns and std_col in res_pivoted.columns:
            res_pivoted[metric] = res_pivoted.apply(
                lambda row: f"{row[mean_col]:.4f} ± {row[std_col]:.4f}" 
                if pd.notna(row[mean_col]) and pd.notna(row[std_col]) else "N/A", 
                axis=1
            )

    # 5. Reset index and reorder
    res_pivoted = res_pivoted.reset_index()

    cols_order = []
    for metric in metrics:
        cols_order.extend([f"{metric}_mean", f"{metric}_std"])
    cols_order.extend(metrics)

    final_cols = index_cols + [c for c in cols_order if c in res_pivoted.columns]
    res_expanded = res_pivoted[final_cols]

    # 把 `identity_bin` 拆成两列，便于论文表格阅读：
    # - bin：基础 identity bin，例如 Overall / (90,100] / No BLAST Hit
    # - enzymeGroup：all / enzyme / non_enzyme / enzyme_vs_non_enzyme
    # Sort: Overall -> Overall_enzyme -> Overall_non_enzyme -> Overall_enzyme_vs_non_enzyme ->
    #       90/80/... bins with各 group 挨在一起 ->
    #       No BLAST Hit_* 最后
    if 'identity_bin' in res_expanded.columns:
        res_expanded = res_expanded.copy()

        def split_identity_bin(bin_val: str) -> tuple[str, str]:
            bin_str = str(bin_val)
            if bin_str.endswith('_enzyme_vs_non_enzyme'):
                return bin_str[:-21], 'enzyme_vs_non_enzyme'
            if bin_str.endswith('_non_enzyme'):
                return bin_str[:-11], 'non_enzyme'
            if bin_str.endswith('_enzyme'):
                return bin_str[:-7], 'enzyme'
            return bin_str, 'all'

        split_cols = res_expanded['identity_bin'].apply(split_identity_bin)
        res_expanded['bin'] = split_cols.apply(lambda x: x[0])
        res_expanded['enzymeGroup'] = split_cols.apply(lambda x: x[1])

        def get_bin_sort_key(bin_val):
            """把 identity_bin 映射成可读顺序。

            排序规则分 3 层：
            1. 大类顺序：Overall -> 普通 identity bins -> No BLAST Hit -> 其他；
            2. 普通 bins 内部：按 lower-bound 从高到低；
            3. 同一个 base bin 内：无后缀 -> _enzyme -> _non_enzyme。

            这样展示时每一组会挨在一起，更适合论文表格阅读。
            """
            bin_str = str(bin_val)

            if bin_str == 'Overall':
                return (0, 0, 0)
            if bin_str == 'Overall_enzyme':
                return (0, 1, 0)
            if bin_str == 'Overall_non_enzyme':
                return (0, 2, 0)
            if bin_str == 'Overall_enzyme_vs_non_enzyme':
                return (0, 3, 0)

            if bin_str == 'No BLAST Hit':
                return (2, 0, 0)
            if bin_str == 'No BLAST Hit_enzyme':
                return (2, 0, 1)
            if bin_str == 'No BLAST Hit_non_enzyme':
                return (2, 0, 2)
            if bin_str == 'No BLAST Hit_enzyme_vs_non_enzyme':
                return (2, 0, 3)

            suffix_rank = 0
            base_str = bin_str
            if bin_str.endswith('_enzyme_vs_non_enzyme'):
                suffix_rank = 3
                base_str = bin_str[:-21]
            elif bin_str.endswith('_non_enzyme'):
                suffix_rank = 2
                base_str = bin_str[:-11]
            elif bin_str.endswith('_enzyme'):
                suffix_rank = 1
                base_str = bin_str[:-7]

            try:
                cleaned = (
                    base_str
                    .replace('(', '')
                    .replace(')', '')
                    .replace('[', '')
                    .replace(']', '')
                )
                parts = [float(x.strip()) for x in cleaned.split(',')]
                if len(parts) == 2:
                    lower_bound = parts[0]
                    return (1, -lower_bound, suffix_rank)
            except Exception:
                pass

            return (9, 0, suffix_rank)

        res_expanded['_sort_key'] = res_expanded['identity_bin'].apply(get_bin_sort_key)
        res_expanded = (
            res_expanded
            .sort_values(by=['baselineName', '_sort_key'])
            .drop(columns=['_sort_key'])
            .reset_index(drop=True)
        )

        # 调整列顺序：第一列 bin，第二列 enzymeGroup，identity_bin 原始列移除。
        metric_cols = [c for c in res_expanded.columns if c not in {'identity_bin', 'bin', 'enzymeGroup'}]
        if 'bin_percentage' in metric_cols:
            metric_cols.remove('bin_percentage')
            leading_percentage_cols = ['bin_percentage']
            for extra_col in ['bin_percentage_mean', 'bin_percentage_std', 'bin_percentage_overall']:
                if extra_col in metric_cols:
                    metric_cols.remove(extra_col)
                    leading_percentage_cols.append(extra_col)
            metric_cols = leading_percentage_cols + metric_cols
        if 'baselineName' in metric_cols:
            metric_cols.remove('baselineName')
            metric_cols = ['baselineName'] + metric_cols
        if 'avgType' in metric_cols:
            metric_cols.remove('avgType')
            insert_pos = 2 if 'baselineName' in metric_cols and 'bin_percentage' in metric_cols else len(metric_cols)
            metric_cols.insert(insert_pos, 'avgType')
        if 'methodGroup' in metric_cols:
            metric_cols.remove('methodGroup')
            insert_pos = metric_cols.index('avgType') + 1 if 'avgType' in metric_cols else len(metric_cols)
            metric_cols.insert(insert_pos, 'methodGroup')

        leading_cols = []
        if 'baselineName' in metric_cols:
            leading_cols.append('baselineName')
            metric_cols.remove('baselineName')
        leading_cols.extend(['bin', 'enzymeGroup'])
        if 'bin_percentage' in metric_cols:
            leading_cols.append('bin_percentage')
            metric_cols.remove('bin_percentage')
        for extra_col in ['bin_percentage_mean', 'bin_percentage_std', 'bin_percentage_overall']:
            if extra_col in metric_cols:
                leading_cols.append(extra_col)
                metric_cols.remove(extra_col)
        if 'avgType' in metric_cols:
            leading_cols.append('avgType')
            metric_cols.remove('avgType')
        if 'methodGroup' in metric_cols:
            leading_cols.append('methodGroup')
            metric_cols.remove('methodGroup')

        res_expanded = res_expanded[leading_cols + metric_cols]

    return res_expanded


def get_combined_enzyme_binned_table(
    avg_type: str = "micro",
    save_dir: str | Path | None = None,
    force_recompute: bool = False,
    bin_workers: int | None = 9,
    num_folds: int = 10,
    direct_baselines: list[str] | None = None,
    ec_baselines: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, pd.DataFrame]:
    """统一计算并保存按 enzyme/non-enzyme 拆分的 identity-bin 指标大表。

    这个函数用于 notebook 末尾的汇总步骤，不修改已经校验过的“不分酶/非酶”
    结果。它会对 direct 方法和 EC 方法分别调用：

        get_eval_results_by_identity(..., split_enzyme=True)

    得到每个方法在如下分组中的 10-fold 指标：
    - 每个 identity bin 的 `_enzyme` 和 `_non_enzyme` 子集；
    - `No BLAST Hit_enzyme` 和 `No BLAST Hit_non_enzyme`；
    - `Overall_enzyme` 和 `Overall_non_enzyme`。

    保存结果：
    - direct_methods_enzyme_split_{avg_type}.feather/csv/tsv
    - ec_methods_enzyme_split_{avg_type}.feather/csv/tsv
    - all_methods_enzyme_split_{avg_type}.feather/csv/tsv

    返回一个 dict，包含 `direct`、`ec`、`all` 三个 DataFrame，方便 notebook
    继续展示、核对或绘图。
    """


    if direct_baselines is None:
        direct_baselines = [
            "blast",
            "unirep_euclidean",
            "unirep_cosine",
            "esm_euclidean",
            "esm_cosine",
            "t5_euclidean",
            "t5_cosine",
            "RXNRECer",
            "tdit5_euclidean",
            "tdit5_cosine",
        ]
    if ec_baselines is None:
        ec_baselines = ["blast", "deepec", "clean", "ecrecer", "catfam", "priam"]

    if save_dir is None:
        save_dir = Path(cfg.DIR_PROJECT_ROOT) / "benchmarks" / "evaluation" / "NC_R1"
    else:
        save_dir = Path(save_dir)
    os.makedirs(save_dir, exist_ok=True)

    def evaluate_group(baselines: list[str], method_type: str) -> pd.DataFrame:
        expanded_tables = []
        for baseline in baselines:
            _log_progress(
                verbose,
                f"{baseline}/{method_type}: enzyme-split identity-bin summary started",
            )
            try:
                std_df, _, _ = get_eval_results_by_identity(
                    baseline_name=baseline,
                    method_type=method_type,
                    num_folds=num_folds,
                    split_enzyme=True,
                    force=force_recompute,
                    bin_workers=bin_workers,
                    verbose=verbose,
                )
            except Exception as e:
                print(
                    f"Warning: Failed to compute enzyme-split identity bins for "
                    f"{baseline}/{method_type}: {e}"
                )
                continue

            df_filtered = std_df[std_df.avgType == avg_type].copy()
            if df_filtered.empty:
                print(f"Warning: No {avg_type} rows for {baseline}/{method_type}")
                continue

            expanded_tables.append(expand_binned_table(df_filtered, method_type=method_type))

        if not expanded_tables:
            return pd.DataFrame()
        return pd.concat(expanded_tables, axis=0).reset_index(drop=True)

    direct_table = evaluate_group(direct_baselines, "direct")
    ec_table = evaluate_group(ec_baselines, "ec")
    all_parts = [df for df in [direct_table, ec_table] if not df.empty]
    all_table = pd.concat(all_parts, axis=0).reset_index(drop=True) if all_parts else pd.DataFrame()

    outputs = {
        "direct": direct_table,
        "ec": ec_table,
        "all": all_table,
    }
    for name, df in outputs.items():
        if df.empty:
            continue
        feather_file = save_dir / f"{name}_methods_enzyme_split_{avg_type}.feather"
        csv_file = save_dir / f"{name}_methods_enzyme_split_{avg_type}.csv"
        tsv_file = save_dir / f"{name}_methods_enzyme_split_{avg_type}.tsv"
        df.to_feather(feather_file)
        df.to_csv(csv_file, index=False)
        df.to_csv(tsv_file, sep="\t", index=False)
        _log_progress(verbose, f"saved {name} enzyme-split table -> {feather_file}")

    return outputs


def get_combined_binned_table(
    avg_type: str = "macro",
    save_dir: str | Path | None = None,
    force_recompute: bool = False,
) -> pd.DataFrame:
    """Generate a combined large table of identity-binned metrics without enzyme splitting.

    This function aggregates cross-validation results binned by sequence identity
    for all 10 direct methods and 6 EC methods, formats them using expand_binned_table,
    renames the 'blast' baseline based on its group to avoid duplication, and sorts
    them in a logical order (EC methods first, then direct methods, with overall/bins descending).

    Parameters:
    -----------
    avg_type : str, default 'macro'
        The averaging method to filter by ('macro', 'micro', 'weighted', 'samples').
    save_dir : str or Path, optional
        The directory to save the output files. If None, saves to the results directory.
    force_recompute : bool, default False
        Whether to force recomputation of the binned results.

    Returns:
    --------
    pd.DataFrame
        The combined and sorted binned results table.
    """
    direct_baselines = [
        "RXNRECer", "blast", "esm_cosine", "esm_euclidean", 
        "t5_cosine", "t5_euclidean", "tdit5_cosine", "tdit5_euclidean", 
        "unirep_cosine", "unirep_euclidean"
    ]
    ec_baselines = [
        "blast", "catfam", "clean", "deepec", "ecrecer", "priam"
    ]

    all_expanded = []

    # Process Direct methods
    for b in direct_baselines:
        p_dir = _method_output_dir("direct")
        p = p_dir / f"{b}_identity_bins_std_mean.feather"
        if force_recompute or not p.exists():
            try:
                get_eval_results_by_identity(
                    baseline_name=b,
                    method_type="direct",
                    num_folds=10,
                    force=force_recompute,
                    split_enzyme=False,
                )
            except Exception as e:
                print(f"Warning: Failed to compute identity bins for direct baseline {b}: {e}")
                continue

        if p.exists():
            df = pd.read_feather(p)
            df_filtered = df[df.avgType == avg_type]
            if not df_filtered.empty:
                expanded = expand_binned_table(df_filtered, method_type="direct")
                all_expanded.append(expanded)

    # Process EC methods
    for b in ec_baselines:
        p_dir = _method_output_dir("ec")
        p = p_dir / f"{b}_identity_bins_std_mean.feather"
        if force_recompute or not p.exists():
            try:
                get_eval_results_by_identity(
                    baseline_name=b,
                    method_type="ec",
                    num_folds=10,
                    force=force_recompute,
                    split_enzyme=False,
                )
            except Exception as e:
                print(f"Warning: Failed to compute identity bins for EC baseline {b}: {e}")
                continue

        if p.exists():
            df = pd.read_feather(p)
            df_filtered = df[df.avgType == avg_type]
            if not df_filtered.empty:
                expanded = expand_binned_table(df_filtered, method_type="ec")
                all_expanded.append(expanded)

    if not all_expanded:
        raise ValueError(f"No binned data found or computed for avg_type: {avg_type}")

    combined = pd.concat(all_expanded, axis=0, ignore_index=True)

    # Rename blast to distinguish direct vs EC
    combined.loc[(combined.baselineName == "blast") & (combined.methodGroup == "direct"), "baselineName"] = "blast_via_rxn"
    combined.loc[(combined.baselineName == "blast") & (combined.methodGroup == "ec"), "baselineName"] = "blast_via_ec"

    # Define a logical sorting order for methods
    preferred_order = [
        "priam", "deepec", "clean", "catfam", "blast_via_ec", "ecrecer",
        "blast_via_rxn", "unirep_euclidean", "unirep_cosine",
        "esm_euclidean", "esm_cosine", "t5_euclidean", "t5_cosine",
        "tdit5_euclidean", "tdit5_cosine", "RXNRECer"
    ]
    method_order_dict = {name: idx for idx, name in enumerate(preferred_order)}

    # Sort key function for sorting methods logically and sequence identity bins descending
    def get_sort_keys(row):
        method = row["baselineName"]
        method_key = method_order_dict.get(method, 999)

        bin_str = str(row["bin"])
        enzyme_group = str(row.get("enzymeGroup", "all"))
        if bin_str == "Overall":
            bin_key = -9999.0
        elif bin_str in ("noPrediction", "No BLAST Hit"):
            bin_key = 9999.0
        else:
            try:
                cleaned = bin_str.replace("(", "").replace("]", "").replace("[", "").replace(")", "")
                parts = [float(x.strip()) for x in cleaned.split(",")]
                if len(parts) == 2:
                    bin_key = -parts[0]
                else:
                    bin_key = 500.0
            except Exception:
                bin_key = 500.0

        enzyme_key = {"all": 0, "enzyme": 1, "non_enzyme": 2}.get(enzyme_group, 9)
        return (method_key, bin_key, enzyme_key)

    keys = combined.apply(get_sort_keys, axis=1)
    combined["_sort_method_idx"] = [k[0] for k in keys]
    combined["_sort_bin_val"] = [k[1] for k in keys]
    combined["_sort_enzyme_group"] = [k[2] for k in keys]

    combined = combined.sort_values(
        by=["_sort_method_idx", "_sort_bin_val", "_sort_enzyme_group"]
    ).reset_index(drop=True)
    combined = combined.drop(columns=["_sort_method_idx", "_sort_bin_val", "_sort_enzyme_group"])

    # Determine save directory
    if save_dir is None:
        save_dir = Path(cfg.RESULTS_DIR)
    else:
        save_dir = Path(save_dir)

    os.makedirs(save_dir, exist_ok=True)
    
    # Save combined table to CSV and TSV
    csv_file = save_dir / f"combined_identity_bins_no_enzyme_{avg_type}.csv"
    tsv_file = save_dir / f"combined_identity_bins_no_enzyme_{avg_type}.tsv"
    
    combined.to_csv(csv_file, index=False)
    combined.to_csv(tsv_file, sep="\t", index=False)
    print(f"Successfully saved combined table to {csv_file} and {tsv_file}")

    return combined
