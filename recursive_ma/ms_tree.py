import logging
from collections import defaultdict

import pandas as pd


def _build_tree(data, level=1, acc=None, parent=None, max_level=3):
    max_level = min(max_level, max(data))
    if level == 1:
        acc = {}
        for peak in data[1]["mz"]:
            acc[peak] = {}
            _build_tree(data, level=2, acc=acc[peak], parent=peak, max_level=max_level)
        return acc
    if level > max_level:
        return
    level_df = data[level]
    parent_df = data[level - 1].drop(columns=["parent_id"], errors="ignore")
    level_df = level_df.join(parent_df, on="parent_id", rsuffix="_parent")
    child_peaks = level_df[level_df["mz_parent"] == parent]["mz"].unique()
    for peak in child_peaks:
        acc[peak] = None if level == max_level else {}
        _build_tree(
            data, level=level + 1, acc=acc[peak], parent=peak, max_level=max_level
        )


def build_tree(data: dict, max_level=3):
    """
    Build a tree from a dictionary of {level: DataFrame, ...}.

    Args:
        data (dict): {level: DataFrame, ...}.
        max_level (int): maximum MS level to include in the tree, e.g. 3 for MS3.
    """
    return _build_tree(data, max_level=max_level)


def tree_depth(tree: dict):
    """
    Calculate the level of nesting in a tree.

    >>> tree_depth({})
    0
    >>> tree_depth({1: None})
    1
    >>> tree_depth({1: {}})
    1
    >>> tree_depth({1: {2: {3: None}}})
    3
    >>> tree_depth({1: {2: {3: None}, 4: None}})
    3
    >>> tree_depth({1: {2: {3: None}, 4: {5: None}}})
    4
    """
    if isinstance(tree, dict) and len(tree) > 0:
        return 1 + max(tree_depth(v) for v in tree.values())
    else:
        return 0


def _process_df(
    level,
    ms_df: pd.DataFrame,
    max_num_peaks,
    min_abs_intensity,
    min_rel_intensity,
    n_digits,
):
    original_len = len(ms_df)

    if "parent" not in ms_df:
        ms_df["parent"] = 10**6
    ms_df = ms_df[ms_df.mz < ms_df.parent - 1]

    ms_df = ms_df.assign(
        mz_bin=(ms_df.mz.round(n_digits) * 10**n_digits).astype(int),
        parent_bin=(ms_df.parent.round(n_digits) * 10**n_digits).astype(int),
    )

    min_intensity_fn = lambda df: max(
        min_abs_intensity[level], df["intensity"].max() * min_rel_intensity
    )
    filter_fn = (
        lambda g: g[g["intensity"] > min_intensity_fn(g)]
        .sort_values("intensity")
        .tail(max_num_peaks)
    )
    ms_df = (
        ms_df.groupby(["mz_bin", "parent_bin"])
        .agg({"intensity": "sum", "mz": "median", "parent": "median"})
        .reset_index()
        .groupby("parent_bin")
        .apply(filter_fn)
        .reset_index()
    )

    result = (
        ms_df.groupby("parent_rounded")
        .apply(lambda g: g.sort_values("intensity").tail(max_num_peaks))
        .reset_index(drop=True)
    )

    logging.debug(f"Level {level}: {len(result)} out of {original_len} peaks retained")
    return result


def process(
    sample: dict[int, pd.DataFrame],
    max_num_peaks: int = 200,
    min_abs_intensity: dict[int, float] = defaultdict(lambda: 0.0),
    min_rel_intensity: float = 0.0,
    n_digits: int = 3,
) -> dict[int, pd.DataFrame]:
    """
    Process a sample by binning peak intensities and removing low-intensity peaks.

    Args:
        sample (dict[int, pd.DataFrame]): {level: DataFrame, ...}.
        max_num_peaks (int): maximum number of peaks to keep per parent.
        min_abs_intensity (dict[int, float]): minimum absolute child peak intensity by level.
        min_rel_intensity (float): minimum relative child peak intensity, relative to parent's most intense child peak.
        n_digits (int): Number of digits to round m/z values for intensity binning.
            Full precision retained in final output by taking median m/z.

    Returns:
        dict[int, pd.DataFrame]: {level: DataFrame, ...}.
    """
    sample = {
        level: _process_df(
            level,
            df.reset_index(),
            max_num_peaks,
            min_abs_intensity,
            min_rel_intensity,
            n_digits,
        )
        for level, df in sample.items()
    }
    if 1 not in sample:
        # Generate placeholder MS1 if only MS2+ present
        parent_peak = (
            sample[2].groupby("parent")["intensity"].sum().sort_values().index[-1]
        )
        sample[1] = pd.DataFrame(
            {
                "mz": [parent_peak],
                "intensity": 100000.0,
                "mz_bin": [int(parent_peak * 10**n_digits)],
            }
        )
    return sample


def identify_parents(dataset, mass_tol: float, ms_n_digits: int):
    """
    Link MSn+1 peaks to MSn peaks.

    Args:
        dataset (dict[int, pd.DataFrame]): {level: DataFrame, ...}.
            Each DataFrame must have columns "mz_bin" and "parent_bin", created by the process() function.
            These values are used to link parent and child peaks.
        mass_tol (float): mass tolerance in Da.
        ms_n_digits (int): number of decimal places used in binning m/z values.
            **The same value used in process() must be specified here.**
    """
    new_dataset = {}
    new_dataset[min(dataset)] = dataset[min(dataset)]
    for level in sorted(dataset)[:-1]:
        new_dataset[level + 1] = (
            pd.merge_asof(
                dataset[level + 1].sort_values("parent_bin"),
                new_dataset[level][["mz_bin"]]
                .sort_values("mz_bin")
                .reset_index(),
                left_on="parent_bin",
                right_on="mz_bin",
                suffixes=("", "_x"),
                tolerance=int(MASS_TOL * 10**MS_N_DIGITS),
                direction="nearest",
            )
            .rename(columns={"index": "parent_id"})
            .dropna(subset=["parent_id"])
            .astype({"parent_id": int})
            .drop(columns=["mz_bin_x"])
        )
    return new_dataset
