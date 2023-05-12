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
