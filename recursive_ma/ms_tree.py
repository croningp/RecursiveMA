
def _build_tree(data, level=1, acc=None, parent=None, max_level=3):
    max_level = min(max_level, max(data))
    if level == 1:
        acc = {}
        for peak in data[1]['mz']:
            acc[peak] = {}
            _build_tree(data, level=2, acc=acc[peak], parent=peak, max_level=max_level)
        return acc
    level_df = data[level]
    parent_df = data[level - 1].drop(columns=['parent_id'], errors='ignore')
    level_df = level_df.join(parent_df, on='parent_id', rsuffix='_parent')
    for peak in level_df[level_df[f'mz_parent'] == parent][f'mz']:
        acc[peak] = None if level == max_level else {}
        if level < max_level:
            _build_tree(data, level=level+1, acc=acc[peak], parent=peak, max_level=max_level)

def build_tree(data, max_level=3):
    return _build_tree(data, max_level=max_level)