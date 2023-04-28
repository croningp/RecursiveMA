# Recursive MS2MA
## Installation
```sh
# clone the repo then
pip install -e recursive_ma
```

## Usage
First convert your MSn data to the tree format:
```python
{ms1_mz1: {ms2_mz1: {ms3_mz1: {}, ms3_mz2: {}}, ms2_mz2: {ms3_mz3: {}}}, ...}
```

There is a utility function to convert data of the form `{ms_level: pd.DataFrame}` to the tree format:
```python
from recursive_ma import build_tree
tree = build_tree(data, max_level=<n>)
```

Note that this function assume that each dataframe has the columns `mz` and `intensity`, as well as `parent` for MS2 and higher.

Then you can run the recursive MS2MA algorithm:
```python
from recursive_ma import estimate_MA
# parent_mz is the m/z of the parent (MS1) ion
estimate_MA(tree, parent_mz)
```