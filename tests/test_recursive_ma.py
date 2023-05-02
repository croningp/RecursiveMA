import pickle
from pathlib import Path

import pytest

from recursive_ma import build_tree, estimate_MA, estimate_by_MW, common_MA

HERE = Path(__file__).parent


@pytest.fixture
def pickle_files():
    data_dict = {
        # filename: 'Sample..._ms<n>.pkl' => n
        int(pickle_file.name.split("_")[1][2]): pickle.load(pickle_file.open("rb"))
        for pickle_file in sorted(HERE.glob("*.pkl"))
    }

    # rename ms<n>_mz columns to mz and ms<n>_intensity to intensity
    for level, data in data_dict.items():
        data.rename(
            columns={f"ms{level}_mz": "mz", f"ms{level}_intensity": "intensity"},
            inplace=True,
        )

    return data_dict


@pytest.fixture
def test_data(pickle_files, request):
    level = request.param
    return build_tree(pickle_files, max_level=level)


@pytest.fixture
def mock_data():
    # Mock MS3 data (could go higher).
    # {} means ion did not fragment; None means didn't try fragmenting ion.
    return {371.2: {150.1: {72.3: None, 89.1: None}, 221.3: {72.3: None, 99.7: None}}}


def test_mock_data(mock_data):
    parent_ma = estimate_MA(mock_data, 371.2)
    child1_ma = estimate_MA(mock_data[371.2], 150.1)
    child2_ma = estimate_MA(mock_data[371.2], 221.3)
    common_ma = common_MA(mock_data[371.2], 150.1, 221.3)
    assert parent_ma <= child1_ma + child2_ma


# test_data parameterised on ms_level
@pytest.mark.parametrize("test_data", [1, 2, 3, 4], indirect=True)
def test_real_data(test_data, request):
    mw = list(test_data)[0]
    ma = estimate_MA(test_data, mw)
    assert estimate_by_MW(mw) >= ma > 0.0
