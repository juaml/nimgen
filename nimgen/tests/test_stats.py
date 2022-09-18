"""Tests for the nimgen/stats.py module."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import numpy as np
import pytest

import nimgen


def test_empirical_pval():
    """Test empirical_pval."""
    stat_one = [1, 5, 3, 7, 8, 10, 11]
    stat_two = [4, 5, 1, 8, 1, 19, 12]

    pval = nimgen.statistics.empirical_pval(stat_one, stat_two)

    # empirical p-val is construed to be 0.375 --> * 1000 == 375
    assert int(pval * 1000) == 375

    stat_one = np.array([[1, 5, 3, 7, 8, 10, 11], [4, 5, 1, 8, 1, 19, 12]]).T

    stat_two = np.array([[3, 1, 4, 5, 6, 1, 1], [2, 3, 1, 87, 12, 190, 2]]).T

    pvals = nimgen.statistics.empirical_pval(stat_one, stat_two)
    assert len(pvals) == 2


def test_winsorized_mean():
    """Test winsorized_mean."""
    data = np.array([[3, 1, 4, 5, 6, 1, 1], [2, 3, 1, 87, 12, 190, 2]]).T

    result = nimgen.statistics.winsorized_mean(data)
    assert isinstance(result, float)


def test_get_funcbyname_correct():
    """Test get_funcbyname with correct name."""
    data = np.array([[3, 1, 4, 5, 6, 1, 1], [2, 3, 1, 87, 12, 190, 2]]).T
    for name in ["mean", "winsorized_mean", "std", "median"]:
        if name == "winsorized_mean":
            param = {"limits": [0.1, 0.1]}
        else:
            param = {}
        func = nimgen.statistics._get_funcbyname(name, param)
        res = func(data)
        assert isinstance(res, float)


def test_get_funcbyname_incorrect():
    """Test get_funcbyname with incorrect name."""
    incorrect_name = "sjkfhjsdfk';alsd'fsdgfsd"
    with pytest.raises(ValueError, match=f"Function {incorrect_name} unknown"):
        nimgen.statistics._get_funcbyname(incorrect_name, {})
