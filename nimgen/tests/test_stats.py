"""Tests for the nimgen/stats.py module."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import numpy as np
import pytest

import nimgen


def test_empirical_pval():
    """Test empirical_pval."""
    stat_real = 5
    stat_null = [4, 5, 1, 8, 1, 19, 12, 2, 3]

    pval = nimgen.statistics.empirical_pval(stat_null, stat_real)
    # empirical p-val is construed to be 0.5 --> * 10 == 5
    assert int(pval * 10) == 5


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


test_empirical_pval()
