"""Tests for the nimgen/stats.py module."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import numpy as np
import pytest

import nimgen


def test_empirical_pval():
    """Test empirical_pval."""
    stat_real = np.array([4, 5, 1, 8, 1, 19, 12, 2, 3])
    stat_null = np.arange(90).reshape(10, 9)
    pval = nimgen.statistics.empirical_pval(stat_null, stat_real)
    assert pval.shape[0] == 9


def test_winsorized_mean():
    """Test winsorized_mean."""
    data = np.array([[3, 1, 4, 5, 6, 1, 1], [2, 3, 1, 87, 12, 190, 2]]).T

    result = nimgen.statistics.winsorized_mean(data)
    assert isinstance(result, float)


def test_dcorr():
    """Test dcorr."""
    x = np.arange(9)
    y = np.flip(np.arange(9))
    dist_corr = nimgen.statistics.dcorr(x, y)
    assert isinstance(dist_corr[0], float)


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
