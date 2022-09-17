""" Tests for the nimgen/stats.py module"""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import numpy as np
import nimgen


def test_empirical_pval():

    stat_one = [1, 5, 3, 7, 8, 10, 11]
    stat_two = [4, 5, 1, 8, 1, 19, 12]

    pval = nimgen.statistics.empirical_pval(stat_one, stat_two)

    # empirical p-val is construed to be 0.375 --> * 1000 == 375
    assert int(pval * 1000) == 375
    
    stat_one = np.array([
        [1, 5, 3, 7, 8, 10, 11],
        [4, 5, 1, 8, 1, 19, 12]
    ]).T

    stat_two = np.array([
        [3, 1, 4, 5, 6, 1, 1],
        [2, 3, 1, 87, 12, 190, 2]
    ]).T

    pvals = nimgen.statistics.empirical_pval(stat_one, stat_two)
    assert len(pvals) == 2