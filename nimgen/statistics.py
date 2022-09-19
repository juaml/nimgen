"""Perform domain general statistical operations."""

from functools import partial

import numpy as np
from scipy.stats.mstats import winsorize

from .utils import logger


def empirical_pval(stat_null, stat):
    """Calculate empirical p-value.

    Calculate empirical p-value based on the observed (surrogate maps)
    and expected (reference map) correlation scores.

    Parameters
    ----------
    stat_null: numpy.array or list
        A vector or matrix of simulated or data-resampled null test statistics.
    stat: numpy.array or list
        A vector of calculated test statistics.

    Returns
    -------
    p-values : numpy.array
        Calculated empirical pvalues.
    """

    logger.info("Empirical p-value calculation...")
    check = np.sum(np.abs(stat_null) >= np.abs(stat), axis=0)
    pvalues = (check + 1) / (len(stat_null) + 1)
    return pvalues


def winsorized_mean(data, axis=None, **win_params):
    """Chain winsorization and mean to compute winsorize mean.

    Parameters
    ----------
    data : array
        Data to calculate winsorized mean on.
    win_params : dict
        Dictionary containing the keyword arguments for the winsorize function.
        E.g. {'limits': [0.1, 0.1]}

    Returns
    -------
    Winsorized mean of the inputted data with the winsorize settings applied
    as specified in win_params.

    """
    win_dat = winsorize(data, axis=axis, **win_params)
    win_mean = win_dat.mean(axis=axis)

    return win_mean


def _get_funcbyname(name, func_params):
    """Apply any function by name.

    Parameters
    ----------
    name : str
        Name to identify the function. Currently supported names and
        corresponding functions are:
        'winsorized_mean' -> scipy.stats.mstats.winsorize
        'mean' -> np.mean
        'std' -> np.std

    func_params : dict
        Dictionary containing functions that need further parameter
        specifications. Keys are the function and values are dictionaries
        with the parameter specifications.
        E.g. 'winsorized_mean': func_params = {'limits': [0.1, 0.1]}

    Returns
    -------
    respective function with inputted (partial) parameters.
    """

    # check validity of names
    _valid_func_names = {"winsorized_mean", "mean", "std"}

    # apply functions
    if name == "winsorized_mean":
        # check validity of func_params
        limits = func_params.get("limits")
        if all((lim >= 0.0 and lim <= 1) for lim in limits):
            logger.info(f"Limits for winsorized mean are set to {limits}.")
        else:
            raise ValueError(
                "Limits for the winsorized mean must be between 0 and 1."
            )
        # partially interpret func_params
        return partial(winsorized_mean, **func_params)
    if name == "mean":
        return np.mean  # No func_params
    if name == "std":
        return np.std
    if name == "median":
        return np.median

    else:
        raise ValueError(
            f"Function {name} unknown. Please provide any of "
            f"{_valid_func_names}"
        )
