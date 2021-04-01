# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
#          Kaustubh Patil <k.patil@fz-juelich.de>
# License: AGPL
import numpy as np
import pandas as pd

from pathlib import Path
import os

from scipy import stats
from statsmodels.stats.multitest import multipletests

import nibabel as nib
import abagen
from .utils import logger


def _save_expressions(exp, atlas):
    logger.info(f'Trying to save expressions')
    save_path = atlas.parent
    atlas_name = atlas.stem
    if os.access(save_path, os.W_OK):
        exp_fname = save_path / f'nimgen_{atlas_name}_expressions.csv'
        logger.info(f'Saving expressions to {exp_fname.as_posix()}')
        exp.to_csv(exp_fname, sep=';', index=False)
    else:
        logger.info('User does not have permissions to save '
                    f'to {save_path.as_posix()}')


def _get_cached_results(atlas):
    exp_path = atlas.parent
    atlas_name = atlas.stem
    exp_fname = exp_path / f'nimgen_{atlas_name}_expressions.csv'
    exp = None
    if exp_fname.exists():
        logger.info(f'Reading expressions from {exp_fname.as_posix()}')

        exp = pd.read_csv(exp_fname, sep=';')
    return exp


def get_gene_expression(weights, atlas, allen_data_dir=None,
                        force_recompute=False, save_expressions=False,
                        multiple_correction='fdr_bh'):
    # WARNING: If this changes, then all the cached results
    # must be invalidated. TODO: Allow for multiple parameters
    # and save parameters values.
    abagen_params = {'probe_selection': 'diff_stability'}
    if not isinstance(weights, pd.DataFrame):
        weights = pd.DataFrame(weights)

    if not isinstance(atlas, Path):
        atlas = Path(atlas)

    if not atlas.exists():
        raise ValueError(f'Atlas file does not exist: {atlas.as_posix()}')

    logger.info('Checking atlas and weights dimensions')
    atlas_img = nib.load(atlas)
    nrois = np.unique(atlas_img.get_fdata()).astype(np.int).shape[0] - 1
    if nrois != len(weights):
        raise ValueError(f'Number of weights ({len(weights)}) does not match '
                         f'the number of ROIs in the atlas ({nrois}).')

    exp = None
    if force_recompute is False:
        # Check if we have this results
        exp = _get_cached_results(atlas)

    if exp is None or force_recompute is True:
        # If not, compute the expressions
        logger.info('Computing expressions, please wait')
        exp = abagen.get_expression_data(
            atlas, data_dir=allen_data_dir, **abagen_params)
        if save_expressions is True:
            _save_expressions(exp, atlas)

    # use only ROIs without NaNs (ROIs with samples)
    good_rois = ~exp.iloc[:, 0].isna().values
    exp = exp[good_rois]
    weights = weights[good_rois]
    pval = exp.apply(
        lambda ser: stats.pearsonr(ser.values, weights.squeeze().values)[1])

    reject, *_ = multipletests(pval, method=multiple_correction)
    genes = pval[reject].index.values.tolist()
    return genes
