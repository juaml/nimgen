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
                        save_expressions=True, force_recompute=False,
                        multiple_correction='fdr_bh', alpha=0.05):
    """Get the genes expressed in the atlas that correlate with the
    specified weights.

    Parameters
    ----------
    weights : list(float) or np.ndarray or pandas.Series or pandas.DataFrame
        Weights for each ROI in the atlas
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
    allen_data_dir : pathlib.Path or string
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded.
    save_expressions : bool
        If True (default), the expressions on the atlas will be saved for later
        calls to the function with the same atlas. It will create a file named
        `nimgen_{atlas}_expressions.csv` next to the atlas file.
    force_recompute : bool
        If True, disregard the previously saved expressions and recompute the
        gene expressions using abagen. Defaults to False.
    multiple_correction : str
        Method to use for the correction for multiple comparisons. Check
        `statsmodels.stats.multitest.mutipletests` for a list of available
        options. Defaults to 'fdr_bh'
    alpha : float
        FWER, family-wise error rate, e.g. 0.1. Defaults to 0.05.

    Returns
    -------
    all_genes : object(DataFrame)
        An object will all the genes and p-values that are expressed in the 
        atlas.
    sign_genes : object(DataFrame)
        An object will significant genes and p-values that are expressed in the 
        atlas and correlate with the weights.

    Raises
    ------
    ValueError
        If there is a problem with the Atlas and/or the weights.

    """
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
    pearson_result = exp.apply(lambda col: stats.pearsonr (col.values, 
                            weights.squeeze().values), result_type='expand').T
    pval, r_score = pearson_result[1], pearson_result[0]

    all_genes = pd.DataFrame({ 'genes': pval.index, 'pval': pval.values, 'r_score': r_score}).set_index('genes')
    sign_genes = None
    if multiple_correction is not None:
        reject, pvals_corrected, *_ = multipletests(pval, alpha=alpha, method=multiple_correction)
        all_genes['pvals_corrected'] = pvals_corrected
        all_genes = all_genes.sort_values(by=["r_score"])    
        genes = pval[reject].index.values.tolist()
        sign_genes = all_genes[all_genes.index.isin(genes)]
        sign_genes.index.name = 'sign_genes'
    else:
        sign_genes = all_genes[ all_genes.pval < 0.05 ]
        
    return all_genes, sign_genes
