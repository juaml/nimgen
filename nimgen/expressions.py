"""Fetch gene expression data from AHBA and run mass-univariate analysis."""

# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Leonard Sasse <l.sasse@fz-juelich.de>
#          Yasir Demirtaş <tyasird@gmail.com>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
#          Kaustubh Patil <k.patil@fz-juelich.de>
# License: AGPL

import os
import warnings
from pathlib import Path

import abagen
import nibabel as nib
import numpy as np
import pandas as pd
import pingouin as pg
from nilearn import image, masking
from scipy import stats
from sklearn.decomposition import PCA

from .statistics import _get_funcbyname
from .utils import (
    _read_sign_genes,
    covariates_to_nifti,
    logger,
    remove_nii_extensions,
)


def _save_expressions(exp, atlas):
    logger.info("Trying to save expressions")
    save_path, atlas_name = os.path.split(remove_nii_extensions(atlas))
    if os.access(save_path, os.W_OK):
        exp_fname = os.path.join(
            save_path, f"nimgen_{atlas_name}_expressions.csv"
        )
        logger.info(
            f"Saving expressions to nimgen_{atlas_name}_expressions.csv"
        )
        exp.to_csv(exp_fname, sep=";", index=False)
    else:
        logger.info(
            "User does not have permissions to save "
            f"to {save_path.as_posix()}"
        )


def _get_cached_results(atlas):
    exp_path, atlas_name = os.path.split(remove_nii_extensions(atlas))
    exp_fname = os.path.join(exp_path, f"nimgen_{atlas_name}_expressions.csv")
    exp = None
    if os.path.isfile(exp_fname):
        logger.info(f"Reading expressions from {exp_fname}")
        exp = pd.read_csv(exp_fname, sep=";")
    else:
        logger.info("No cached results...")
    return exp


def apply_pca(exp, pca_dict=None):
    """
    Apply Principal component analysis to given expressions.

    Parameters
    ----------
    exp : np.ndarray or pandas.DataFrame
        Gene expression values ROIxGenes
    pca_dict : dict
        PCA parameters i.e. n_components

    Returns
    -------
    covariates_df : object(DataFrame)
        An object will all the covariates after PCA analysis.
    """
    if pca_dict is None:
        pca_dict = {"n_components": 5}

    if isinstance(exp, pd.DataFrame):
        exp = exp.values

    pca = PCA(**pca_dict)
    covariates = pca.fit_transform(exp.copy())
    covariates_df = pd.DataFrame(covariates)
    covariates_df.columns = [
        f"covariate-{(int(i)+1)}" for i, _ in enumerate(covariates_df)
    ]
    return covariates_df


def correlated_gene_expression(parcellation, sign_genes, metric="spearman"):
    """Perform correlated gene expression analysis (ROIxROI matrix).

    Parameters
    ----------
    parcellation : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
    sign_genes : pandas.DataFrame, or path to .csv, .tsv, or text file.
        if it's a pd.DataFrame, the index should contain names of genes,
        if its a csv or tsv file, the first column should contain names of
        genes,
        if its a txt file, it should simply be a column of the gene names
        Significant gene list.

        alternatively this can be the string "all" meaning all genes will be
        included

    Returns
    -------
    pandas.DataFrame
        A dataframe (ROIxROI matrix) contains correlation score for each ROI.
    """
    exp_with_nans = _get_cached_results(parcellation)
    if isinstance(sign_genes, str) and not (os.path.isfile(sign_genes)):
        if sign_genes in ["all"]:
            sign_genes_expression = exp_with_nans
        else:
            raise ValueError(
                "if 'sign_genes' is a str it should be 'all' or path to file!"
            )
    else:
        sign_genes = _read_sign_genes(sign_genes)

        for g in sign_genes:
            assert g in exp_with_nans.columns, (
                f"{g} not in the gene expression data set!"
                "Make sure significant genes are provided in "
                "the correct format!"
            )

        sign_genes_expression = exp_with_nans[sign_genes]

    local_exp = sign_genes_expression.T.copy()

    return local_exp.corr(method=metric)


def gene_coexpression(parcellation, sign_genes, metric="spearman"):
    """Perform gene co-expression analysis for (Gene x Gene matrix).

    Parameters
    ----------
    parccellation : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
    sign_genes : pandas.DataFrame
        Significant gene list.

        if it's a pd.DataFrame, the index should contain names of genes,
        if its a csv or tsv file, the first column should contain names of
        genes,
        if its a txt file, it should simply be a column of the gene names
        Significant gene list.

        alternatively this can be the string "all" meaning all genes will be
        included


    Returns
    -------
    pandas.DataFrame
        A dataframe (Gene x Gene matrix) contains correlation score for each
        gene.
    """
    exp_with_nans = _get_cached_results(parcellation)
    if isinstance(sign_genes, str) and not os.path.isfile(sign_genes):
        if sign_genes in ["all"]:
            sign_genes_expression = exp_with_nans
        else:
            raise ValueError("if 'sign_genes' is a str it should be 'all'!")
    else:
        sign_genes = _read_sign_genes(sign_genes)

        for g in sign_genes:
            assert g in exp_with_nans.columns, (
                f"{g} not in the gene expression data set!"
                "Make sure significant genes are provided in "
                "the correct format!"
            )

        sign_genes_expression = exp_with_nans[sign_genes]

    local_exp = sign_genes_expression.copy()
    return local_exp.corr(method=metric)


def correlation_analysis(
    exp, markers, correlation_method, partial_correlation, covariates_df
):
    """Apply correlation analysis for given gene expressions and markers.

    If partial correlation is False, performs normal correlation based on
    correlation_method. If partial correlation is True, covariates_df should
    be given.

    Parameters
    ----------
    exp : pandas.DataFrame
        Gene expression values ROIxGenes
    markers : list(float) or np.ndarray or pandas.Series or pandas.DataFrame
        Markers for each ROI in the atlas
    correlation_method : str
        Correlation method, 'pearson' or 'spearman'.
    partial_correlation : bool
        Performs partial correlation with given covariates and markers using
        pg.partial_corr function. If set True, covariates_df should be given.
    covariates_df : pandas.Dataframe
        A dataframe object contains covariates after PCA analysis.
        Works only with partial_correlation=True.

    Returns
    -------
    pval, r_score : tuple
        [0] p-values and [1] correlation scores.
    """
    if partial_correlation:
        covariates = covariates_df.columns
        correlation_results_list = []
        for gene, gene_vector in exp.iteritems():
            gene_spec_data = covariates_df.copy()
            gene_spec_data.index = exp.index
            gene_spec_data[gene] = gene_vector
            gene_spec_data["marker"] = np.array(markers.squeeze())

            correlation_results_list.append(
                pg.partial_corr(
                    gene_spec_data,
                    x=gene,
                    y="marker",
                    x_covar=covariates,
                    method=correlation_method,
                )
            )
        correlation_result = pd.concat(correlation_results_list)
        correlation_result.index = exp.columns

        pval, r_score = correlation_result["p-val"], correlation_result["r"]
    else:
        corr_funcs = {"spearman": stats.spearmanr, "pearson": stats.pearsonr}
        corr_func = corr_funcs[correlation_method]

        correlation_result = exp.apply(
            lambda col: corr_func(col.values, np.array(markers.squeeze())),
            result_type="expand",
        ).T

        pval, r_score = correlation_result[1], correlation_result[0]

    return pval, r_score


def get_gene_expression(
    marker,
    atlas,
    aggregation_method="mean",
    allen_data_dir=None,
    save_expressions=True,
    force_recompute=False,
    correlation_method="spearman",
    alpha=0.05,
    perform_pca=False,
    pca_dict=None,
    partial_correlation=False,
    custom_covariates_df=None,
):
    """Get the genes expressed that correlate with the specified markers.

    Parameters
    ----------
    marker : str or os.PathLike or niimg
        Can be a path to a parcellated marker nifti-file or the nifti-object
        already loaded
    atlas : str or os.PathLike
        path to a parcellation image in MNI space, where each parcel is
        identified by a unique integer ID.
    aggregation_method : str
        method to aggregate the marker given the parcellation. Can be
        'winsorized_mean', 'mean', or 'std'. Default is 'mean'.
    allen_data_dir : pathlib.Path or string
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded.
    save_expressions : bool, default = True
        If True, the expressions on the atlas will be saved for later
        calls to the function with the same atlas. It will create a file named
        `nimgen_{atlas}_expressions.csv` next to the atlas file.
    force_recompute : bool
        If True, disregard the previously saved expressions and recompute the
        gene expressions using abagen. Defaults to False.
    correlation_method : str, default = 'spearman'
        Method to use for the correlation.
    alpha : float
        FWER, family-wise error rate, e.g. 0.1. Defaults to 0.05.
    perform_pca : bool, default = False
        Performs Principal component analysis to gene expression values.
    pca_dict : dict, default = None
        Dictinoary for PCA variables i.e. {"n_components": 5}
    partial_correlation : bool, default = False
        Performs partial correlation to gene expression, marker and covarite.
    custom_covariates_df : dict, default = None
        If perform_pca is False and partial_correlation is True,
        custom_covariates_df should be defined.

    Returns
    -------
    all_genes : object(DataFrame)
        An object will all the genes and p-values that are expressed in the
        atlas.
    sign_genes : object(DataFrame)
        An object will significant genes and p-values that are expressed in the
        atlas and correlate with the markers.
    pca_components : dict
        A dictionary contains niimg-like object for each component.

    Raises
    ------
    ValueError
        If there is a problem with the Atlas and/or the markers.

    """

    marker_aggregated, _ = _aggregate_marker(atlas, marker)
    marker = pd.DataFrame(marker_aggregated[aggregation_method])

    # parcellate gene expression data and extract ROI's where
    # gene expression levels are `NaN`
    expressions, good_rois, bad_rois = _prepare_expressions(
        marker,
        atlas,
        force_recompute=force_recompute,
        allen_data_dir=allen_data_dir,
        save_expressions=save_expressions,
    )
    exp_no_nan, marker_no_nan = expressions[good_rois], marker[good_rois]

    # prepare covariates if partial correlation is desired
    covariates_dict_of_niftis = None
    if partial_correlation:
        assert perform_pca or custom_covariates_df, (
            "If partial_correlation is True you should provide some covariates"
            " or opt to perform a pca!"
        )
        covariates_df, covariates_dict_of_niftis = _prepare_covariates(
            expressions,
            atlas,
            good_rois,
            bad_rois,
            perform_pca,
            pca_dict,
            custom_covariates_df,
        )
    elif perform_pca and (custom_covariates_df is not None):
        warnings.warn(
            "partial_correlation is set to False, but either perform_pca is "
            "set to True or custom_covariates_df is not None!"
            " Partial correlation will not be performed, "
            "but pca will be performed and pc's will bereturned as niftis."
        )
        covariates_df, covariates_dict_of_niftis = _prepare_covariates(
            expressions,
            atlas,
            good_rois,
            bad_rois,
            perform_pca,
            pca_dict,
            custom_covariates_df,
        )
    elif perform_pca:
        warnings.warn(
            "perform_pca is set to True, but partial_correlation"
            " is set to False. Partial correlation will not be performed, "
            "but pca will be performed and pc's will bereturned as niftis."
        )
        covariates_df, covariates_dict_of_niftis = _prepare_covariates(
            expressions,
            atlas,
            good_rois,
            bad_rois,
            perform_pca,
            pca_dict,
            custom_covariates_df,
        )
    else:
        covariates_df = None

    # perform mass-univariate correlation analysis
    pval, r_score = correlation_analysis(
        exp_no_nan,
        marker_no_nan,
        correlation_method,
        partial_correlation,
        covariates_df,
    )

    all_genes = pd.DataFrame(
        {"genes": pval.index, "pval": pval.values, "r_score": r_score}
    ).set_index("genes")

    sign_genes = all_genes[all_genes.pval < alpha]

    return all_genes, sign_genes, covariates_dict_of_niftis


def _aggregate_marker(atlas, vbm, aggregation=None, limits=None):
    """Construct a masker based on the input atlas_nifti.

    Applies resampling of the atlas if necessary and applies the masker to
    the vbm_nifti to extract brain-imaging based vbm markers.
    So far the aggregation methods "winsorized mean", "mean" and
    "std" are supported.

    Parameters
    ----------
    atlas_nifti : niimg-like object
        Nifti of atlas to use for parcellation.
    vbm_nifti: niimg-like object
        Nifti of voxel based morphometry as e.g. outputted by CAT.
    aggregation: list
        List with strings of aggregation methods to apply. Defaults to
        aggregation = ['winsorized_mean', 'mean', 'std'].
    limits: array
        Array with lower and upper limit for the calculation of the winsorized
        mean. Only needed when 'winsorized_mean' was specified
        in aggregation. If wasn't specified defaults to [0.1, 0.1].

    Returns
    -------
    marker_aggregated : dict
        Dictionary with keys being each of the chosen aggregation methods
        and values the corresponding array with the calculated marker based on
        the provided atlas. The array therefore as the shape of the chosen
        number of ROIs (granularity).
    marker_func_params: dict
        Dictionary with parameters used for the aggregation function. Keys:
        respective aggregation function, values: dict with responding
        parameters
    """

    atlas_nifti = image.load_img(atlas)
    vbm_nifti = image.load_img(vbm)

    # defaults (validity is checked in _get_funcbyname())
    if aggregation is None:  # Don't put mutables as defaults, use None instead
        aggregation = ["winsorized_mean", "mean", "std", "median"]
    if limits is None:
        limits = [0.1, 0.1]

    # aggregation function parameters (validity is checked in
    # _get_funcbyname())
    agg_func_params = {"winsorized_mean": {"limits": limits}}

    # definitions
    # sort rois to be related to the order of i_roi (and get rid of 0 entry)
    rois = sorted(np.unique(image.get_data(atlas_nifti)))[1:]  # roi numbering
    n_rois = len(rois)  # granularity
    marker_aggregated = {
        x: np.ones(shape=(n_rois)) * np.nan for x in aggregation
    }

    # resample atlas if needed
    if not atlas_nifti.shape == vbm_nifti.shape:
        atlas_nifti = image.resample_to_img(
            atlas_nifti, vbm_nifti, interpolation="nearest"
        )
        logger.info("Atlas nifti was resampled to resolution of VBM nifti.")

    logger.info("make masker and apply")

    # make masker and apply
    for i_roi, roi in enumerate(rois):
        mask = image.math_img(f"img=={roi}", img=atlas_nifti)
        marker = masking.apply_mask(
            imgs=vbm_nifti, mask_img=mask
        )  # gmd per roi
        # logger.info(f'Mask applied for roi {roi}.')
        # aggregate (for all aggregation options in list)
        for agg_name in aggregation:
            # logger.info(f'Aggregate GMD in roi {roi} using {agg_name}.')
            agg_func = _get_funcbyname(
                agg_name, agg_func_params.get(agg_name, None)
            )
            marker_aggregated[agg_name][i_roi] = agg_func(marker)
    logger.info(f"{aggregation} was computed for all {n_rois} ROIs.\n")

    return marker_aggregated, agg_func_params


def _prepare_expressions(
    marker,
    atlas,
    force_recompute=False,
    allen_data_dir=None,
    save_expressions=True,
):
    """Prepare parcellated gene expressions and marker.

    Parameters
    ----------
    allen_data_dir : pathlib.Path or string
        Directory where expression data should be downloaded (if it does not
        already exist) / loaded.
    save_expressions : bool, default = True
        If True, the expressions on the atlas will be saved for later
        calls to the function with the same atlas. It will create a file named
        `nimgen_{atlas}_expressions.csv` next to the atlas file.
    force_recompute : bool
        If True, disregard the previously saved expressions and recompute the
        gene expressions using abagen. Defaults to False.

    Returns
    -------
    exp : pd.DataFrame
        all parcellated gene expressions
    marker : pd.DataFrame
        parcellated marker
    """

    abagen_params = {"probe_selection": "diff_stability"}

    if not isinstance(atlas, Path):
        atlas = Path(atlas)

    if not atlas.exists():
        raise ValueError(f"Atlas file does not exist: {atlas.as_posix()}")

    logger.info("Checking atlas and markers dimensions")
    atlas_img = nib.load(atlas)
    nrois = np.unique(atlas_img.get_fdata()).astype(np.int).shape[0] - 1
    if nrois != len(marker):
        raise ValueError(
            f"Number of markers ({len(marker)}) does not match "
            f"the number of ROIs in the atlas ({nrois})."
        )

    exp = None
    if force_recompute is False:
        # Check if we have these results
        exp = _get_cached_results(atlas)

    if exp is None or force_recompute is True:
        # If not, compute the expressions
        logger.info("Computing expressions, please wait")
        exp = abagen.get_expression_data(
            atlas, data_dir=allen_data_dir, **abagen_params
        )
        if save_expressions is True:
            _save_expressions(exp, atlas)

    # use only ROIs without NaNs (ROIs with samples)
    good_rois = ~exp.iloc[:, 0].isna().values
    bad_rois = exp.iloc[:, 0].isna().values

    return exp, good_rois, bad_rois


def _prepare_covariates(
    expression_data,
    atlas_img,
    good_rois,
    bad_rois,
    perform_pca,
    pca_dict=None,
    custom_covariates_df=None,
):
    if perform_pca:
        logger.info("Principal component analysis (PCA) started...")
        covariates_df = apply_pca(expression_data[good_rois], pca_dict)
        covariates_df_with_nans = pd.DataFrame(
            np.zeros((expression_data.shape[0], covariates_df.shape[1]))
        )
        covariates_df_with_nans.columns = covariates_df.columns
        covariates_df_with_nans.iloc[good_rois] = covariates_df
        covariates_df_with_nans.iloc[bad_rois] = np.nan
        covariates_nifti = covariates_to_nifti(
            atlas_img, covariates_df_with_nans
        )
        if custom_covariates_df is not None:
            raise NotImplementedError(
                "Use of both custom covariates and pca covariates not"
                "implemented yet!"
            )
    else:

        covariates_df = custom_covariates_df

    return covariates_df, covariates_nifti
