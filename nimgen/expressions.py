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
import nibabel as nib
import abagen
from .utils import logger
import pingouin as pg
from sklearn.decomposition import PCA
from nilearn import image


def _save_expressions(exp, atlas):
    logger.info(f"Trying to save expressions")
    save_path = atlas.parent
    atlas_name = atlas.stem
    if os.access(save_path, os.W_OK):
        exp_fname = save_path / f"nimgen_{atlas_name}_expressions.csv"
        logger.info(
            f"Saving expressions to nimgen_{atlas_name}_expressions.csv")
        exp.to_csv(exp_fname, sep=";", index=False)
    else:
        logger.info(
            "User does not have permissions to save "
            f"to {save_path.as_posix()}")


def _get_cached_results(atlas):
    if not isinstance(atlas, Path):
        atlas = Path(atlas)
    exp_path = atlas.parent
    atlas_name = atlas.stem
    exp_fname = exp_path / f"nimgen_{atlas_name}_expressions.csv"
    exp = None
    if exp_fname.exists():
        logger.info(f"Reading expressions from {exp_fname.as_posix()}")

        exp = pd.read_csv(exp_fname, sep=";")
    return exp


def apply_pca(exp, pca_dict=None):
    """
    Apply Principal component analysis to given expressions.

    Parameters
    ----------
    exp : np.ndarray or pandas.Series or pandas.DataFrame
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

    pca = PCA(**pca_dict)
    covariates = pca.fit_transform(exp.values.copy())
    covariates_df = pd.DataFrame(covariates)
    covariates_df.columns = [
        f"covariate-{(int(i)+1)}" for i, _ in enumerate(covariates_df)
    ]
    return covariates_df


def covariates_to_nifti(parcellation, covariates_df):
    """
    Creates nifti images for given PCA covariates.

    Parameters
    ----------
    parccellation : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
    covariates_df : dict
        PCA covariates. Each key represents different covariate.

    Returns
    -------
    covariate_niftis : dict
        A dictionary contains niimg-like object for each covariate.
    """
    parcellation_array = np.array(parcellation.dataobj)
    covariate_niftis = {}
    for covariate_label, covariate in covariates_df.items():
        null_mat = np.zeros(parcellation_array.shape)
        null_mat[parcellation_array == 0] = -20000
        marker_min = np.min(covariate)
        marker_max = np.max(covariate)
        for label, value in covariate.iteritems():
            null_mat[parcellation_array == label + 1] = value

        pc_nii = image.new_img_like(
            parcellation, null_mat
        )
        pc_nii.header["cal_min"] = marker_min
        pc_nii.header["cal_max"] = marker_max
        covariate_niftis[covariate_label] = pc_nii

    return covariate_niftis


def _read_sign_genes(sign_genes):
    if isinstance(sign_genes, pd.DataFrame):
        sign_genes = sign_genes.index
    elif not isinstance(sign_genes, pd.DataFrame):
        _, ext = os.path.splitext(sign_genes)
        if ext in [".csv", ".tsv"]:
            extensions = {".csv": ",", ".tsv": "\t"}
            sign_genes = pd.read_csv(
                sign_genes,
                header=None,
                index_col=0,
                dtype=str,
                sep=extensions[ext]
            ).index
        elif ext in [".txt"]:
            sign_genes = list(np.loadtxt(sign_genes, dtype=str))
        else:
            raise ValueError(
                "'sign_genes' should be a pd.DataFrame,"
                " .csv/.tsv or .txt file!"
            )

    return list(sign_genes)


def correlated_gene_expression(parcellation, sign_genes, metric="spearman"):
    """
    Performs correlated gene expression analysis for original parcellation
    file using significant genes. (ROIxROI matrix)

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
    if isinstance(sign_genes, str):
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

    local_exp = sign_genes_expression.T.copy()

    return local_exp.corr(method=metric)


def gene_coexpression(parcellation, sign_genes, metric="spearman"):
    """
    Performs gene co-expression analysis for original parcellation file using
    significant genes. (Gene x Gene matrix)

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
    if isinstance(sign_genes, str):
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
        exp,
        markers,
        correlation_method,
        partial_correlation,
        covariates_df):
    """
    Applies correlation analysis for given gene expressions and markers.
    If partial correlation is False, performs normal correlation based on
    correlation_method. If partial correlation is True, covariates_df should
    be given.

    Parameters
    ----------
    exp : np.ndarray or pandas.Series or pandas.DataFrame
        Gene expression values ROIxGenes
    markers : list(float) or np.ndarray or pandas.Series or pandas.DataFrame
        Markers for each ROI in the atlas
    correlation_method : string
        Correlation method, pearson or spearman.
    partial_correlation : bool
        Performs partial correlation with given covariates and markers using
        pg.partial_corr function. If set True, covariates_df should be given.
    covariates_df : pandas.Dataframe
        A dataframe object contains covariates after PCA analysis.
        Works only with partial_correlation=True.

    Returns
    -------
    pval, r_score : dict
        A dictionary contains p-value and correlation score.
    """

    if partial_correlation:
        correlation_results_list = []
        for gene, gene_vector in exp.iteritems():
            gene_spec_data = covariates_df.copy()
            gene_spec_data.index = exp.index
            gene_spec_data[gene] = gene_vector
            gene_spec_data["marker"] = markers.squeeze().values

            correlation_results_list.append(
                pg.partial_corr(
                    gene_spec_data,
                    x=gene, y="marker",
                    x_covar=covariates_df.columns,
                    method=correlation_method
                )
            )
        correlation_result = pd.concat(correlation_results_list)
        correlation_result.index = exp.columns

        pval, r_score = correlation_result["p-val"], correlation_result["r"]
    else:
        corr_funcs = {
            'spearman': stats.spearmanr,
            'pearson': stats.pearsonr
        }
        corr_func = corr_funcs[correlation_method]

        correlation_result = exp.apply(
            lambda col: corr_func(col.values, markers.squeeze().values),
            result_type="expand",
        ).T

        pval, r_score = correlation_result[1], correlation_result[0]

    return pval, r_score


def get_gene_expression(
        markers,
        atlas,
        allen_data_dir=None,
        save_expressions=True,
        force_recompute=False,
        correlation_method='spearman',
        alpha=0.05,
        perform_pca=False,
        pca_dict=None,
        partial_correlation=False,
        custom_covariates_df=None):
    """Get the genes expressed in the atlas that correlate with the
    specified markers.

    Parameters
    ----------
    markers : list(float) or np.ndarray or pandas.Series or pandas.DataFrame
        Markers for each ROI in the atlas
    atlas : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
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

    abagen_params = {"probe_selection": "diff_stability"}
    if not isinstance(markers, pd.DataFrame):
        markers = pd.DataFrame(markers)

    if not isinstance(atlas, Path):
        atlas = Path(atlas)

    if not atlas.exists():
        raise ValueError(f"Atlas file does not exist: {atlas.as_posix()}")

    logger.info("Checking atlas and markers dimensions")
    atlas_img = nib.load(atlas)
    nrois = np.unique(atlas_img.get_fdata()).astype(np.int).shape[0] - 1
    if nrois != len(markers):
        raise ValueError(
            f"Number of markers ({len(markers)}) does not match "
            f"the number of ROIs in the atlas ({nrois})."
        )

    exp = None
    if force_recompute is False:
        # Check if we have this results
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
    exp_with_nan = exp.copy()
    exp = exp[good_rois]
    markers = markers[good_rois]

    covariates_df = None
    custom_covariates_df = None
    covariates_nifti = None

    if perform_pca:
        logger.info("Principal component analysis (PCA) started..")
        covariates_df = apply_pca(exp, pca_dict)
        covariates_df_with_nans = pd.DataFrame(
            np.zeros((exp_with_nan.shape[0], covariates_df.shape[1]))
        )
        covariates_df_with_nans.columns = covariates_df.columns
        covariates_df_with_nans.iloc[good_rois] = covariates_df
        covariates_df_with_nans.iloc[bad_rois] = np.nan
        covariates_nifti = covariates_to_nifti(
            atlas_img, covariates_df_with_nans)
    else:
        # if perform_pca False but partial_correlation True
        # then use custom_covarites_df
        if partial_correlation:
            covariates_df = custom_covariates_df

    pval, r_score = correlation_analysis(exp, markers, correlation_method,
                                         partial_correlation, covariates_df
                                         )

    all_genes = pd.DataFrame(
        {"genes": pval.index, "pval": pval.values, "r_score": r_score}
    ).set_index("genes")

    sign_genes = all_genes[all_genes.pval < 0.05]

    return all_genes, sign_genes, covariates_nifti
