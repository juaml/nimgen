import os
from pathlib import Path
import glob

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from ..expressions import (
    get_gene_expression,
    correlated_gene_expression,
    gene_coexpression
)
from ..stats import empirical_pval
from ..web import run_webgestalt
from ..smash import (
    generate_distance_matrices,
    export_voxel_coordinates,
    generate_surrogate_map
)
from ..utils import remove_nii_extensions


def _save_correlation_matrices(
    parcellation_file, significant_genes, output_path, metric="spearman"
):

    # calculate gene co-expression based on significant genes
    corr_gene_exp_matrix = correlated_gene_expression(
        parcellation_file, significant_genes, metric=metric
    )
    coexp_matrix = gene_coexpression(
        parcellation_file, significant_genes, metric=metric
    )

    # calculate gene co-expression based on all genes
    corr_all_gene_exp_matrix = correlated_gene_expression(
        parcellation_file, "all", metric=metric
    )
    coexp_all_matrix = gene_coexpression(
        parcellation_file, "all", metric=metric
    )

    # save as csv in the output path
    coexp_matrix.to_csv(
        os.path.join(
            output_path, f'significant-genes_gene_by_gene_'
            f'correlation_matrix_{metric}.tsv'
        ), sep="\t"
    )
    corr_gene_exp_matrix.to_csv(
        os.path.join(
            output_path, f'significant-genes_region_by_region_'
            f'correlation_matrix_{metric}.tsv'
        ), sep="\t"
    )

    # save matrices based on all genes
    coexp_all_matrix.to_csv(
        os.path.join(
            output_path,
            f'all-genes_gene_by_gene_correlation_matrix_{metric}.tsv'
        ), sep="\t"
    )
    corr_all_gene_exp_matrix.to_csv(
        os.path.join(
            output_path,
            f'all-genes_region_by_region_correlation_matrix_{metric}.tsv'
        ), sep="\t"
    )


def step_1(parcellation_file):
    """ Step 1 in HTCondor-based pipeline.
    Create a distance matrix for a given parcellation scheme so that surrogate
    parcellations can be created using BrainSmash. This step is quite time
    intensive and to speed up the following steps and their (potential)
    re-computation it is recommended to cache the results from this step on
    disk. However, depending on size and resolution of the parcellation nifti
    file this can take a lot of space (~100GB).

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        path to the parcellation file within the nimgen pipeline project
        directory


    Returns
    -------
    """
    if not os.path.isfile(parcellation_file):
        raise ValueError('Input file not found.')

    path_to_parc, _ = os.path.split(parcellation_file)
    # `generate_distance_matrix` already checks if the distance matrix exists,
    # so no further checks necessary here.
    # Return values are not needed here, brainsmash already saves these files
    # in specified path
    export_voxel_coordinates(parcellation_file, path_to_parc)
    generate_distance_matrices(path_to_parc)


def step_2(
    parcellation_file,
    marker_file,
    name_markers_dir,
    name_output_dir,
    smap_id,
    allen_data_dir,
    correlation_method="spearman",
    alpha=0.05,
    n_pca_covariates=None,
    partial_correlation=False,
):
    """ Step 2 in HTCondor-based pipeline.
    Create a single instance of a surrogate map using the distance matrix
    generated in step 1. Depending on the number of permutations desired for
    statistical testing and calculation of empirical p-values in the following
    mass-univariate correlation analysis, this step should be repeated n times,
    each time with a unique ID, i.e. for n jobs use str(i) for i in range(n).
    The surrogate map will then be used for a correlation analysis and results
    will be saved for later use.

    Parameters
    ----------

    Returns
    -------

    """
    if not os.path.isfile(parcellation_file):
        raise ValueError('Input file not found.')

    path_to_parc, name_parc_ext = os.path.split(parcellation_file)
    name_parc = remove_nii_extensions(name_parc_ext)

    path_to_marker_dir, marker_structure = Path(
        marker_file
    ).absolute().as_posix().split(name_markers_dir)
    head, _ = os.path.split(marker_structure)
    path_to_specific_marker_output = os.path.join(
        path_to_marker_dir, name_output_dir, name_parc, head
    )

    voxel_parcel_file = os.path.join(path_to_parc, "brain_map.txt")
    matrix_files = {
        'D': os.path.join(path_to_parc, 'distmat.npy'),
        'index': os.path.join(path_to_parc, 'index.npy')
    }

    smap_id_str = (
        f"smapid-{smap_id}-correlationmethod-{correlation_method}"
        f"_alpalevel-{alpha}_npcacovariates-{n_pca_covariates}"
        f"_partialcorrelation-{partial_correlation}"
    )

    # generate surrogate map for given atlas
    surrogate_map = generate_surrogate_map(
        parcellation_file,
        smap_id_str,
        path_to_parc,
        voxel_parcel_file,
        matrix_files
    )

    if n_pca_covariates is None:
        perform_pca = False
        pca_dict = None
    else:
        assert isinstance(n_pca_covariates, int), (
            "n_pca_covariates has to be an integer!"
        )
        pca_dict = {"n_components": n_pca_covariates}

    # perform correlation analysis for given surrogate map and marker
    all_genes_corr_scores, _, _ = get_gene_expression(
        marker=marker_file,
        atlas=surrogate_map,
        aggregation_method="mean",
        allen_data_dir=allen_data_dir,
        correlation_method=correlation_method,
        alpha=alpha,
        perform_pca=perform_pca,
        pca_dict=pca_dict,
        partial_correlation=partial_correlation,
    )

    # save correlation results for all genes for this surrogate map and these
    # settings
    all_genes_corr_scores.to_csv(
        os.path.join(
            path_to_specific_marker_output,
            'smap_corr_scores', f'{smap_id_str}.tsv'
        ), sep="\t"
    )


def step_3(
    parcellation_file,
    marker_file,
    name_markers_dir,
    name_output_dir,
    allen_data_dir,
    r_path,
    correlation_method="spearman",
    alpha=0.05,
    n_pca_covariates=None,
    partial_correlation=False,
):
    """ Step 3 in HTCondor-based pipeline.
    Perform mass-univariate correlation analysis, calculate empirical p-values
    using surrogate results, and export significant genes, as well as gene
    co-expression matrix, region-wise gene expression correlation matrix.

    """
    smap_corr_files = glob.glob(
        f"smapid-*-correlationmethod-{correlation_method}"
        f"_alpalevel-{alpha}_npcacovariates-{n_pca_covariates}_"
        f"_partialcorrelation-{partial_correlation}"
    )

    # read, concat, delete p-val column from the smashed correlation df
    smashed_data = []
    for f in smap_corr_files:
        smashed_data.append(pd.read_csv(f, sep="\t", index_col=0))
    smashed_corr_df = pd.concat(smashed_data, axis=1)

    # prepare output directories
    if not os.path.isfile(parcellation_file):
        raise ValueError('Input file not found.')

    path_to_parc, name_parc_ext = os.path.split(parcellation_file)
    name_parc = remove_nii_extensions(name_parc_ext)

    path_to_marker_dir, marker_structure = Path(
        marker_file
    ).absolute().as_posix().split(name_markers_dir)
    head, _ = os.path.split(marker_structure)
    path_to_specific_marker_output = os.path.join(
        path_to_marker_dir, name_output_dir, name_parc, head
    )

    # prepare potential pca
    if n_pca_covariates is None:
        perform_pca = False
        pca_dict = None
    else:
        assert isinstance(n_pca_covariates, int), (
            "n_pca_covariates has to be an integer!"
        )
        pca_dict = {"n_components": n_pca_covariates}

    # perform correlation analysis for given surrogate map and marker
    all_genes_corr_scores, _, covariates_dict_of_niftis = get_gene_expression(
        marker=marker_file,
        atlas=parcellation_file,
        aggregation_method="mean",
        allen_data_dir=allen_data_dir,
        correlation_method=correlation_method,
        alpha=alpha,
        perform_pca=perform_pca,
        pca_dict=pca_dict,
        partial_correlation=partial_correlation,
    )
    real_correlations = all_genes_corr_scores["r_score"].T.values
    smashed_correlations = smashed_corr_df["r_score"].T.values

    empirical_pvalues = empirical_pval(smashed_correlations, real_correlations)

    all_genes_corr_scores["empirical_pvals"] = empirical_pvalues
    reject, corrected, *_ = multipletests(
        all_genes_corr_scores["empirical_pvals"],
        alpha=alpha,
        method="fdr_bh",
        is_sorted=False,
        returnsorted=False
    )
    all_genes_corr_scores["fdr_bh_corrected_empirical_pvals"] = corrected
    all_genes_corr_scores[f"reject_at_alpha-{alpha}"] = reject

    significant_genes_df = all_genes_corr_scores[reject]
    significant_genes = significant_genes_df["genes"]

    if partial_correlation and perform_pca:
        output_path = os.path.join(
            path_to_specific_marker_output, "pca_covariates",
            f"{n_pca_covariates}_component_pca"
        )
    else:
        output_path = path_to_specific_marker_output

    genes_file = os.path.join(
        output_path,
        'significant-empirical-pvalue-fdr-corrected_genes.txt'
    )
    np.savetxt(
        genes_file, significant_genes, fmt="%s"
    )

    # Run gene set enrichment analysis via webgestalt R package
    run_webgestalt(r_path=r_path, genes_file=genes_file)

    if isinstance(covariates_dict_of_niftis, dict):
        for key, value in covariates_dict_of_niftis.items():
            value.to_filename(
                os.path.join(output_path, f"{key}.nii.gz")
            )

    for metric in ["spearman", "pearson"]:
        _save_correlation_matrices(
            parcellation_file, significant_genes, output_path, metric=metric
        )
