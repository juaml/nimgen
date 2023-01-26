"""Basic steps to run the mass-univariate nimgen correlation analysis."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from neuromaps.parcellate import Parcellater
from statsmodels.stats.multitest import multipletests

from ..expressions import (
    correlated_gene_expression,
    gene_coexpression,
    gene_expression_correlations,
)
from ..smash import cached_distance_matrix, cached_null_maps
from ..statistics import empirical_pval
from ..utils import remove_nii_extensions
from ..web import run_webgestalt


logger = logging.getLogger(__name__)


def _save_correlation_matrices(
    parcellation_file,
    significant_genes,
    output_path,
    specific_marker_output,
    metric="spearman",
):

    path_all_genes_based_output = (
        specific_marker_output.parent.parent.parent / "all_gene_outputs"
    )
    logger.info(
        "Calculate gene co-expression and region-by-region gene"
        " expression correlation profiles based on significant genes."
    )
    corr_gene_exp_matrix = correlated_gene_expression(
        parcellation_file, significant_genes, metric=metric
    )
    coexp_matrix = gene_coexpression(
        parcellation_file, significant_genes, metric=metric
    )

    # save as csv in the output path
    coexp_matrix.to_csv(
        os.path.join(
            output_path,
            f"significant-genes_gene_by_gene_"
            f"correlation_matrix_{metric}.tsv",
        ),
        sep="\t",
    )
    corr_gene_exp_matrix.to_csv(
        os.path.join(
            output_path,
            f"significant-genes_region_by_region_"
            f"correlation_matrix_{metric}.tsv",
        ),
        sep="\t",
    )

    # save matrices based on all genes
    all_gene_coexp = (
        path_all_genes_based_output
        / f"all-genes_gene_by_gene_correlation_matrix_{metric}.tsv"
    )
    if not all_gene_coexp.is_file():
        logger.info("Calculate gene co-expression for all genes!")
        coexp_all_matrix = gene_coexpression(
            parcellation_file, "all", metric=metric
        )
        coexp_all_matrix.to_csv(all_gene_coexp, sep="\t")

    all_genes_roixroi = (
        path_all_genes_based_output
        / f"all-genes_region_by_region_correlation_matrix_{metric}.tsv"
    )
    if not all_genes_roixroi.is_file():
        logger.info(
            "Calculate region-by-region gene expression "
            "correlation profiles for all genes!"
        )
        corr_all_gene_exp_matrix = correlated_gene_expression(
            parcellation_file, "all", metric=metric
        )
        corr_all_gene_exp_matrix.to_csv(all_genes_roixroi, sep="\t")


def _step_1(parcellation_file):
    """Run step 1 in nimgen pipeline.

    Create a distance matrix for a given parcellation scheme so that surrogate
    parcellations can be created using Neuromaps/BrainSmash.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        path to the parcellation file within the nimgen pipeline project
        directory

    Returns
    -------
    None; saves distance matrices in the appropriate parcellation directory
    """
    if not os.path.isfile(parcellation_file):  # TODO: change to pathlb
        raise ValueError("Input file not found.")

    logger.info("Checking for cached distance matrix...")
    cached_distance_matrix(parcellation_file)
    logger.info("Done!")


def _step_2(parcellation_file, marker_file, n_perm, seed):
    """Run step 2 in nimgen pipeline.

    Generate null maps for a combination of parcellation and marker.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        path to the parcellation file within the nimgen pipeline project
        directory.
    marker_file : str or os.PathLike
        path to the marker file within the nimgen pipeline project
        directory.
    n_perm : int
        number of null maps to generate.
    seed : int
        random seed for null map generation.
    """

    logger.info("Checking for cached distance null maps...")
    dist_mat = cached_distance_matrix(parcellation_file)
    cached_null_maps(
        parcellation_file, marker_file, dist_mat, n_perm, seed=seed
    )
    logger.info("Done!")


def _step_3(
    parcellation_file,
    marker_file,
    n_perm,
    smap_id,
    seed,
    allen_data_dir,
    correlation_method,
    n_pca_covariates=None,
):
    """Run step 3 in nimgen pipeline.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Path to the parcellation nifti file
    marker_file : str or os.PathLike
        Path to the marker nifti file
    n_perm : int
        number of null maps that were created in the first place.
    smap_id : int
        unique number identifying the surrogate map
    seed : int
        seed for operations involving randomness
    allen_data_dir : str or os.PathLike
        root directory of AHBA data
    correlation_method : str
        'spearman' or 'pearson'
    n_pca_covariates : int or None
        number of components gene expression components (after pca) to include
        as covariates in the partial correlation between marker and genes.

    Returns
    -------
    None; saves correlation scores for surrogate maps in appropriate output
    directory
    """
    logger.info("Starting surrogate correlation analysis...")
    for key, value in locals().items():
        logger.info(f"{key}     ==================      {value}")

    partial_correlation = False if n_pca_covariates is None else True
    if not os.path.isfile(parcellation_file):
        raise ValueError("Input file not found.")

    if n_pca_covariates is None:
        perform_pca = False
        pca_dict = None
    else:
        assert isinstance(
            n_pca_covariates, int
        ), "n_pca_covariates has to be an integer!"
        pca_dict = {"n_components": n_pca_covariates}
        perform_pca = True

    dist_mat = cached_distance_matrix(parcellation_file)
    null_maps = cached_null_maps(
        parcellation_file, marker_file, dist_mat, n_perm, seed=seed
    )

    specific_null_map = null_maps[:, smap_id]
    nm_parcellater = Parcellater(parcellation_file, space="MNI152")
    null_map_nifti = nm_parcellater.inverse_transform(specific_null_map)

    logger.info("Running the correlation analysis for the null map.")
    # perform correlation analysis for given surrogate map and marker
    all_genes_corr_scores, _, _ = gene_expression_correlations(
        marker=null_map_nifti,
        atlas=parcellation_file,
        aggregation_method="mean",
        allen_data_dir=allen_data_dir,
        correlation_method=correlation_method,
        perform_pca=perform_pca,
        pca_dict=pca_dict,
        partial_correlation=partial_correlation,
    )
    logger.info("Done!")

    # save correlation results for all genes for this surrogate map and these
    # settings

    name_parc = remove_nii_extensions(Path(parcellation_file).name)
    outpath = (
        Path(marker_file).parent / "nullmaps" / name_parc / "nullmaps_results"
    )
    smap_id_corr_score_str = (
        f"smapid-{smap_id}-correlationmethod-{correlation_method}"
        f"_npcacovariates-{n_pca_covariates}_seed-{seed}"
    )
    if not outpath.exists():
        outpath.mkdir(parents=True)

    outfile = outpath / f"{smap_id_corr_score_str}.tsv"
    logger.info(f"Saving null map results at {outfile}")
    all_genes_corr_scores.to_csv(outfile, sep="\t")


def _step_4(
    parcellation_file,
    marker_file,
    n_perm,
    seed,
    allen_data_dir,
    correlation_method,
    alpha,
    n_pca_covariates=None,
    r_path="Rscript",
):
    """Run step 4 in HTCondor-based pipeline.

    Perform mass-univariate correlation analysis, calculate empirical p-values
    using surrogate results, and export significant genes, as well as gene
    co-expression matrix, region-wise gene expression correlation matrix.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Path to the parcellation nifti file
    marker_file : str or os.PathLike
        Path to the marker nifti file
    n_perm : int
        number of null maps to use
    seed : int
        seed for randomness
    allen_data_dir : str or os.PathLike | None (default)
        root directory of AHBA data
    r_path : str or os.PathLike
        Rscript path at which to execute r files
    correlation_method : str
        'spearman' or 'pearson'
    alpha : float
        alpha level at which to reject the null
    n_pca_covariates : int or None
        number of components gene expression components (after pca) to include
        as covariates in the partial correlation between marker and genes.

    """
    marker_folder_output = Path(marker_file).parent / "outputs"
    perm_specific = marker_folder_output / f"nperm-{n_perm}_seed-{seed}"
    if not perm_specific.exists():
        perm_specific.mkdir(parents=True)

    # i do not use glob here because i want to be explicit about the number of
    # null maps to use
    logger.info("Loading null map results.")
    name_parc = remove_nii_extensions(Path(parcellation_file).name)
    glob_files = [
        Path(marker_file).parent
        / "nullmaps"
        / name_parc
        / "nullmaps_results"
        / (
            f"smapid-{smap_id}-correlationmethod-{correlation_method}"
            f"_npcacovariates-{n_pca_covariates}_seed-{seed}.tsv"
        )
        for smap_id in range(n_perm)
    ]

    partial_correlation = False if n_pca_covariates is None else True
    assert (
        len(glob_files) == n_perm
    ), f"{n_perm=}, and there are only {len(glob_files)} null map results!"

    # read, concat, from the smashed correlation df
    smashed_data = []
    for f in glob_files:
        smashed_data.append(pd.read_csv(f, sep="\t", index_col=0))

    smashed_corr_df = pd.concat(smashed_data, axis=1)

    # prepare output directories
    if not os.path.isfile(parcellation_file):
        raise ValueError("Input file not found.")

    # prepare potential pca
    if n_pca_covariates is None:
        perform_pca = False
        pca_dict = None
    else:
        assert isinstance(
            n_pca_covariates, int
        ), "n_pca_covariates has to be an integer!"
        pca_dict = {"n_components": n_pca_covariates}
        perform_pca = True

    logger.info(f"{perform_pca=}")
    logger.info(f"{n_pca_covariates=}")
    logger.info(f"{pca_dict=}")

    logger.info(
        "Running correlation analysis between marker and gene expressions."
    )
    (
        all_genes_corr_scores,
        _,
        covariates_dict_of_niftis,
    ) = gene_expression_correlations(
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
    real_correlations = all_genes_corr_scores[
        f"{correlation_method}_value"
    ].T.values
    smashed_correlations = smashed_corr_df[
        f"{correlation_method}_value"
    ].T.values

    logger.info("Calculating empirical p-values...")
    empirical_pvalues = empirical_pval(smashed_correlations, real_correlations)
    all_genes_corr_scores["empirical_pvals"] = empirical_pvalues
    reject, corrected, *_ = multipletests(
        all_genes_corr_scores["empirical_pvals"],
        alpha=alpha,
        method="fdr_bh",
        is_sorted=False,
        returnsorted=False,
    )
    all_genes_corr_scores["fdr_bh_corrected_empirical_pvals"] = corrected
    all_genes_corr_scores[f"reject_at_alpha-{alpha}"] = reject

    significant_genes_df = all_genes_corr_scores[reject]
    significant_genes = significant_genes_df.index.to_list()
    logger.info("Done")

    if partial_correlation and perform_pca:
        output_path_comps = (
            perm_specific
            / "pca_covariates"
            / f"{n_pca_covariates}_component_pca"
        )
        output_path = output_path_comps / correlation_method / f"alpha-{alpha}"
    else:
        output_path = perm_specific / correlation_method / f"alpha-{alpha}"

    if not os.path.isdir(output_path):
        output_path.mkdir(parents=True)

    genes_file = (
        output_path / "significant-empirical-pvalue-fdr-corrected_genes.txt"
    )
    logger.info(f"Significant genes can be found at {output_path}")
    np.savetxt(genes_file, significant_genes, fmt="%s")

    # Run gene set enrichment analysis via webgestalt R package
    run_webgestalt(r_path=r_path, genes_file=genes_file)

    if isinstance(covariates_dict_of_niftis, dict):
        for key, value in covariates_dict_of_niftis.items():
            compfile = output_path_comps / f"{key}.nii.gz"
            if not os.path.isfile(compfile):
                value.to_filename(compfile)

    genes_tsv = output_path / f"genes_{correlation_method}_and_pvalues.tsv"
    all_genes_corr_scores.to_csv(genes_tsv, sep="\t")
    for metric in ["spearman", "pearson"]:
        _save_correlation_matrices(
            parcellation_file,
            significant_genes,
            output_path,
            marker_folder_output,
            metric=metric,
        )
