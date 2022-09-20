"""Provide tests for the nimgen/expressions.py module."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile
from itertools import product

import numpy as np
import pandas as pd
import pytest
from nilearn import image
from nilearn.datasets import fetch_atlas_schaefer_2018

from nimgen import expressions


def test_apply_pca():
    """Test apply_pca."""

    exp_data = np.random.randint(low=0, high=100, size=(100, 10))
    result = expressions.apply_pca(exp_data)
    assert result.shape == (100, 5)
    assert isinstance(result, pd.DataFrame)

    exp_data = pd.DataFrame(np.random.randint(low=0, high=100, size=(100, 10)))
    result = expressions.apply_pca(exp_data, {"n_components": 9})
    assert result.shape == (100, 9)
    assert isinstance(result, pd.DataFrame)

    for i in range(1, 10):
        assert f"covariate-{i}" in result.columns


def test_correlated_gene_expression():
    """Test correlated_gene_expression."""

    with tempfile.TemporaryDirectory() as tmp:

        image.load_img(fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"])
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        # lets create some fake gene expression data and cache it so we dont
        # have to fetch abagen data in test
        exp_data = pd.DataFrame(
            np.random.randint(low=0, high=100, size=(100, 10))
        )
        exp_data.columns = [f"gene-{x}" for x in range(10)]
        expressions._save_expressions(exp_data, parc_path)

        # test case all genes
        result = expressions.correlated_gene_expression(parc_path, "all")
        assert result.shape == (100, 100)

        # test case sign genes list
        sign_genes = [f"gene-{x}" for x in range(0, 8, 2)]
        result = expressions.correlated_gene_expression(parc_path, sign_genes)
        assert result.shape == (100, 100)

        # test case txt file
        filename = os.path.join(tmp, "file.txt")
        np.savetxt(filename, sign_genes, fmt="%s")
        result = expressions.correlated_gene_expression(parc_path, filename)
        assert result.shape == (100, 100)


def test_gene_coexpression():
    """Test gene_coexpression to obtain gene co-expression matrix."""

    with tempfile.TemporaryDirectory() as tmp:

        image.load_img(fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"])
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        # lets create some fake gene expression data and cache it so we dont
        # have to fetch abagen data in test
        exp_data = pd.DataFrame(
            np.random.randint(low=0, high=100, size=(100, 10))
        )
        exp_data.columns = [f"gene-{x}" for x in range(10)]
        expressions._save_expressions(exp_data, parc_path)

        # test case all genes
        result = expressions.gene_coexpression(parc_path, "all")
        assert result.shape == (10, 10)

        # test case sign genes list
        sign_genes = [f"gene-{x}" for x in range(0, 8, 2)]
        result = expressions.gene_coexpression(parc_path, sign_genes)
        assert result.shape == (len(sign_genes), len(sign_genes))

        # test case txt file
        filename = os.path.join(tmp, "file.txt")
        np.savetxt(filename, sign_genes, fmt="%s")
        result = expressions.gene_coexpression(parc_path, filename)
        assert result.shape == (len(sign_genes), len(sign_genes))


def test_correlation_analysis():
    """Test main mass-univariate correlation analysis."""

    with tempfile.TemporaryDirectory() as tmp:

        image.load_img(fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"])
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        marker = np.arange(100)

        # lets create some fake gene expression data and cache it so we dont
        # have to fetch abagen data in test
        exp_data = pd.DataFrame(
            np.random.randint(low=0, high=100, size=(100, 10))
        )
        exp_data.columns = [f"gene-{x}" for x in range(10)]
        expressions._save_expressions(exp_data, parc_path)

        for metric, part_corr in product(
            ["pearson", "spearman"], [True, False]
        ):

            if part_corr:
                covariates = pd.DataFrame(
                    np.random.randint(low=1, high=20, size=(100, 5))
                )

            pvals, r_scores = expressions.correlation_analysis(
                exp_data,
                marker,
                correlation_method=metric,
                partial_correlation=part_corr,
                covariates_df=covariates,
            )
            assert pvals.shape == r_scores.shape
            assert pvals.shape == (10,)

            for p_idx, r_idx in zip(pvals.index, r_scores.index):
                assert p_idx in exp_data.columns
                assert r_idx in exp_data.columns


def test_get_gene_expression():
    """Test get_gene_expression."""

    with tempfile.TemporaryDirectory() as tmp:

        parc_img = image.load_img(
            fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"]
        )
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        marker = image.new_img_like(
            parc_img,
            np.random.randint(low=1, high=40, size=parc_img.dataobj.shape),
        )

        # lets create some fake gene expression data and cache it so we dont
        # have to fetch abagen data in test
        exp_data = pd.DataFrame(
            np.random.randint(low=0, high=100, size=(100, 10))
        )
        exp_data.columns = [f"gene-{x}" for x in range(10)]
        expressions._save_expressions(exp_data, parc_path)
        all_genes, _, _ = expressions.get_gene_expression(marker, parc_path)
        assert all_genes.shape == (10, 2)
        assert "pval" in all_genes.columns
        assert "r_score" in all_genes.columns

        # test with pca covariates
        pca_dict = {"n_components": 5}
        all_genes, _, _ = expressions.get_gene_expression(
            marker,
            parc_path,
            perform_pca=True,
            pca_dict=pca_dict,
            partial_correlation=True,
        )
        assert all_genes.shape == (10, 2)
        assert "pval" in all_genes.columns
        assert "r_score" in all_genes.columns

        with pytest.warns(UserWarning, match="partial_correlation is set to "):
            # remove when implemented 'raises'
            with pytest.raises(
                NotImplementedError, match="Use of both custom"
            ):
                # test with pca covariates and custom covariates
                pca_dict = {"n_components": 5}
                all_genes, _, covar_niftis = expressions.get_gene_expression(
                    marker,
                    parc_path,
                    perform_pca=True,
                    pca_dict=pca_dict,
                    partial_correlation=False,
                    custom_covariates_df=exp_data,
                )

        with pytest.warns(UserWarning, match="perform_pca is set to "):
            # test with pca covariates
            pca_dict = {"n_components": 5}
            all_genes, _, covar_niftis = expressions.get_gene_expression(
                marker,
                parc_path,
                perform_pca=True,
                pca_dict=pca_dict,
                partial_correlation=False,
            )
            assert isinstance(covar_niftis, dict)
            assert all_genes.shape == (10, 2)
            assert "pval" in all_genes.columns
            assert "r_score" in all_genes.columns


def test_aggregate_marker():
    """Test aggregate_marker."""

    with tempfile.TemporaryDirectory() as tmp:

        parc_img = image.load_img(
            fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"]
        )
        marker = image.new_img_like(
            parc_img,
            np.random.randint(low=1, high=40, size=parc_img.dataobj.shape),
        )

        mark_agg, func_params = expressions._aggregate_marker(parc_img, marker)
        assert isinstance(mark_agg, dict)
        assert isinstance(func_params, dict)

        for key in mark_agg.keys():
            assert key in ["winsorized_mean", "mean", "std", "median"]
