"""Provide tests for the nimgen/expressions.py module."""

import os
import tempfile

import numpy as np
import pandas as pd
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
        fname_fake_cached_expressions = os.path.join(
            tmp,
            "schaefer_2018",
            "nimgen_Schaefer2018_100Parcels_"
            "7Networks_order_FSLMNI152_1mm_expressions.csv",
        )

        exp_data = pd.DataFrame(
            np.random.randint(low=0, high=100, size=(100, 10))
        )
        exp_data.columns = [f"gene-{x}" for x in range(10)]
        exp_data.to_csv(fname_fake_cached_expressions, sep=";", index=False)

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
