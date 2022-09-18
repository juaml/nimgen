"""Tests for the nimgen/utils.py module."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile
import numpy as np
import pandas as pd
import pytest

from nimgen import utils


def test_remove_nii_extensions():
    """Test remove_nii_extensions."""
    filenames = ["file1.nii", "file1.nii.gz"]
    for f in filenames:
        result = utils.remove_nii_extensions(f)
        assert result == "file1"


def test_read_csv_tsv():
    """Test read_csv_tsv."""
    mock_df = pd.DataFrame(
        {
            "col1": [x for x in range(5)],
            "col2": [x for x in range(5)]
        }
    )

    # test csv
    with tempfile.TemporaryDirectory() as tmp:
        csv_file = os.path.join(tmp, "file.csv")
        mock_df.to_csv(csv_file)

        df = utils.read_csv_tsv(csv_file)
        assert df.shape == (5, 2)

    # test tsv
    with tempfile.TemporaryDirectory() as tmp:
        tsv_file = os.path.join(tmp, "file.tsv")
        mock_df.to_csv(tsv_file, sep="\t")

        df = utils.read_csv_tsv(tsv_file)
        assert df.shape == (5, 2)


def test_read_sign_genes():
    """Test read_sign_genes."""
    sign_genes = [f"gene{x}" for x in range(5)]
    mock_df = pd.DataFrame(
        {
            "col1": [x for x in range(5)],
            "col2": [x for x in range(5)]
        }
    )
    mock_df.index = sign_genes

    result = utils._read_sign_genes(mock_df)
    assert result == sign_genes

    # test csv
    with tempfile.TemporaryDirectory() as tmp:
        csv_file = os.path.join(tmp, "file.csv")
        mock_df.to_csv(csv_file, header=None)
        result = utils._read_sign_genes(csv_file)
        assert result == sign_genes

    # test tsv
    with tempfile.TemporaryDirectory() as tmp:
        tsv_file = os.path.join(tmp, "file.tsv")
        mock_df.to_csv(tsv_file, sep="\t", header=None)
        result = utils._read_sign_genes(tsv_file)
        assert result == sign_genes

    # test txt
    with tempfile.TemporaryDirectory() as tmp:
        txt_file = os.path.join(tmp, "file.txt")
        np.savetxt(txt_file, sign_genes, fmt="%s")
        result = utils._read_sign_genes(txt_file)
        assert result == sign_genes

    with pytest.raises(ValueError, match="'sign_genes' should be a pd."):
        result = utils._read_sign_genes(4)
