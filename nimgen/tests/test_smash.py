"""Provide tests for the nimgen/smash.py module interfacing with brainsmash."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import tempfile
from pathlib import Path

import numpy as np
from nibabel import Nifti1Image

from nimgen import smash


# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL


def _make_atlas(size):
    return Nifti1Image(
        np.random.default_rng(seed=5).integers(low=1, high=5, size=size),
        np.eye(4),
        dtype="int32",
    )


def test_vox_dist():
    """Test faster vox_dist implementation."""

    atlas = _make_atlas(size=(5, 5, 2))
    original_implementation = smash._vox_dist_original(atlas)
    fast_implementation = smash.vox_dist(atlas)

    rows, cols = original_implementation.shape
    assert rows == cols

    off_diag_mask = ~np.eye(rows, dtype=bool)
    original_off_diag = original_implementation[off_diag_mask]
    fast_off_diag = fast_implementation[off_diag_mask]

    np.testing.assert_almost_equal(original_off_diag, fast_off_diag, decimal=6)


def test_cache_distance_matrix():
    """Test cache_distance_matrix."""
    parcellation = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        parcellation_path = Path(tmp) / "parcellation.nii.gz"
        parcellation.to_filename(parcellation_path)
        new = smash.cached_distance_matrix(parcellation_path)

        # check if caching worked correctly
        assert (Path(tmp) / "parcellation_dist_mat.npy").is_file()
        cached = smash.cached_distance_matrix(parcellation_path)
        np.testing.assert_array_equal(new, cached)
