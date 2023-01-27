"""Provide tests for the nimgen/smash.py module interfacing with brainsmash."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import tempfile
from pathlib import Path

import numpy as np
from nibabel import Nifti1Image

from nimgen import smash


def _make_atlas(size, seed=5):
    return Nifti1Image(
        np.random.default_rng(seed=seed).integers(low=1, high=5, size=size),
        np.eye(4),
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


def test_cached_distance_matrix():
    """Test cached_distance_matrix."""
    parcellation = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        parcellation_path = Path(tmp) / "parcellation.nii.gz"
        parcellation.to_filename(parcellation_path)
        new = smash.cached_distance_matrix(parcellation_path)

        # check if caching worked correctly
        assert (Path(tmp) / "parcellation_dist_mat.npy").is_file()
        cached = smash.cached_distance_matrix(parcellation_path)
        np.testing.assert_array_equal(new, cached)


def test_cached_null_maps():
    """Test cached_null_maps."""
    parcellation = _make_atlas(size=(5, 5, 2))
    marker = _make_atlas(size=(5, 5, 2), seed=10)

    with tempfile.TemporaryDirectory() as tmp:
        parcellation_path = Path(tmp) / "parcellation.nii.gz"
        marker_path = Path(tmp) / "marker.nii.gz"
        null_maps_dir = Path(tmp) / "nullmaps" / "parcellation"
        null_maps_file = (
            null_maps_dir / "marker-marker_perms-100_seed-10_desc-nullmaps.npy"
        )

        parcellation.to_filename(parcellation_path)
        marker.to_filename(marker_path)

        # generating null maps with my made up marker data
        # fails
        # TODO: make generating null maps work fast with some
        # made up data for testing
        null_maps = np.random.default_rng(seed=10).random((100, 5))

        null_maps_dir.mkdir(parents=True, exist_ok=True)
        np.save(null_maps_file, null_maps)

        dist_mat = np.random.default_rng(seed=10).random((5, 5))
        dist_mat = (dist_mat + dist_mat.T) / 2

        null_maps_cached = smash.cached_null_maps(
            parcellation_path,
            marker_path,
            dist_mat,
            n_perm=100,
            seed=10,
        )
        np.testing.assert_array_equal(null_maps, null_maps_cached)


test_cached_null_maps()
