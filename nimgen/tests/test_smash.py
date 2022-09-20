"""Provide tests for the nimgen/smash.py module interfacing with brainsmash."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile

import numpy as np
from nibabel import Nifti1Image
from nilearn._utils import check_niimg

from nimgen import smash


def _make_atlas(size):
    return Nifti1Image(
        np.random.default_rng(seed=5).integers(low=1, high=5, size=size),
        np.eye(4),
    )


def test_export_voxel_coordinates():
    """Test export_voxel_coordinates."""

    atlas = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        coord_file, parcel_file = smash.export_voxel_coordinates(atlas, tmp)
        assert os.path.isfile(coord_file)
        assert os.path.isfile(parcel_file)

        _, name_coord_file = os.path.split(coord_file)
        _, name_parcel_file = os.path.split(parcel_file)

        assert name_coord_file == "voxel_coordinates.txt"
        assert name_parcel_file == "brain_map.txt"


def test_generate_distance_matrices():
    """Test generate_distance_matrices."""

    atlas = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        _, _ = smash.export_voxel_coordinates(atlas, tmp)
        filenames = smash.generate_distance_matrices(tmp)
        for key, value in filenames.items():
            assert key in ["D", "index"]
            assert os.path.isfile(value)
            _, filename = os.path.split(value)
            assert filename in ["index.npy", "distmat.npy"]

        # test case where distance matrix already exists
        filenames = smash.generate_distance_matrices(tmp)


def test_generate_surrogate_map():
    """Test generate_surrogate_map."""

    atlas = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        _, parcel_file = smash.export_voxel_coordinates(atlas, tmp)
        filenames = smash.generate_distance_matrices(tmp)
        smap_file = smash.generate_surrogate_map(
            atlas, 5, tmp, parcel_file, filenames, knn=20, ns=5, pv=50
        )
        assert os.path.isfile(smap_file)
        check_niimg(smap_file)

        # test case where smaps_file already exists
        smap_file = smash.generate_surrogate_map(
            atlas, 5, tmp, parcel_file, filenames, knn=20, ns=5, pv=50
        )
