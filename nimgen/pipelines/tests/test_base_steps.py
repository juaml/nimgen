"""Provide tests for basic pipeline steps."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile

import numpy as np
import pandas as pd
from nibabel import Nifti1Image

from nimgen.expressions import _save_expressions
from nimgen.pipelines import base_steps, htcondor


def _make_atlas(size):
    return Nifti1Image(
        np.arange(50).reshape(5, 5, 2),
        np.eye(4),
    )


def _make_marker_nii(size):
    return Nifti1Image(
        np.random.default_rng(seed=10).integers(low=1, high=5, size=size),
        np.eye(4),
    )


def _make_marker_dir(tmp):
    marker_dir = f"{tmp}/markers"
    marker_file = f"{marker_dir}/marker.nii"
    os.mkdir(marker_dir)
    marker_nii = _make_marker_nii(size=(5, 5, 2))
    marker_nii.to_filename(marker_file)
    return marker_dir, marker_file


def _make_pipeline_configs(tmp):
    marker_dir, marker_file = _make_marker_dir(tmp)
    config_dict = {
        "project_path": tmp,
        "marker_dir": marker_dir,
        "parcellation_files": {f"{tmp}/test_atlas.nii": ["marker.nii"]},
        "r_path": "path/to/Rscript",
        "pipeline": {
            "type": "HTCondor",
            "step_1": {
                "CPU": 8,
                "MEMORY": 64,
                "DISK": 100,
            },
            "step_2": {"CPU": 1, "MEMORY": 64, "DISK": 10},
            "step_3": {"CPU": 1, "MEMORY": 64, "DISK": 10},
        },
        "path_to_venv": "path/to/venv",
        "n_surrogate_maps": 10,
        "n_pca_covariates": ["None", 1, 3, 5, 10],
        "correlation_method": ["spearman", "pearson"],
        "alpha": [0.05, 0.01],
    }
    return config_dict, marker_dir, marker_file


def test_steps():
    """Test each step together since earlier steps set up for later steps."""
    # set up and save a fake atlas as nifti
    atlas = _make_atlas(size=(5, 5, 2))
    with tempfile.TemporaryDirectory() as tmp:
        atlas_path = os.path.join(tmp, "test_atlas.nii")
        atlas.to_filename(atlas_path)
        configs, marker_dir, marker_file = _make_pipeline_configs(tmp)
        pipeline = htcondor.HTCondor(**configs)
        pipeline.create()

        parcellation = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "test_atlas",
            "test_atlas.nii",
        )
        output_dir = os.path.join(pipeline.project_path, pipeline.output_dir)

        # we also need to cache some fake gene expression data so it wont take
        # forever
        smaps_dir = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "test_atlas",
            "smaps",
        )
        base_steps.step_1(parcellation)

        for i in range(1, 4):

            smap_file = os.path.join(smaps_dir, f"{i}_smap.nii")
            exp = pd.DataFrame(
                np.random.randint(low=5, high=100, size=(49, 10))
            )
            _save_expressions(exp, smap_file)
            base_steps.step_2(
                parcellation,
                marker_file,
                marker_dir,
                output_dir,
                smap_id=i,
                allen_data_dir=tmp,
                knn=20,
                ns=5,
                pv=50,
            )

        # webgestalt will fail due to fake r path but the python program will
        # continue after
        exp = pd.DataFrame(np.random.randint(low=5, high=100, size=(49, 10)))
        _save_expressions(exp, parcellation)
        base_steps.step_3(
            parcellation,
            marker_file,
            marker_dir,
            output_dir,
            r_path="this/is/not/real",
        )
