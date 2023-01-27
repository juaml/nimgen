"""Provide tests for the HTCondor pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import shutil
import tempfile
from pathlib import Path

from nilearn import datasets

from nimgen.pipelines import htcondor


def _make_marker_dir(tmp):
    marker_dir = f"{tmp}/markers"
    marker_file = f"{marker_dir}/marker.nii"
    os.mkdir(marker_dir)
    os.system(f"touch {marker_file}")
    return marker_file


def _make_pipeline_configs(tmp):
    parc = datasets.fetch_atlas_schaefer_2018(100)["maps"]
    config_dict = {
        "name": "mockup",
        "verbosity": "WARNING",
        "seed": 100,
        "markers": [{"path": _make_marker_dir(tmp), "parcellation": [parc]}],
        "pipeline": {
            "type": "HTCondor",
            "step_1": {"CPU": 1, "MEMORY": 64, "DISK": 10},
            "step_2": {"CPU": 1, "MEMORY": 64, "DISK": 10},
            "step_3": {"CPU": 1, "MEMORY": 64, "DISK": 10},
            "step_4": {"CPU": 1, "MEMORY": 64, "DISK": 10},
        },
        "conda_env": "conda_env_name_mockup",
        "n_surrogate_maps": 10,
        "n_pca_covariates": ["None", 1, 3, 5, 10],
        "correlation_method": ["spearman", "pearson"],
        "alpha": [0.05, 0.01],
    }
    return config_dict


def test_pipeline_create():
    """Test create method for any failure."""

    with tempfile.TemporaryDirectory() as tmp:

        config_dict = _make_pipeline_configs(tmp)
        pipeline = htcondor.HTCondor(**config_dict)
        pipeline.create()
        shutil.rmtree(Path.cwd() / config_dict["name"])
