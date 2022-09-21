"""Provide tests for nimgen/nimgen.py."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import sys
import tempfile

import pytest
from nilearn import datasets

from nimgen import nimgen


YAML_STRING = """project_path: "{tmp}"

pipeline:
    type: "HTCondor"
    step_1:
        CPU: 8
        MEMORY: 64
        DISK: 100
    step_2:
        CPU: 1
        MEMORY: 64
        DISK: 10
    step_3:
        CPU: 1
        MEMORY: 64
        DISK: 10

marker_dir: "{tmp}/markers"
r_path: "/r/path"
path_to_venv: "venv"

parcellation_files:
    "{parcellation_file}":
        - "{marker_file}"

n_surrogate_maps: 1000
n_pca_covariates:
    - None
    - 1
    - 3
    - 5
    - 10

correlation_method:
    - spearman
    - pearson

alpha:
    - 0.05
    - 0.01
"""


def test_cli():
    """Test CLI."""
    with tempfile.TemporaryDirectory() as tmp:
        sys.argv = [sys.argv[0]]
        sys.argv.append("--create")
        sys.argv.append(f"{tmp}/my_yaml.yaml")
        args = nimgen.parse_args()

        with pytest.raises(FileNotFoundError, match="my_yaml.yaml not found!"):
            args.create = os.path.join(tmp, "my_yaml.yaml")
            nimgen.validate_args(args)

        datasets.fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"]
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        os.mkdir(f"{tmp}/markers")
        os.system(f"touch {tmp}/markers/marker.nii")

        with open(args.create, "w") as f:
            f.write(
                YAML_STRING.format(
                    tmp=tmp,
                    parcellation_file=parc_path,
                    marker_file="marker.nii",
                )
            )

        nimgen.validate_args(args)
        nimgen.main()
