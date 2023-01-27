"""Provide tests for nimgen/nimgen.py."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import shutil
import sys
import tempfile
from pathlib import Path

from nilearn import datasets

from nimgen import nimgen


YAML_STRING = """name: {name}
verbosity: INFO
seed: 100
pipeline:
    type: "HTCondor"
    step_1:
        CPU: 1
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
    step_4:
        CPU: 1
        MEMORY: 64
        DISK: 10

conda_env: "nimgen_venv"

markers:
    - path: {marker_file}
      parcellation:
        - {parcellation_file}

n_surrogate_maps: 1000

correlation_method:
    - spearman
    - pearson

alpha:
    - 0.05
"""


def test_cli():
    """Test CLI."""
    name_of_pipeline = "mockup"
    with tempfile.TemporaryDirectory() as tmp:
        sys.argv = [sys.argv[0]]
        sys.argv.append("create")
        sys.argv.append("-f")  # force overwriting previous test runs
        sys.argv.append(f"{tmp}/my_yaml.yaml")
        args = nimgen.parse_args()

        datasets.fetch_atlas_schaefer_2018(100, data_dir=tmp)["maps"]
        parc_path = os.path.join(
            tmp,
            "schaefer_2018",
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )

        marker_file = Path(tmp) / "marker.nii"
        marker_file.touch()

        with open(args.config_yaml, "w") as f:
            f.write(
                YAML_STRING.format(
                    name=name_of_pipeline,
                    tmp=tmp,
                    parcellation_file=parc_path,
                    marker_file=marker_file,
                )
            )

        nimgen.main()

        mockup_path = Path(name_of_pipeline)
        all_gene_output_path = mockup_path / "all_gene_outputs"
        marker_path = mockup_path / "markers" / "marker" / "marker.nii"
        parc_path = (
            mockup_path
            / "parcellations"
            / "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm"
            / "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz"
        )
        yaml_path = mockup_path / "my_yaml.yaml"
        log_path = mockup_path / "submit_files" / "logs"

        assert mockup_path.exists()
        assert all_gene_output_path.exists()
        assert marker_path.exists()
        assert parc_path.exists()
        assert yaml_path.exists()
        assert log_path.exists()

        shutil.rmtree(name_of_pipeline)
