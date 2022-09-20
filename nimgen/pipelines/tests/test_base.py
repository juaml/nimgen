"""Provide tests for nimgen/pipelines/base.py."""

import os
import tempfile

from nilearn import datasets

from nimgen.pipelines import base


def test_pipeline():
    """Test base pipeline."""
    with tempfile.TemporaryDirectory() as tmp:
        marker_dir = f"{tmp}/markers"
        marker_file = f"{marker_dir}/marker.nii"
        os.mkdir(marker_dir)
        os.system(f"touch {marker_file}")
        config_dict = {
            "project_path": tmp,
            "marker_dir": marker_dir,
            "parcellation_files": {
                datasets.fetch_atlas_schaefer_2018(100)["maps"]: ["marker.nii"]
            },
            "r_path": "path/to/Rscript",
            "pipeline": "local",
            "path_to_venv": "path/to/venv",
            "n_surrogate_maps": 10,
            "n_pca_covariates": ["None", 1, 3, 5, 10],
            "correlation_method": ["spearman", "pearson"],
            "alpha": [0.05, 0.01],
        }
        pipeline = base.Pipeline(**config_dict)
        parcfile = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )
        pipeline.create_output_dirs(parcfile, "marker.nii")

        param_list = [
            [parcfile],
            [parcfile, "marker.nii", 100, "spearman", 5],
            [parcfile, "marker.nii", "spearman", 0.05, 5],
        ]
        for step, params in enumerate(param_list):
            assert not pipeline._output_exists(step + 1, *params)
