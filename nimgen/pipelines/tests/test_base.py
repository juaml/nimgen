"""Provide tests for nimgen/pipelines/base.py."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile

from nilearn import datasets

from nimgen.pipelines import base


def _make_marker_dir(tmp):
    marker_dir = f"{tmp}/markers"
    marker_file = f"{marker_dir}/marker.nii"
    os.mkdir(marker_dir)
    os.system(f"touch {marker_file}")
    return marker_dir, marker_file


def _make_pipeline_configs(tmp):
    marker_dir, _ = _make_marker_dir(tmp)
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
    return config_dict


def _make_pipeline_configs_defaults_and_nones(tmp):
    marker_dir, _ = _make_marker_dir(tmp)
    config_dict = {
        "project_path": tmp,
        "marker_dir": marker_dir,
        "parcellation_files": {
            datasets.fetch_atlas_schaefer_2018(100)["maps"]: ["marker.nii"]
        },
        "r_path": "path/to/Rscript",
        "pipeline": "local",
        "path_to_venv": "path/to/venv",
    }
    return config_dict


def test_pipeline_all_params_set():
    """Test base pipeline with all params set by users."""
    with tempfile.TemporaryDirectory() as tmp:

        config_dict = _make_pipeline_configs(tmp)
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

        parc_dir = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm",
        )
        pipeline_dir = os.path.join(
            pipeline.project_path, pipeline.pipeline_dir
        )
        assert os.path.isdir(parc_dir)
        assert os.path.isdir(pipeline_dir)

        # test again but this time test if pipeline directory already exists
        pipeline = base.Pipeline(**config_dict)
        parcfile = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz",
        )
        output_path = pipeline.create_output_dirs(parcfile, "marker.nii")
        param_list = [
            [parcfile],
            [parcfile, "marker.nii", 100, "spearman", 5],
            [parcfile, "marker.nii", "spearman", 0.05, 5],
        ]
        for step, params in enumerate(param_list):
            assert not pipeline._output_exists(step + 1, *params)

        # check in case step 3 output partly exists
        specific_marker_output = base._specific_marker_output(
            "marker.nii",
            pipeline.marker_dir,
            output_path,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm",
        )
        genes_file = os.path.join(
            specific_marker_output,
            "spearman",
            f"alpha-{0.05}",
            "significant-empirical-pvalue-fdr-corrected_genes.txt",
        )
        os.system(f"touch {genes_file}")
        param_list = [
            [parcfile],
            [parcfile, "marker.nii", 100, "spearman", 5],
            [parcfile, "marker.nii", "spearman", 0.05, 5],
        ]
        for step, params in enumerate(param_list):
            assert not pipeline._output_exists(step + 1, *params)

        # check that pipeline dir etc exist
        parc_dir = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm",
        )
        pipeline_dir = os.path.join(
            pipeline.project_path, pipeline.pipeline_dir
        )
        assert os.path.isdir(parc_dir)
        assert os.path.isdir(pipeline_dir)


def test_pipeline_with_mainly_defaults():
    """Test pipeline with minimum cofiguration and mainly default values."""
    with tempfile.TemporaryDirectory() as tmp:

        config_dict = _make_pipeline_configs_defaults_and_nones(tmp)
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

        parc_dir = os.path.join(
            pipeline.project_path,
            pipeline.parcellations_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm",
        )
        pipeline_dir = os.path.join(
            pipeline.project_path, pipeline.pipeline_dir
        )
        assert os.path.isdir(parc_dir)
        assert os.path.isdir(pipeline_dir)
