"""Provide tests for the HTCondor pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import tempfile

from nilearn import datasets

from nimgen.pipelines import htcondor


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
    return config_dict


def test_pipeline_individual_methods():
    """Test base pipeline with all params set by users."""
    with tempfile.TemporaryDirectory() as tmp:

        config_dict = _make_pipeline_configs(tmp)
        pipeline = htcondor.HTCondor(**config_dict)
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
        submit_files_dir = os.path.join(
            pipeline.project_path, pipeline.submit_files_dir
        )
        assert os.path.isdir(parc_dir)
        assert os.path.isdir(pipeline_dir)
        assert os.path.isdir(submit_files_dir)

        for step in range(1, 4):
            pipeline.prepare_step(step)
            check_file = os.path.join(
                pipeline.project_path, pipeline.pipeline_dir, f"step_{step}.py"
            )
            assert os.path.isfile(check_file)

        pipeline.prepare_run_in_venv()
        check_file = os.path.join(
            pipeline.project_path, pipeline.pipeline_dir, "run_in_venv.sh"
        )
        assert os.path.isfile(check_file)
        pipeline.prepare_submit_files()

        submit_run_dir = os.path.join(
            pipeline.project_path,
            pipeline.submit_files_dir,
            "Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm",
            "run_1",
        )
        for step in range(1, 4):
            check_file = os.path.join(submit_run_dir, f"step_{step}.submit")
            assert os.path.isfile(check_file)


def test_pipeline_create():
    """Test create method for any failure."""

    with tempfile.TemporaryDirectory() as tmp:

        config_dict = _make_pipeline_configs(tmp)
        pipeline = htcondor.HTCondor(**config_dict)
        pipeline.create()
