"""Define base pipeline based on an abstract base class."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
from itertools import product
from ast import literal_eval
from pathlib import Path
import shutil
from ..utils import remove_nii_extensions


def _specific_marker_output(
    marker_file, marker_dir, output_dir, name_parc
):
    # step 1 is taken care of by the step itsel
    # step 2:
    marker_file = os.path.join(marker_dir, marker_file)
    _, marker_structure = Path(
        marker_file
    ).absolute().as_posix().split(marker_dir)
    head, tail = os.path.split(marker_structure[1:])
    marker_name, _ = os.path.splitext(tail)
    return os.path.join(
        output_dir, name_parc, head, marker_name
    )


class Pipeline:
    """Object to create and run a nimgen pipeline.

    Parameters
    -----------
    config_dict : dict
        dictionary with valid pipeline configurations.
    """

    def __init__(
        self,
        project_path,
        marker_dir,
        pipeline,
        r_path,
        path_to_venv,
        parcellation_files,
        n_surrogate_maps=100,
        correlation_method=None,
        alpha=0.05,
        n_pca_covariates=None,
        pipeline_dir="pipeline",
        output_dir="output",
        parcellations_dir="parcellations",
        config_dict=None,

    ):
        """Initialise pipeline directory."""
        self.project_path = os.path.abspath(project_path)

        self.pipeline = pipeline
        self.marker_dir = marker_dir
        self.pipeline_dir = pipeline_dir
        self.output_dir = output_dir
        self.parcellations_dir = "parcellations"
        self.parcellation_files = parcellation_files
        self.config_dict = config_dict
        self.parcellation_marker_dict = {}
        self.r_path = os.path.abspath(r_path)
        self.path_to_venv = path_to_venv
        self.n_surrogate_maps = n_surrogate_maps

        if correlation_method is None:
            self.correlation_method = ["spearman"]
        else:
            if isinstance(correlation_method, str):
                correlation_method = [correlation_method]
            assert isinstance(correlation_method, list), (
                "Correlation method should be provided as a list or a string!"
            )
            self.correlation_method = correlation_method

        if not isinstance(alpha, list):
            alpha = [alpha]

        self.alpha = alpha

        if not isinstance(n_pca_covariates, list):
            n_pca_covariates = [n_pca_covariates]
        elif n_pca_covariates == 0:
            n_pca_covariates = None
        self.n_pca_covariates = n_pca_covariates

        self.allen_data_dir = os.path.join(self.project_path, "allen_data")
        if not os.path.isdir(self.project_path):
            raise FileNotFoundError(f"{self.project_path} not found!")

        if not os.path.isdir(self.marker_dir):
            raise FileNotFoundError(f"{self.marker_dir} not found!")

        for dir in [
            self.output_dir, self.parcellations_dir, self.pipeline_dir
        ]:
            dir = os.path.join(self.project_path, dir)
            if not os.path.isdir(dir):
                print(f"Creating {dir}")
                os.mkdir(dir)
            else:
                print(f"{dir} already exists! Skipping creation of dir {dir}")

        for parcellation_file in self.parcellation_files.keys():
            if not os.path.isfile(parcellation_file):
                FileNotFoundError(f"{parcellation_file} not found!")
            else:
                _, parc_file_tail = os.path.split(parcellation_file)
                parcellation_name = remove_nii_extensions(parc_file_tail)
                self.parcellation_marker_dict[parcellation_name] = (
                    parc_file_tail,
                    self.parcellation_files[parcellation_file]
                )
                dir_this_parc = os.path.join(
                    self.parcellations_dir, parcellation_name
                )
                smaps_dir = os.path.join(dir_this_parc, "smaps")
                if not os.path.isdir(dir_this_parc):
                    print(f"Creating {dir_this_parc}")
                    os.mkdir(dir_this_parc)
                    os.mkdir(smaps_dir)
                else:
                    print(
                        f"{dir_this_parc} already exists!"
                        f" Skipping creation of dir {dir}"
                    )
                    if not os.path.isdir(smaps_dir):
                        os.mkdir(smaps_dir)

                parc_file_new = os.path.join(dir_this_parc, parc_file_tail)
                if not os.path.isfile(parc_file_new):
                    print(
                        f"Copying from {parcellation_file} to {dir_this_parc}"
                    )
                    shutil.copyfile(parcellation_file, parc_file_new)
                else:
                    print(f"{parc_file_new} already exists. Skipping copy.")

        for _, (_, markers) in self.parcellation_marker_dict.items():
            for marker_file in markers:
                current_marker = os.path.join(self.marker_dir, marker_file)
                if not os.path.isfile(current_marker):
                    raise FileNotFoundError(
                        "Marker files are interpreted relative to marker_dir."
                        f"({self.marker_dir})"
                    )

    def create_output_dirs(self, name_parc, marker_file):
        """Create output directory structure for a marker and parcellation."""
        # step 1 is taken care of by the step itsel
        # step 2:
        output_dir = os.path.join(self.project_path, self.output_dir)
        specific_marker_output = _specific_marker_output(
            marker_file, self.marker_dir, output_dir, name_parc
        )
        smap_corr_scores = os.path.join(
            specific_marker_output, "smap_corr_scores"
        )
        if not os.path.isdir(smap_corr_scores):
            os.makedirs(smap_corr_scores)

        # step 3:
        for n_pca_cov, correlation_method, alpha in product(
            self.n_pca_covariates, self.correlation_method, self.alpha
        ):
            if isinstance(n_pca_cov, int):
                output_path_comps = os.path.join(
                    specific_marker_output, "pca_covariates",
                    f"{n_pca_cov}_component_pca"
                )
                output_path = os.path.join(
                    output_path_comps, correlation_method, f"alpha-{alpha}"
                )
            elif literal_eval(n_pca_cov) is None:
                output_path = os.path.join(
                    specific_marker_output,
                    correlation_method, f"alpha-{alpha}"
                )
            if not os.path.isdir(output_path):
                os.makedirs(output_path)

    def _output_exists(self, step, *step_params):
        check_funcs = [
            self._output_step_1_exists,
            self._output_step_2_exists,
            self._output_step_3_exists,
        ]
        return check_funcs[int(step) - 1](*step_params)

    def _output_step_1_exists(self, parcfile):

        path, _ = os.path.split(parcfile)
        files_exist = [os.path.split(x)[1] for x in os.listdir(path)]
        files_needed = [
            "brain_map.txt", "distmat.py", "index.npy", "voxel_coordinates.txt"
        ]
        for f in files_needed:
            if f not in files_exist:
                return False
        return True

    def _output_step_2_exists(
        self, parc_file, marker_file, smapid, corr_method, n_pca_cov
    ):

        _, parc_file_tail = os.path.split(parc_file)
        name_parc, _ = os.path.splitext(parc_file_tail)
        output_dir = os.path.join(self.project_path, self.output_dir)
        specific_marker_output = _specific_marker_output(
            marker_file, self.marker_dir, output_dir, name_parc
        )
        smap_corr_scores = os.path.join(
            specific_marker_output, "smap_corr_scores"
        )

        fname = (
            f"smapid-{smapid}-correlationmethod-{corr_method}_"
            f"npcacovariates-{n_pca_cov}.tsv"
        )
        check_file = os.path.join(smap_corr_scores, fname)
        return os.path.isfile(check_file)

    def _output_step_3_exists(
        self, parc_file, marker_file, corr_method, alpha, n_pca_cov
    ):

        _, parc_file_tail = os.path.split(parc_file)
        name_parc, _ = os.path.splitext(parc_file_tail)

        output_dir = os.path.join(self.project_path, self.output_dir)
        specific_marker_output = _specific_marker_output(
            marker_file, self.marker_dir, output_dir, name_parc
        )

        if isinstance(n_pca_cov, int):
            output_path_comps = os.path.join(
                specific_marker_output, "pca_covariates",
                f"{n_pca_cov}_component_pca"
            )
            output_path = os.path.join(
                output_path_comps, corr_method, f"alpha-{alpha}"
            )
        elif literal_eval(n_pca_cov) is None:
            output_path = os.path.join(
                specific_marker_output, corr_method, f"alpha-{alpha}"
            )

        genes_file = os.path.join(
            output_path,
            'significant-empirical-pvalue-fdr-corrected_genes.txt'
        )
        if not os.path.isfile(genes_file):
            return False

        for metric in ["spearman", "pearson"]:
            matrix_files = [
                (
                    'significant-genes_gene_by_gene_'
                    f'correlation_matrix_{metric}.tsv'
                ),
                (
                    'significant-genes_region_by_region_'
                    f'correlation_matrix_{metric}.tsv'
                ),
                f'all-genes_gene_by_gene_correlation_matrix_{metric}.tsv',
                f'all-genes_region_by_region_correlation_matrix_{metric}.tsv',
            ]
            for f in matrix_files:
                filename = os.path.join(output_path, f)
                if not os.path.isfile(filename):
                    return False

        return True
