"""Define base pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import shutil
from pathlib import Path

from ..utils import remove_nii_extensions


class Pipeline:
    """Object to create and run a nimgen pipeline.

    Parameters
    -----------
    config_dict : dict
        dictionary with valid pipeline configurations.
    """

    def __init__(
        self,
        pipeline,
        markers,
        seed=1,
        name="nimgen",
        r_path=None,
        path_to_venv=None,
        conda_env=None,
        n_surrogate_maps=100,
        correlation_method=None,
        alpha=0.05,
        n_pca_covariates=None,
        output_dir="output",
        config_dict=None,
        ahba_dir=None,
    ):
        """Initialise pipeline directory."""
        self.pipeline = pipeline
        self.output_dir = output_dir
        self.config_dict = config_dict
        self.r_path = r_path
        self.path_to_venv = path_to_venv
        self.conda_env = conda_env
        self.n_surrogate_maps = n_surrogate_maps
        self._init_dir = Path(os.getcwd())
        self.seed = int(seed)
        self.name = name

        # take care of virtual environment to run in
        if self.path_to_venv is not None:
            self.conda = False
            self.path_to_venv = os.path.abspath(self.path_to_venv)
        else:
            assert conda_env is not None
            self.conda_env = self.conda_env
            self.conda = True

        if r_path is None:
            self.r_path = "/usr/bin/env Rscript"

        # correlation analysis parameters
        if correlation_method is None:
            self.correlation_method = ["spearman"]
        else:
            if isinstance(correlation_method, str):
                correlation_method = [correlation_method]
            assert isinstance(
                correlation_method, list
            ), "Correlation method should be provided as a list or a string!"
            self.correlation_method = correlation_method

        if not isinstance(alpha, list):
            alpha = [alpha]

        self.alpha = alpha

        # potential (still experintal) pca removal
        if not isinstance(n_pca_covariates, list):
            n_pca_covariates = [n_pca_covariates]
        elif n_pca_covariates == 0 or n_pca_covariates is None:
            n_pca_covariates = [None]
        self.n_pca_covariates = n_pca_covariates

        # prepare markers (and parcellations)
        self.markers = markers

        self.jobs_dir = Path(".") / f"{name}_jobs"
        self.markers_dir = self.jobs_dir / "markers"
        self.parcellations_dir = self.jobs_dir / "parcellations"
        self.all_gene_outputs = self.jobs_dir / "all_gene_outputs"
        self.ahba_dir = ahba_dir
        if self.ahba_dir is None:
            self.ahba_dir = self.jobs_dir / "AHBA"

    def _prepare_marker_dir(self, marker):
        """Prepare the marker directory."""
        marker_path = Path(marker["path"])
        marker_name = Path(remove_nii_extensions(marker_path)).stem
        marker_dir = self.markers_dir / marker_name
        marker_dir.mkdir()
        marker_dst = marker_dir / marker_path.name
        marker_path = marker_path.resolve()

        shutil.copy(marker_path.as_posix(), marker_dst.as_posix())
        marker_parcellations = marker["parcellation"]

        if isinstance(marker_parcellations, str):
            marker_parcellations = [marker_parcellations]

        parc_names = [remove_nii_extensions(x) for x in marker_parcellations]

        # null maps directories
        null_maps_dir = marker_dir / "nullmaps"
        null_maps_dir.mkdir()

        for name in parc_names:
            parc_dir = null_maps_dir / Path(name).name
            parc_dir.mkdir()
            output_nullmaps = parc_dir / "nullmaps_results"
            output_nullmaps.mkdir()

        # outpath
        outpath = marker_dir / "outputs"
        outpath.mkdir()

        return marker_dst, marker_parcellations

    def _prepare_parcellation_dir(self, parcellation_path):
        """Prepare the parcellation output directory."""
        name = remove_nii_extensions(parcellation_path)
        parcellation_dir = self.parcellations_dir / Path(name).name
        parcellation_dir.mkdir()
        shutil.copy(parcellation_path, parcellation_dir)
        return parcellation_dir / Path(parcellation_path).name

    def create_jobs_dir(self):
        """Create the pipeline directory."""

        if self.jobs_dir.exists():
            raise FileExistsError(
                "Pleave remove or rename existing nimgen_jobs directories."
            )
        self.jobs_dir.mkdir()
        self.markers_dir.mkdir()
        self.parcellations_dir.mkdir()
        self.all_gene_outputs.mkdir()

        # create the marker dirs and obtain list of associated parcellations
        all_parcellations = []
        for idx, marker in enumerate(self.markers):
            new_marker_path, marker_parcellations = self._prepare_marker_dir(
                marker
            )

            # create parcellation-only dirs (i.e. for distance matrices)
            marker_specific_parcs = []
            for parc in marker_parcellations:
                marker_specific_parcs += [
                    self._prepare_parcellation_dir(parc).as_posix()
                ]

            all_parcellations += marker_specific_parcs.copy()
            self.markers[idx]["parcellation"] = marker_specific_parcs
            self.markers[idx]["path"] = new_marker_path

        self.all_parcellations = list(set(all_parcellations))

    def create(self):
        """Create the base pipeline."""
        self.create_jobs_dir()
