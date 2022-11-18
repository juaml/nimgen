"""Define HTCondor pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
from itertools import product

from nimgen.utils import remove_nii_extensions

from ._htcondor_python_strings import (
    RUN_CONDA,
    RUN_IN_VENV,
    STEP_ONE_FSTRING,
    STEP_THREE_FSTRING,
    STEP_TWO_FSTRING,
)
from ._htcondor_submit_strings import (
    TEMPLATE_DAG,
    TEMPLATE_JOB,
    TEMPLATE_QUEUED_JOB,
)
from .base import Pipeline


class HTCondor(Pipeline):
    """Object to create and run a nimgen pipeline on a HTCondor cluster.

    Parameters
    -----------
    config_dict : dict
        dictionary with valid pipeline configurations.
    """

    def __init__(self, submit_files_dir="submit_files", *args, **kwargs):
        """Initialise HTCondor pipeline."""
        super().__init__(*args, **kwargs)
        self.submit_files_dir = submit_files_dir
        directory = os.path.join(self.project_path, self.submit_files_dir)
        if not os.path.isdir(directory):
            print(f"Creating {directory}")
            os.mkdir(directory)
        else:
            print(
                f"{directory} already exists!"
                f"Skipping creation of dir {directory}"
            )

    def prepare_run_in_venv(self):
        """Prepare a bash file to source a venv or conda env."""

        run_in_venv = os.path.join(
            self.project_path, self.pipeline_dir, "run_in_venv.sh"
        )

        if self.conda:
            with open(run_in_venv, "w") as f:
                f.write(RUN_CONDA.format(self.path_to_conda_env))
        else:
            with open(run_in_venv, "w") as f:
                f.write(RUN_IN_VENV.format(self.path_to_venv))

        os.system(f"chmod +x {run_in_venv}")

    def prepare_step(self, step):
        """Prepare a python3 file for a given step in the pipeline."""
        _, name_marker_dir = os.path.split(self.marker_dir)
        output_dir = os.path.join(self.project_path, self.output_dir)
        allen_data_dir = self.allen_data_dir
        marker_dir = os.path.join(self.project_path, self.marker_dir)

        step_args = [
            (STEP_ONE_FSTRING, ["placeholder"]),
            (STEP_TWO_FSTRING, [marker_dir, output_dir, allen_data_dir]),
            (
                STEP_THREE_FSTRING,
                [marker_dir, output_dir, allen_data_dir, self.r_path],
            ),
        ]
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_file = os.path.join(pipeline_dir, f"step_{step}.py")
        if not os.path.isfile(step_file):
            with open(step_file, "w") as f:
                fstring, args = step_args[step - 1]
                f.write(fstring.format(*args))

    def prepare_submit_files(self):
        """Prepare submit files to submit the pipeline to HTCondor as a DAG."""

        for parcellation, (
            _,
            markers,
        ) in self.parcellation_marker_dict.items():

            submit_parc_dir = self._prepare_submit_parc_dir(parcellation)
            run_dir = self._prepare_run_dir(submit_parc_dir)
            self._prepare_dag(run_dir)
            for marker in markers:
                print(parcellation, marker)
                self.create_output_dirs(parcellation, marker)

            for step in range(1, 4):
                job_count = self._prepare_submit_file(
                    step, run_dir, parcellation
                )
                print(f"{job_count} jobs registered for step {step}.")

    def _prepare_submit_file(self, step, run_dir, parcellation):

        logs = os.path.abspath(os.path.join(run_dir, "logs"))
        initial_dir = os.path.join(self.project_path, self.pipeline_dir)
        n_cpus = self.pipeline[f"step_{step}"]["CPU"]
        memory = self.pipeline[f"step_{step}"]["MEMORY"]
        disk = self.pipeline[f"step_{step}"]["DISK"]
        parc_file_name = self.parcellation_marker_dict[parcellation][0]
        marker_files = [
            x for x in self.parcellation_marker_dict[parcellation][1]
        ]
        parcellation_file = os.path.join(
            self.project_path,
            self.parcellations_dir,
            parcellation,
            parc_file_name,
        )

        submit_file = os.path.join(run_dir, f"step_{step}.submit")
        with open(submit_file, "w") as f:
            f.write(TEMPLATE_JOB.format(initial_dir, n_cpus, memory, disk))

        step_params = [
            [[parcellation_file]],
            [
                [parcellation_file],
                marker_files,
                range(self.n_surrogate_maps),
                self.correlation_method,
                self.n_pca_covariates,
            ],
            [
                [parcellation_file],
                marker_files,
                self.correlation_method,
                self.alpha,
                self.n_pca_covariates,
            ],
        ]

        with open(submit_file, "a") as f:

            job_count = 0
            args = step_params[step - 1]
            for arg_tuple in product(*args):
                job_id = "_".join([f"step_{step}", parc_file_name])
                arguments = " ".join(
                    [f"./run_in_venv.sh step_{step}.py"]
                    + [str(x) for x in arg_tuple]
                )
                check_params = [x for x in arg_tuple]
                if not self._output_exists(step, *check_params):
                    f.write(
                        TEMPLATE_QUEUED_JOB.format(
                            arguments, logs=logs, job_id=job_id
                        )
                    )
                    job_count += 1

        return job_count

    def _prepare_submit_parc_dir(self, parcellation):

        _, parcellation_head = os.path.split(parcellation)
        parcellation_name = remove_nii_extensions(parcellation_head)
        submit_parc_dir = os.path.join(
            self.project_path, self.submit_files_dir, parcellation_name
        )

        if not os.path.isdir(submit_parc_dir):
            os.mkdir(submit_parc_dir)

        return submit_parc_dir

    def _prepare_run_dir(self, submit_parc_dir):

        runs = [x for x in os.listdir(submit_parc_dir) if "run" in x]
        run = len(runs) + 1
        run_dir = os.path.join(submit_parc_dir, f"run_{run}")
        logs_dir = os.path.join(run_dir, "logs")
        for directory in [run_dir, logs_dir]:
            os.mkdir(directory)

        return run_dir

    def _prepare_dag(self, run_dir):
        nimgen_dag = os.path.join(run_dir, "nimgen.dag")
        with open(nimgen_dag, "w") as f:
            f.write(TEMPLATE_DAG)

    def create(self):
        """Create the pipeline."""
        self.prepare_run_in_venv()
        for step in range(1, 4):
            self.prepare_step(step)

        self.prepare_submit_files()
