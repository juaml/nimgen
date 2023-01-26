"""Define HTCondor pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
from itertools import product
from pathlib import Path

from ._htcondor_python_strings import RUN_CONDA, RUN_IN_VENV
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

        self.submit_files_dir = self.jobs_dir / submit_files_dir
        self.logs_dir = self.submit_files_dir / "logs"

    def create(self):
        """Create a HTCondor-specific pipeline."""
        super().create()
        self.submit_files_dir.mkdir()
        self.logs_dir.mkdir()
        self.prepare_run_in_venv()
        self.prepare_submit_files()

    def prepare_run_in_venv(self):
        """Prepare a bash file to source a venv or conda env."""

        run_in_venv = self.submit_files_dir / "run_in_venv.sh"

        if self.conda:
            with open(run_in_venv, "w") as f:
                f.write(RUN_CONDA.format(conda_env=self.conda_env))
        else:
            with open(run_in_venv, "w") as f:
                f.write(RUN_IN_VENV.format(venv=self.path_to_venv))

        os.system(f"chmod +x {run_in_venv}")

    def prepare_submit_files(self):
        """Prepare submit files to submit the pipeline to HTCondor as a DAG."""

        for step in range(1, 5):
            job_count = self._prepare_submit_file(step)
            print(f"{job_count} jobs registered for step {step}.")

        self._prepare_dag()

    def _prepare_submit_file(self, step):

        parent_dir = Path("..") / ".."
        initial_dir = parent_dir / self.submit_files_dir.as_posix()
        n_cpus = self.pipeline[f"step_{step}"]["CPU"]
        memory = self.pipeline[f"step_{step}"]["MEMORY"]
        disk = self.pipeline[f"step_{step}"]["DISK"]

        submit_file = self.submit_files_dir / f"step_{step}.submit"
        with open(submit_file, "w") as f:
            f.write(TEMPLATE_JOB.format(initial_dir, n_cpus, memory, disk))

        step_args = {
            "step_1": {
                "parcellation_file": [
                    parent_dir / x for x in self.all_parcellations
                ],
            },
        }

        job_count = 0

        logs_dir = (Path("..") / ".." / self.logs_dir).as_posix()

        for marker in self.markers:

            step_args["step_2"] = {
                "parcellation_file": [
                    parent_dir / x for x in marker["parcellation"]
                ],
                "marker_file": [parent_dir / marker["path"]],
                "n_perm": [self.n_surrogate_maps],
                "seed": [self.seed],
            }
            step_args["step_3"] = {
                "parcellation_file": [
                    parent_dir / x for x in marker["parcellation"]
                ],
                "marker_file": [parent_dir / marker["path"]],
                "n_perm": [self.n_surrogate_maps],
                "smap_id": range(self.n_surrogate_maps),
                "seed": [self.seed],
                "allen_data_dir": [parent_dir / self.ahba_dir],
                "correlation_method": self.correlation_method,
            }
            step_args["step_4"] = {
                "parcellation_file": [
                    parent_dir / x for x in marker["parcellation"]
                ],
                "marker_file": [parent_dir / marker["path"]],
                "n_perm": [self.n_surrogate_maps],
                "seed": [self.seed],
                "allen_data_dir": [parent_dir / self.ahba_dir],
                "correlation_method": self.correlation_method,
                "alpha": self.alpha,
            }

            with open(submit_file, "a") as f:

                args = step_args[f"step_{step}"]
                arg_vals = [arg_val for _, arg_val in args.items()]

                for arg_tuple in product(*arg_vals):
                    arguments = " ".join(
                        [f"nimgen -v {self.verbosity} step_{step}"]
                        + [str(x) for x in arg_tuple]
                    )

                    job_id = f"step-{step}_job-{job_count}_"
                    f.write(
                        TEMPLATE_QUEUED_JOB.format(
                            arguments, logs=logs_dir, job_id=job_id
                        )
                    )
                    job_count += 1
        return job_count

    def _prepare_dag(self):
        nimgen_dag = self.submit_files_dir / "nimgen.dag"
        with open(nimgen_dag, "w") as f:
            f.write(TEMPLATE_DAG)
