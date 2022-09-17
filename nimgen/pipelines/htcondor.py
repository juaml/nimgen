"""Define HTCondor pipeline."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
from ._htcondor_python_strings import (
    STEP_ONE_FSTRING,
    STEP_TWO_FSTRING,
    STEP_THREE_FSTRING,
)
from .base import Pipeline


class HTCondor(Pipeline):
    """Object to create and run a nimgen pipeline on a HTCondor cluster.

    Parameters
    -----------
    config_dict : dict
        dictionary with valid pipeline configurations.
    """

    def __init__(self, config_dict):
        """Initialise HTCondor pipeline."""
        super().__init__(config_dict)
        self.submit_files_dir = config_dict["submit_files_dir"]
        directory = os.path.join(self.project_path, self.submit_files_dir)
        if not os.path.isdir(directory):
            print(f"Creating {directory}")
            os.mkdir(directory)
        else:
            print(
                f"{directory} already exists! Skipping creation of dir {dir}"
            )

    def prepare_run_in_venv(self):
        """Prepare a bash file to source a python venv for the pipeline."""
        pass

    def prepare_step(self, step):
        """Prepare a python3 file for a given step in the pipeline."""
        _, name_marker_dir = os.path.split(self.marker_dir)
        name_output_dir = self.output_dir
        allen_data_dir = os.path.join(self.project_path, "allen_data_dir")

        step_args = [
            (
                STEP_ONE_FSTRING,
                ["placeholder"]
            ),
            (
                STEP_TWO_FSTRING,
                [name_marker_dir, name_output_dir, allen_data_dir]
            ),
            (
                STEP_THREE_FSTRING,
                [name_marker_dir, name_output_dir, allen_data_dir, self.r_path]
            )
        ]
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_file = os.path.join(pipeline_dir, f"step_{step}.py")
        if not os.path.isfile(step_file):
            with open(step_file, "w") as f:
                fstring, args = step_args[step - 1]
                f.write(fstring.format(*args))

    def prepare_submit_files(self):
        """Prepare submit files to submit the pipeline to HTCondor as a DAG."""
        pass

    def create(self):
        """Create the pipeline."""
        self.prepare_run_in_venv()
        for step in range(1, 4):
            self.prepare_step(step)

        self.prepare_submit_files()
