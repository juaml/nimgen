import os
import shutil
from ._htcondor_file_strings import (
    STEP_ONE_FSTRING,
    STEP_TWO_FSTRING,
    STEP_THREE_FSTRING,
)
from ..utils import remove_nii_extensions


class HTCondor:
    """ Object to create and run a nimgen pipeline on a HTCondor cluster.

    """

    def __init__(self, config_dict):
        """ Initialise HTCondor pipeline.

        Parameters
        -----------
        config_dict : dict
            dictionary with valid pipeline configurations.

        """
        self.project_path = os.path.abspath(config_dict["project_path"])

        self.pipeline_type = config_dict["pipeline"]
        self.submit_files_dir = config_dict["submit_files_dir"]
        self.pipeline_dir = config_dict["pipeline_dir"]
        self.marker_dir = config_dict["marker_dir"]
        self.output_dir = config_dict["output_dir"]
        self.parcellations_dir = "parcellations"
        self.config_dict = config_dict
        self.parcellation_marker_dict = {}
        self.r_path = os.path.abspath(config_dict["r_path"])

        if not os.path.isdir(self.project_path):
            raise FileNotFoundError(f"{self.project_path} not found!")

        if not os.path.isdir(self.marker_dir):
            alternatively = os.path.join(self.project_path, self.marker_dir)
            if not os.path.isdir(alternatively):
                raise FileNotFoundError(f"{self.marker_dir} not found!")
            else:
                self.marker_dir = os.path.abspath(alternatively)
        else:
            self.marker_dir = os.path.abspath(self.marker_dir)

        for dir in [
            self.submit_files_dir, self.output_dir,
            self.pipeline_dir, self.parcellations_dir
        ]:
            dir = os.path.join(self.project_path, dir)
            if not os.path.isdir(dir):
                print(f"Creating {dir}")
                os.mkdir(dir)
            else:
                print(f"{dir} already exists! Skipping creation of dir {dir}")

        for parcellation_file in config_dict["parcellation_files"].keys():
            if not os.path.isfile(parcellation_file):
                FileNotFoundError(f"{parcellation_file} not found!")
            else:
                _, this_parc = os.path.split(
                    remove_nii_extensions(parcellation_file)
                )
                self.parcellation_marker_dict[this_parc] = config_dict[
                    "parcellation_files"
                ][parcellation_file]
                dir_this_parc = os.path.join(
                    self.parcellations_dir, this_parc
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

                _, tail = os.path.split(parcellation_file)
                parc_file_new = os.path.join(dir_this_parc, tail)
                if not os.path.isfile(parc_file_new):
                    print(
                        f"Copying from {parcellation_file} to {dir_this_parc}"
                    )
                    shutil.copyfile(parcellation_file, parc_file_new)
                else:
                    print(f"{parc_file_new} already exists. Skipping copy.")

        for _, markers in self.parcellation_marker_dict.items():
            for marker_file in markers:
                current_marker = os.path.join(self.marker_dir, marker_file)
                if not os.path.isfile(current_marker):
                    raise FileNotFoundError(
                        "Marker files are interpreted relative to marker_dir."
                        f"({self.marker_dir})"
                    )

    def prepare_run_in_venv(self):
        pass

    def prepare_step(self, step):
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
        pass

    def create(self):
        self.prepare_run_in_venv()
        for step in range(1, 4):
            self.prepare_step(step)

        self.prepare_submit_files()
