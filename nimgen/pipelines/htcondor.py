import os
from pathlib import Path
import shutil
from ._file_strings import (
    STEP_ONE_FSTRING,
    STEP_TWO_FSTRING,
    STEP_THREE_FSTRING,
    STEP_FOUR_FSTRING
)

def remove_nii_extensions(nii_file):
    nii_name = Path(nii_file)

    while nii_name.suffix in {".nii", ".gz"}:
        nii_name = nii_name.with_suffix("")

    return nii_name


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

    def prepare_step_one(self):
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_one_file = os.path.join(pipeline_dir, "step_one.py")
        if not os.path.isfile(step_one_file):
            with open(step_one_file, "w") as f:
                f.write(
                    STEP_ONE_FSTRING.format(
                        self.config_dict["n_surrogate_maps"],
                        self.project_path,
                        os.path.join(self.project_path, self.parcellations_dir)
                    )
                )

    def prepare_step_two(self):
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_two_file = os.path.join(pipeline_dir, "step_two.py")
        if not os.path.isfile(step_two_file):
            with open(step_two_file, "w") as f:
                f.write(
                    STEP_TWO_FSTRING.format(
                        self.project_path
                    )
                )

    def prepare_step_three(self):
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_three_file = os.path.join(pipeline_dir, "step_three.py")
        if not os.path.isfile(step_three_file):
            with open(step_three_file, "w") as f:
                f.write(
                    STEP_THREE_FSTRING.format(
                        self.project_path,
                        self.marker_dir,
                        os.path.join(self.project_path, self.output_dir),
                        self.r_path
                    )
                )

    def prepare_step_four(self):
        pipeline_dir = os.path.join(self.project_path, self.pipeline_dir)
        step_four_file = os.path.join(pipeline_dir, "step_four.py")
        if not os.path.isfile(step_four_file):
            with open(step_four_file, "w") as f:
                f.write(
                    STEP_FOUR_FSTRING.format(
                    )
                )


    def prepare_submit_files(self):
        pass

    def create(self):
        self.prepare_run_in_venv()
        self.prepare_step_one()
        self.prepare_step_two()
        self.prepare_step_three()
        self.prepare_step_four()
        self.prepare_submit_files()