import os
from pathlib import Path
import shutil


def remove_nii_extensions(nii_file):
    nii_name = Path(nii_file)

    while nii_name.suffix in {".nii", ".gz"}:
        nii_name = nii_name.with_suffix("")

    return nii_name


class Pipeline:
    def __init__(self, config_dict):
        self.project_path = os.path.abspath(config_dict["project_path"])

        self.pipeline_type = config_dict["pipeline"]
        self.submit_files_dir = config_dict["submit_files_dir"]
        self.pipeline_dir = config_dict["pipeline_dir"]
        self.parcellation_files = config_dict["parcellation_files"]
        self.marker_dir = config_dict["marker_dir"]
        self.marker_files = config_dict["marker_files"]
        self.output_dir = config_dict["output_dir"]
        self.parcellations_dir = "parcellations"
        self.config_dict = config_dict

        if not os.path.isdir(self.project_path):
            raise FileNotFoundError(f"{self.project_path} not found!")

        if not os.path.isdir(self.marker_dir):
            alternatively = os.path.join(self.project_path, self.marker_dir)
            if not os.path.isdir(alternatively):
                raise FileNotFoundError(f"{self.marker_dir} not found!")
            else:
                self.marker_dir = os.path.abspath(alternatively)

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

        for parcellation_file in self.parcellation_files:
            if not os.path.isfile(parcellation_file):
                FileNotFoundError(f"{parcellation_file} not found!")
            else:
                this_parc = remove_nii_extensions(parcellation_file)
                dir_this_parc = os.path.join(
                    self.parcellations_dir, this_parc
                )
                if not os.path.isdir(dir_this_parc):
                    print(f"Creating {dir_this_parc}")
                else:
                    print(
                        f"{dir_this_parc} already exists!"
                        f" Skipping creation of dir {dir}"
                    )
                head, tail = os.path.split(parcellation_file)
                parc_file_new = os.path.join(dir_this_parc, tail)
                if not os.path.isfile(parc_file_new):
                    print(
                        f"Copying from {parcellation_file} to {dir_this_parc}"
                    )
                    shutil.copyfile(parcellation_file, dir_this_parc)
                else:
                    print(f"{parc_file_new} already exists. Skipping copy.")

        for marker_file in self.marker_files:
            current_marker = os.path.join(self.marker_dir, marker_file)
            if not os.path.isfile(current_marker):
                raise FileNotFoundError(
                    "Marker files are interpreted relative to marker_dir."
                    f"({self.marker_dir})"
                )


"""
    def step_one(self):
        PARCEL_FILE = sys.argv[1]
        MARKER_FILE = sys.argv[2]
        N_SURROGATE_MAPS = sys.argv[3]
        project_path="/data/project/nimgen/nimgen_sex_diff"

        if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
            raise ValueError('Input file not found.')

        # create sample path specific for marker and atlas
        _, _, parcellation_path = create_sample_path(
            PARCEL_FILE, MARKER_FILE, project_path, True
        )

        # check surrogate maps for given atlas
        smaps_dir = os.path.join(parcellation_path, 'smaps', '*.nii')
        smaps = glob.glob(smaps_dir)

        # if number of surrogate maps are enough, do not create distance matrix
        if len(smaps) >= int(N_SURROGATE_MAPS):
            print(
                f'matrix file creation skipped.'
                f' In the surrogate map folder there are already {len(smaps)}'
                f' surrogate maps.'
            )
            pass
        else:
            coord_file, parcel_file = export_voxel_coordinates(
                PARCEL_FILE, MARKER_FILE)
            matrix_files = generate_distance_matrices(PARCEL_FILE, MARKER_FILE)

    def step_two():
        pass

    def step_three():
        pass

    def step_four():
        pass
"""
