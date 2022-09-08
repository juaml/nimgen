#!/usr/bin/env python3
from nimgen import export_voxel_coordinates, generate_distance_matrices, create_sample_path
import sys
import os
import glob

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]
N_SURROGATE_MAPS = sys.argv[3]

if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
    raise ValueError('Input file not found.')

# create sample path specific for marker and atlas
_, _, parcellation_path = create_sample_path(PARCEL_FILE, MARKER_FILE, True)

# check surrogate maps for given atlas
smaps_dir = os.path.join(parcellation_path, 'smaps', '*.nii')
smaps = glob.glob(smaps_dir)

# if number of surrogate maps are enough, do not create distance matrix
if len(smaps) >= int(N_SURROGATE_MAPS):
    print(f'matrix file creation skipped. In the surrogate map folder there are allready {len(smaps)} surrogate maps.')
    pass
else:
    coord_file, parcel_file = export_voxel_coordinates(PARCEL_FILE, MARKER_FILE)
    matrix_files = generate_distance_matrices(PARCEL_FILE, MARKER_FILE)

