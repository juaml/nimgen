#!/usr/bin/env python3
from nimgen import generate_surrogate_map, get_corr_scores, create_sample_path
import sys
import os

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]
SMAP_ID = sys.argv[3]

if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
    raise ValueError('Input file not found.')

# generate surrogate map for given atlas
surrogate_map = generate_surrogate_map(PARCEL_FILE, MARKER_FILE, SMAP_ID)

# perform correlation analysis for given surrogate map and marker
get_corr_scores(
    parcellation_file=surrogate_map,
    marker_file=MARKER_FILE,
)
