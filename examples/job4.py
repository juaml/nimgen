#!/usr/bin/env python3
from nimgen import create_sample_path, get_corr_scores
import os
import sys

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]
PARTIAL_CORR = bool(sys.argv[3].split('=')[1])
PERFORM_PCA = bool(sys.argv[4].split('=')[1])
N_PCA_COMP = sys.argv[5].split('=')[1]

_, output_path, _ = create_sample_path(PARCEL_FILE, MARKER_FILE)

# for only original atlas perform PCA 
corr_scores, significant_genes, pca_components = get_corr_scores(
        PARCEL_FILE,
        MARKER_FILE,
        partial_correlation=PARTIAL_CORR,
        perform_pca=PERFORM_PCA,
        n_pca_comp=N_PCA_COMP,
        custom_covariates_df=None,
        is_surrogate=False,
)
