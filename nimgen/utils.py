# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
# License: AGPL
import logging
import json
import os

from pathlib import Path

import numpy as np
import pandas as pd

from nilearn import image

logger = logging.getLogger(name='nimgen')
logger.setLevel(logging.INFO)
# console = logging.StreamHandler()
# logger.addHandler(console)
logger.propagate = False


def save_as_json(data, file):
    with open('data.txt', 'w'):
        json.dump(data, file, sort_keys=True, indent=4,
                  ensure_ascii=False)


def raise_error(error):
    logger.error(error)
    raise


def remove_nii_extensions(nii_file):
    nii_name = Path(nii_file)

    while nii_name.suffix in {".nii", ".gz"}:
        nii_name = nii_name.with_suffix("")

    return nii_name


def read_csv_tsv(path):
    _, ext = os.path.splitext(path)
    extensions = {
        ".csv": ",",
        ".tsv": "\t"
    }

    return pd.read_csv(
        path, index_col=0, sep=extensions[ext]
    )


def covariates_to_nifti(parcellation, covariates_df):
    """
    Creates nifti images for given PCA covariates.

    Parameters
    ----------
    parccellation : niimg-like object
        A parcellation image in MNI space, where each parcel is identified by a
        unique integer ID.
    covariates_df : dict
        PCA covariates. Each key represents different covariate.

    Returns
    -------
    covariate_niftis : dict
        A dictionary contains niimg-like object for each covariate.
    """
    parcellation_array = np.array(parcellation.dataobj)
    covariate_niftis = {}
    for covariate_label, covariate in covariates_df.items():
        null_mat = np.zeros(parcellation_array.shape)
        null_mat[parcellation_array == 0] = -20000
        marker_min = np.min(covariate)
        marker_max = np.max(covariate)
        for label, value in covariate.iteritems():
            null_mat[parcellation_array == label + 1] = value

        pc_nii = image.new_img_like(
            parcellation, null_mat
        )
        pc_nii.header["cal_min"] = marker_min
        pc_nii.header["cal_max"] = marker_max
        covariate_niftis[covariate_label] = pc_nii

    return covariate_niftis


def _read_sign_genes(sign_genes):
    if isinstance(sign_genes, pd.DataFrame):
        sign_genes = sign_genes.index
    elif not isinstance(sign_genes, pd.DataFrame):
        _, ext = os.path.splitext(sign_genes)
        if ext in [".csv", ".tsv"]:
            extensions = {".csv": ",", ".tsv": "\t"}
            sign_genes = pd.read_csv(
                sign_genes,
                header=None,
                index_col=0,
                dtype=str,
                sep=extensions[ext]
            ).index
        elif ext in [".txt"]:
            sign_genes = list(np.loadtxt(sign_genes, dtype=str))
        else:
            raise ValueError(
                "'sign_genes' should be a pd.DataFrame,"
                " .csv/.tsv or .txt file!"
            )

    return list(sign_genes)


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
