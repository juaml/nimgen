"""Run permutation tests using brainsmash to remove spatial autocorrelation."""

# Authors: Yasir Demirtaş <tyasird@gmail.com>
#          Leonard Sasse <l.sasse@fz-juelich.de>
#          Kaustubh Patil <k.patil@fz-juelich.de>
# License: AGPL


from itertools import combinations_with_replacement
from pathlib import Path

import neuromaps as nm
import nibabel as nib
import numpy as np
from scipy import ndimage
from scipy.spatial.distance import cdist

from .utils import logger, remove_nii_extensions


def _vox_dist_original(parcellation):
    """Generate a distance matrix for a given parcellation.

    Parameters
    ----------
    parcellation : niimg-like
        Nifti file of parcellation for which to generate the distance matrix.

    Returns
    -------
    dist : numpy.array
        Distance matrix with shape (N_ROI x N_ROI).

    This function is copied from:
    https://netneurolab.github.io/neuromaps/
    _modules/neuromaps/nulls/nulls.html#burt2020
    """

    darr = nm.images.load_data(parcellation)
    affine = nm.images.load_nifti(parcellation).affine
    labels = np.trim_zeros(np.unique(darr))
    mask = np.logical_not(np.logical_or(np.isclose(darr, 0), np.isnan(darr)))
    xyz = nib.affines.apply_affine(affine, np.column_stack(np.where(mask)))

    parcellation = darr[mask]
    row_dist = np.zeros((len(xyz), len(labels)), dtype="float32")
    dist = np.zeros((len(labels), len(labels)), dtype="float32")
    for n, row in enumerate(xyz):
        xyz_dist = cdist(row[None], xyz).astype("float32")
        row_dist[n] = ndimage.mean(xyz_dist, parcellation, labels)
    for n in range(len(labels)):
        dist[n] = ndimage.mean(row_dist[:, n], parcellation, labels)

    return dist


def vox_dist(parcellation):
    """Generate a distance matrix for a given parcellation.

    Parameters
    ----------
    parcellation : niimg-like
        Nifti file of parcellation for which to generate the distance matrix.

    Returns
    -------
    dist : numpy.array
        Distance matrix with shape (N_ROI x N_ROI).

    This function is adapted from code available in:
    https://netneurolab.github.io/neuromaps/
    _modules/neuromaps/nulls/nulls.html#burt2020
    """
    logger.info("Loading parcellation to generate distance matrix...")
    # faster implementation by avoiding extra computation
    darr = nm.images.load_data(parcellation)
    affine = nm.images.load_nifti(parcellation).affine
    labels = np.trim_zeros(np.unique(darr))
    mask = np.logical_not(np.logical_or(np.isclose(darr, 0), np.isnan(darr)))
    logger.info("Applying affine to coordinates of interest...")
    xyz = nib.affines.apply_affine(affine, np.column_stack(np.where(mask)))

    parcellation = darr[mask]
    dist = np.zeros((len(labels), len(labels)), dtype="float32")
    logger.info("Generating distance matrix, this may take a while...")
    for i, j in combinations_with_replacement(range(len(labels)), 2):
        l1, l2 = labels[i], labels[j]
        dist[i, j] = (
            cdist(xyz[parcellation == l1, :], xyz[parcellation == l2, :])
            .mean()
            .astype("float32")
        )
        dist[j, i] = dist[i, j]

    return dist


def cached_distance_matrix(parcellation_path):
    """Load a chached distance matrix for a given parcellation.

    If a distance matrix is not yet cached, it will be generated.

    Parameters
    ----------
    parcellation : niimg-like
        Nifti file of parcellation for which to generate the distance matrix.

    Returns
    -------
    dist : numpy.array
        Distance matrix with shape (N_ROI x N_ROI).
    """
    # get the path to the parcellation
    parcellation_path = Path(parcellation_path)
    parent = parcellation_path.parent

    # get the name of the parcellation
    parc_name = remove_nii_extensions(parcellation_path.stem)

    # path and filename distance matrix
    dist_mat_file = parent / f"{parc_name}_dist_mat.npy"

    if dist_mat_file.is_file():
        logger.info(f"{dist_mat_file} already exists! Loading...")
        dist_mat = np.load(dist_mat_file)
    else:
        logger.info(f"{dist_mat_file} does not exist yet! Generating...")
        dist_mat = vox_dist(parcellation_path)
        logger.info("Caching newly generated distance matrix...")
        np.save(dist_mat_file, dist_mat)

    return dist_mat
