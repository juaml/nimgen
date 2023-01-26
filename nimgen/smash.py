"""Run permutation tests using brainsmash to remove spatial autocorrelation."""

# Authors: Yasir Demirtaş <tyasird@gmail.com>
#          Leonard Sasse <l.sasse@fz-juelich.de>
#          Kaustubh Patil <k.patil@fz-juelich.de>
# License: AGPL


import logging
import shutil
from itertools import combinations_with_replacement
from pathlib import Path

# import matplotlib.pyplot as plt
import neuromaps as nm
import nibabel as nib
import numpy as np
from neuromaps.nulls import burt2020
from neuromaps.parcellate import Parcellater
from scipy import ndimage
from scipy.spatial.distance import cdist  # , pdist

from .utils import remove_nii_extensions


logger = logging.getLogger(__name__)


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
    logger.info("Generating a fresh distance matrix.")
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


def cached_distance_matrix(parcellation_path, force_overwrite=False):
    """Load a chached distance matrix for a given parcellation.

    If a distance matrix is not yet cached, it will be generated.

    Parameters
    ----------
    parcellation_path : path
        Nifti file of parcellation for which to generate the distance matrix.
    force_overwrite : boolean
        Whether or not to force overwriting the cached distance matrix.

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

    if dist_mat_file.is_file() and not force_overwrite:
        logger.info(f"{dist_mat_file} already exists! Loading...")
        dist_mat = np.load(dist_mat_file)
    else:
        logger.info(f"Creating a new distance matrix at {dist_mat_file}")
        dist_mat = vox_dist(parcellation_path)
        np.save(dist_mat_file, dist_mat)

    return dist_mat


# def _plot_spatial_correlations(data, nulls, dist_mat, fname):

# data_array = np.array(data)[np.newaxis]
# n_parcels, n_perm = nulls.shape
# distance_data = pdist(data_array.T)
# dspat = dist_mat[np.triu_indices(n_parcels, k=1)]
# logger.info(f"spatial correlation {spearmanr(distance_data, dspat)}")

# spatcorr = np.zeros((n_perm, 1))
# for i in range(n_perm):
#     d = np.array(nulls[:, i])[np.newaxis]
#     ddata = pdist(d.T)
#     spatcorr[i], _ = spearmanr(distance_data, dspat)

# TODO: Add histograms later
# counts, bin_edges = np.histogram(spatcorr)
# fig = plt.figure()
# fig.hist(counts, bin_edges)
# plt.savefig(fname)


def cached_null_maps(
    parcellation_path,
    marker_path,
    distmat,
    n_perm,
    seed=None,
    force_overwrite=False,
):
    """Load a chached null maps file for a parcellation/marker combination.

    If the null maps are not yet cached, they will be generated.

    Parameters
    ----------
    parcellation_path : path
        Nifti file of parcellation for which to generate the distance matrix.
    marker_path : path
        Path to nifti file of marker of interest.
    distmat : path
        Path to npy file of distance matrix for this parcellation.
    n_perm : int
        How many null maps should be generated.
    seed : int
        Random seed for null map generation.
    force_overwrite : bool
        Whether to force overwriting the cached null maps with new ones.

    Returns
    -------
    null_maps : numpy.array
        Null maps.
    """
    # TODO: ADD SHAPE INFO TO THE null_maps docstring
    # get the path to the marker file
    marker_path = Path(marker_path)
    parcellation_path = Path(parcellation_path)
    marker_parent = marker_path.parent

    # get the name of the marker and parcellation
    marker_name = remove_nii_extensions(marker_path.stem)
    parc_name = remove_nii_extensions(parcellation_path.stem)

    # path and filename null maps
    null_maps_file_name = (
        f"marker-{marker_name}_perms-{n_perm}_seed-{seed}_desc-nullmaps"
    )
    null_maps_dir = marker_parent / "nullmaps" / parc_name

    null_maps_dir.mkdir(parents=True, exist_ok=True)

    null_maps_file = null_maps_dir / f"{null_maps_file_name}.npy"
    # null_maps_plot_file = null_maps_dir / f"{null_maps_file_name}.svg"

    if null_maps_file.is_file() and not force_overwrite:
        logger.info(f"{null_maps_file} already exists! Loading...")
        return np.load(null_maps_file)

    logger.info(f"Creating new null maps file at {null_maps_file}")
    masker = Parcellater(parcellation_path, space="MNI152").fit()
    data = masker.transform(marker_path, space="MNI152")[0]

    # call neuromaps function to generate null maps
    null_maps = burt2020(
        data=data,
        atlas="MNI152",
        density="2mm",  # TODO: Use correct resolution programmatically
        n_perm=n_perm,
        parcellation=parcellation_path,
        seed=seed,
        distmat=distmat,
    )
    logger.info(f"Caching newly generated null maps at {null_maps_file}")
    np.save(null_maps_file, null_maps)
    # _plot_spatial_correlations(
    #     data, null_maps, distmat, null_maps_plot_file
    # )
    shutil.copy(
        parcellation_path, null_maps_dir
    )  # TODO: Check if this is necessary and perhaps remove
    logger.info("Saving parcellation with nullmaps.")

    return null_maps
