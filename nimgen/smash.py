import os
from time import perf_counter

import numpy as np
from nilearn import image

from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled

from .utils import logger


def export_voxel_coordinates(parcellation_file, outpath):
    """
    Extracts XYZ voxel coordinates and parcel numbers from every voxel within
    an ROI image consisting of ones and zeros based on the parcellation file.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    outpath : str or os.PathLike
        directory in which to save voxel coordinate files

    Returns
    -------
    voxel_coord_file : str or os.PathLike
        voxel_coordinate file name `voxel_coordinates.txt`.
    voxel_parcel_file : str or os.PathLike
        voxel parcel file `brain_map.txt`.
    """

    logger.info(f"Trying to export voxel coordinates..")
    voxel_coord_file = os.path.join(outpath, "voxel_coordinates.txt")
    voxel_parcel_file = os.path.join(outpath, "brain_map.txt")

    # check txt files, if exists return filenames.
    if os.path.isfile(voxel_coord_file) and os.path.isfile(voxel_parcel_file):
        return voxel_coord_file, voxel_parcel_file

    # load parcellation file
    niimg = image.load_img(parcellation_file)
    # get image data as a numpy array
    data = niimg.get_fdata()
    # find voxels that are not zeroes
    idx = np.where(data > 0)
    # list of arrays to (voxels, 3) array
    ijk = np.vstack(idx).T
    parcels = [data[tuple(i)] for i in ijk]
    # xyz = nib.affines.apply_affine(niimg.affine, ijk)
    # get mm coords
    coords = image.coord_transform(
        idx[0], idx[1], idx[2], niimg.affine)
    coords = np.vstack(coords).T
    # save files as a txt file
    np.savetxt(voxel_coord_file, coords, fmt="%s")
    np.savetxt(voxel_parcel_file, parcels, fmt="%s")
    return voxel_coord_file, voxel_parcel_file


def generate_distance_matrices(
    path,
    chunk_size=1000,
):
    """
    Generates distance matrices for BrainSMASH based on the XYZ voxel
    coordinates.

    Parameters
    ----------
    path : str or os.PathLike
        path to directory, in which the distance matrix should be saved
    chunck_size : int, default 1000
        The number of voxels to process per chunk. For N voxels, this will
        impose a memory burden of N*`chunk_size` per iteration (in contrast to
        a memory burden of N*N for a single iteration, in the absence of
        chunking).

    Returns
    -------
    dict
        Keys are 'D' and 'index'; values are absolute paths to the
        corresponding files on disk. These files are used as inputs to
        `brainsmash.mapgen.sampled.Sampled`.
    """
    voxel_coordinate_file = os.path.join(path, 'voxel_coordinates.txt')
    print(f"Trying to generate distance matrices..")
    matrix_files = {
        'D': os.path.join(path, 'distmat.npy'),
        'index': os.path.join(path, 'index.npy')
    }

    # check distance_matrix_files, if there isn't any generate.
    if os.path.isfile(
        matrix_files["D"]
    ) and os.path.isfile(matrix_files["index"]):
        logger.info(f"Distance matrix files already exist.")
        return matrix_files

    pc1 = perf_counter()
    filenames = volume(voxel_coordinate_file, path, chunk_size=chunk_size)
    pc2 = perf_counter()
    logger.info(f"generate_distance_matrices: {(pc2 - pc1) / 60:0.0f} minutes")

    return filenames


def generate_surrogate_map(
    parcellation_file,
    smap_id,
    outpath,
    voxel_parcel_file,
    matrix_files
):
    """
    Randomly generates surrogate maps with matched spatial autocorrelation
    based on the parcellation file
    and matrix files.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    smap_id : str
        ID of surrogate maps to randomly generate.
    outpath : str or os.PathLike
        directory under which smaps will be saved
    voxel_parcel_file : str or os.PathLike
    matrix_files : dict
        dict with keys "D" and "index" for brainsmash distance matrix files.

    Returns
    -------
    smap_file : str or os.PathLike
        Surrogate brain map filename.
    """

    logger.info(f"Trying to generate surrogate map..")
    smaps_dir = os.path.join(outpath, "smaps")
    smap_file = os.path.join(smaps_dir, f"{smap_id}_smap.nii")

    if os.path.isfile(smap_file):
        logger.log(
            f'{smap_file} already exists!'
        )
        return smap_file

    pc1 = perf_counter()
    gen = Sampled(
        x=voxel_parcel_file,
        D=matrix_files['D'],
        index=matrix_files["index"],
        resample=True,
        n_jobs=1,
    )
    generated_smap = gen(n=1)
    pc2 = perf_counter()
    logger.info(f"generate_surrogate_maps: {(pc2 - pc1) / 60:0.0f} minutes")
    _create_nifti(generated_smap, parcellation_file, smap_file)
    return smap_file


def _create_nifti(xyz_arr, ref_parcellation_file, output_filename):
    """
    Creates nifti file based on the XYZ coordinates and reference
    parcellation file.

    Parameters
    ----------
    xyz_arr : list
        XYZ coordinates of brain map.
    ref_parcellation_file : str or os.PathLike
        Filename of reference parcellation file.
    output_filename : str or os.PathLike
        New filename for created nifti file.

    Returns
    -------
    nii : list
        Nifti file.
    """
    niimg = image.load_img(ref_parcellation_file)
    data = niimg.get_fdata()
    idx = np.where(data > 0)
    ijk = np.vstack(idx)
    data[ijk[0], ijk[1], ijk[2]] = xyz_arr
    nii = image.new_img_like(ref_parcellation_file, data)
    nii.to_filename(output_filename)
    return nii
