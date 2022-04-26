import os
from black import out
import nilearn
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
from time import perf_counter
import tempfile
from .expressions import *
from .aggregation import *
from .utils import logger


# smashed nifti files
def smash_parcellation(
    parcellation_file, output_dir, chunck_size=1000, smap=10, n_jobs=1, tmp=None
):
    voxel_coord_file, voxel_parcel_file = _export_voxel_coordinates(
        parcellation_file, output_dir
    )
    surrogate_maps, elapsed_time = _generate_surrogate_maps(
        voxel_coord_file,
        voxel_parcel_file,
        output_dir,
        chunck_size=chunck_size,
        smap=smap,
        n_jobs=n_jobs,
        tmp=tmp,
    )
    smaps = np.load(os.path.join(output_dir, "smaps.npy"))
    smashed_atlaeses = []
    smashed_atlases_dir = os.path.join(output_dir, "smashed_atlases")
    if not os.path.isdir(smashed_atlases_dir):
        os.mkdir(smashed_atlases_dir)

    pc1 = perf_counter()
    for key, smap in enumerate(smaps):
        filename = os.path.join(smashed_atlases_dir, f"{key}_smashed.nii")
        _create_nifti(smap, parcellation_file, filename)
        smashed_atlaeses.append(filename)

    pc2 = perf_counter()
    elapsed_time.loc[len(elapsed_time)] = ["smaps_create_nifti", (pc2 - pc1) / 60]

    return smashed_atlaeses, elapsed_time


# voxel_coordinates, brain_map files
def _export_voxel_coordinates(parcellation_file, output_dir):
    voxel_coord_file = os.path.join(output_dir, "voxel_coordinates.txt")
    voxel_parcel_file = os.path.join(output_dir, "brain_map.txt")
    # load parcellation file
    niimg = nilearn.image.load_img(parcellation_file)
    # get image data as a numpy array
    data = niimg.get_fdata()
    # find voxels that are not zeroes
    idx = np.where(data > 0)
    # list of arrays to (voxels, 3) array
    ijk = np.vstack(idx).T
    parcels = [data[tuple(i)] for i in ijk]
    # xyz = nib.affines.apply_affine(niimg.affine, ijk)
    # get mm coords
    coords = nilearn.image.coord_transform(idx[0], idx[1], idx[2], niimg.affine)
    coords = np.vstack(coords).T

    np.savetxt(voxel_coord_file, coords, fmt="%s")
    np.savetxt(voxel_parcel_file, parcels, fmt="%s")
    return voxel_coord_file, voxel_parcel_file


# distance, index, smaps files
def _generate_surrogate_maps(
    voxel_coord,
    voxel_parcel,
    output_dir,
    chunck_size,
    smap,
    n_jobs,
    distance_matrix=None,
    index_matix=None,
    tmp=None,
):
    # tmpdir = tempfile.TemporaryDirectory()

    pc1 = perf_counter()
    # filenames = {}
    # if distance_matrix and index_matix:
    #     filenames["D"] = distance_matrix
    #     filenames["index"] = index_matix
    # else:
    #     filenames = volume(voxel_coord, tmpdir.name, chunk_size=chunck_size)
    tmpdir = output_dir
    if tmp:
        tmpdir = tmp
    logger.info(f"tmp_dir: {tmp}")
    print(f"tmp_dir: {tmp}")
    filenames = volume(voxel_coord, tmpdir, chunk_size=chunck_size)

    pc2 = perf_counter()
    logger.info(f"volume function finished in {(pc2 - pc1) / 60:0.0f} minutes")
    # These are three of the key parameters affecting the variogram fit
    kwargs = {"ns": 500, "knn": 1500, "pv": 70}
    gen = Sampled(
        x=voxel_parcel,
        D=filenames["D"],
        index=filenames["index"],
        resample=True,
        n_jobs=n_jobs,
        **kwargs,
    )
    surrogate_maps = gen(n=smap)
    # tmpdir.cleanup()
    os.remove(os.path.join(tmpdir, "distmat.npy"))
    os.remove(os.path.join(tmpdir, "index.npy"))
    pc3 = perf_counter()
    logger.info(f"Sampled function finished in {(pc3 - pc2) / 60:0.0f} minutes")
    np.save(os.path.join(output_dir, "smaps.npy"), surrogate_maps)
    df = pd.DataFrame(columns=["function", "elapsed"])
    df.loc[len(df)] = ["volume", (pc2 - pc1) / 60]
    df.loc[len(df)] = ["smaps_generate", (pc3 - pc2) / 60]

    return surrogate_maps, df


def _create_nifti(xyz_arr, ref_parcellation_file, output_filename):
    niimg = nilearn.image.load_img(ref_parcellation_file)
    data = niimg.get_fdata()
    idx = np.where(data > 0)
    ijk = np.vstack(idx)
    data[ijk[0], ijk[1], ijk[2]] = xyz_arr
    nii = nilearn.image.new_img_like(ref_parcellation_file, data)
    nii.to_filename(output_filename)
    return nii


def _empirical_pval(stat, stat0):
    """
    stat:A vector of calculated test statistics.
    stat0:A vector or matrix of simulated or data-resampled null test statistics.
    """
    check = np.sum(np.abs(stat) > np.abs(stat0), axis=0)
    pvalues = (check + 1) / (len(stat) + 1)
    return pvalues



def export_significance_genes(reference_data, smashed_data, output_dir, **kwargs):
    if os.path.isfile(smashed_data):
        smashed_data = np.load(smashed_data, allow_pickle=True)

    smashed_concat = pd.concat(smashed_data, axis=1).drop(columns=["pval"])
    reference_data.drop(columns=["pval"], inplace=True)
    stat = smashed_concat.T.values
    stat0 = reference_data.T.values
    empirical_pval = _empirical_pval(stat, stat0)
    df = pd.DataFrame(
        {"genes": reference_data.index, "pval": empirical_pval}
    ).set_index("genes")
    reject, corrected, *_ = multipletests(
        df["pval"], alpha=0.05, method="fdr_bh", is_sorted=False, returnsorted=False
    )
    df["fdr"] = corrected
    df.drop(columns=["pval"], inplace=True)
    np.savetxt( os.path.join(output_dir, 'genes.txt'), df.index, fmt='%s')
    df.to_csv(os.path.join(output_dir, 'genes.tsv'), sep='\t')
    return df[reject]
