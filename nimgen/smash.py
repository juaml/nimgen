import os
from os.path import isfile
from pathlib import Path
import subprocess
import glob
from time import perf_counter

import nilearn
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.sampled import Sampled
from nilearn import image, masking

from functools import partial
from scipy.stats.mstats import winsorize
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd

from .expressions import get_gene_expression
from .utils import logger


def r_script_file():
    """
    Returns webgestalt.r abspath. It is generally inside in the r_script folder
    in the package.
    """
    r_file = os.path.dirname(os.path.abspath(
        __file__)) + '/../r_script/webgestalt.r'
    return r_file


def create_sample_path(
    parcellation_file, marker_file, project_path, create=False
):
    """
    Creates sample_output_path using parcellation and marker filenames.
    If create is False, function will return only abspath without creating.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT. Extracts
        (and returns) measures of region-wise gray matter density (GMD).
        Aggregation method is mean.
    project_path : str or os.PathLike
        path to project directory, in which the input/output folders should be
        available
    create : bool
        If create is False, function will return only abspath without creating
        directories.

    Returns
    -------
    project_path : str or os.PathLike
        Project path
    output_path : str or os.PathLike
        Output path for given spesific parcellation and marker file.
    parcellation_path : str or os.PathLike
        Input path of parcellation file. It is generally using for
        save parcel spesific files i.e.
        distance matrix files, surrogate maps etc.
    """
    project_path = os.path.abspath(project_path)
    m = Path(marker_file.split('input/')[1])
    p = Path(parcellation_file.split('input/')[1])
    output_path = (
        f'{project_path}/output/{m.parts[1]}/{p.parts[1].split(".")[0]}/'
        f'{os.path.join(*m.parts[2:-1])}/{m.parts[-1].split(".")[0]}'
    )

    parcellation_path = os.path.dirname(parcellation_file)

    if create:
        input_smaps_dir = os.path.join(parcellation_path, 'smaps')
        corr_scores_dir = os.path.join(output_path, 'smap_corr_scores')
        pca_covariates_dir = os.path.join(output_path, 'pca_covariates')
        if not os.path.isdir(pca_covariates_dir):
            os.makedirs(pca_covariates_dir)
        if not os.path.isdir(corr_scores_dir):
            os.makedirs(corr_scores_dir)
        if not os.path.isdir(input_smaps_dir):
            os.makedirs(input_smaps_dir)
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        logger.info(f"Sample path created: {output_path}")

    return project_path, output_path, parcellation_path


def export_voxel_coordinates(parcellation_file, marker_file):
    """
    Extracts XYZ voxel coordinates and parcel numbers from every voxel within
    an ROI image consisting of ones and zeros based on the parcellation file.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT.
        Extracts (and returns) measures of region-wise gray matter density
        (GMD). Aggregation method is mean.

    Returns
    -------
    voxel_coord_file : str or os.PathLike
        voxel_coordinate file name (voxel_coordinates.txt).
    voxel_parcel_file : str or os.PathLike
        voxel parcel file name (brain_map.txt).
    """

    _, _, parcellation_path = create_sample_path(
        parcellation_file, marker_file, True)

    logger.info(f"Trying to export voxel coordinates..")
    voxel_coord_file = os.path.join(parcellation_path, "voxel_coordinates.txt")
    voxel_parcel_file = os.path.join(parcellation_path, "brain_map.txt")

    # check txt files, if exists return filenames.
    if isfile(voxel_coord_file) and isfile(voxel_parcel_file):
        print(f"voxel_coord_file & voxel_parcel_file allready exists.")
        return voxel_coord_file, voxel_parcel_file

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
    coords = nilearn.image.coord_transform(
        idx[0], idx[1], idx[2], niimg.affine)
    coords = np.vstack(coords).T
    # save files as a txt file
    np.savetxt(voxel_coord_file, coords, fmt="%s")
    np.savetxt(voxel_parcel_file, parcels, fmt="%s")
    return voxel_coord_file, voxel_parcel_file


def generate_distance_matrices(
    parcellation_file,
    marker_file,
    project_path,
    chunck_size=1000,
):
    """
    Generates distance matrices for BrainSMASH based on the XYZ voxel
    coordinates.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT.
        Extracts (and returns) measures of region-wise gray matter density
        (GMD). Aggregation method is mean.
    project_path : str or os.PathLike
        path to project directory, in which the input/output folders should be
        available.
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

    _, _, parcellation_path = create_sample_path(
        parcellation_file, marker_file, project_path
    )

    voxel_coordinate_file = os.path.join(
        parcellation_path, 'voxel_coordinates.txt')

    print(f"Trying to generate distance matrices..")
    matrix_files = {
        'D': os.path.join(parcellation_path, 'distmat.npy'),
        'index': os.path.join(parcellation_path, 'index.npy')
    }

    # check distance_matrix_files, if there isn't any generate.
    if isfile(matrix_files["D"]) and isfile(matrix_files["index"]):
        print(f"distance_matrix_files allready exists.")
        return matrix_files

    pc1 = perf_counter()
    filenames = volume(
        voxel_coordinate_file,
        parcellation_path,
        chunk_size=chunck_size)
    pc2 = perf_counter()
    print(f"generate_distance_matrices: {(pc2 - pc1) / 60:0.0f} minutes")

    return filenames


def generate_surrogate_map(
    parcellation_file, marker_file, smap_id, project_path
):
    """
    Randomly generates surrogate maps with matched spatial autocorrelation
    based on the parcellation file
    and matrix files. Matrix files should be in the same directory with
    coordinate files.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT.
        Extracts (and returns) measures of region-wise gray matter density
        (GMD). Aggregation method is mean.
    smap_id : int
        ID of surrogate maps to randomly generate.
    project_path : str or os.PathLike
        path to project directory, in which the input/output folders should be
        available

    Returns
    -------
    smap_file : str or os.PathLike
        Surrogate brain map filename.
    """

    _, _, parcellation_path = create_sample_path(
        parcellation_file, marker_file, project_path, create=True)

    logger.info(f"Trying to generate surrogate map..")
    smaps_dir = os.path.join(parcellation_path, "smaps")
    smap_file = os.path.join(smaps_dir, f"{smap_id}_smap.nii")

    if os.path.isfile(smap_file):
        print(f'{smap_file} allready exits, continue to correlation analysis')
        return smap_file

    voxel_parcel_file = os.path.join(parcellation_path, "brain_map.txt")
    matrix_files = {
        'D': os.path.join(parcellation_path, 'distmat.npy'),
        'index': os.path.join(parcellation_path, 'index.npy')
    }
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
    Calculates empirical pvalue based on the observed (surrogate maps) and
    expected (reference map) correlation scores.

    Parameters
    ----------
    stat: numpy.array or list
        A vector of calculated test statistics.
    stat0: numpy.array or list
        A vector or matrix of simulated or data-resampled null test statistics.

    Returns
    -------
    pvalues : numpy.array
        Calculated empirical pvalues.
    """

    logger.info(f"Empirical pvalue calculation..")
    check = np.sum(np.abs(stat) > np.abs(stat0), axis=0)
    pvalues = (check + 1) / (len(stat) + 1)
    return pvalues


def get_corr_scores(
    parcellation_file,
    marker_file,
    project_path,
    output_path,
    correlation_method="spearman",
    alpha=0.05,
    partial_correlation=False,
    perform_pca=False,
    n_pca_comp=None,
    custom_covariates_df=None,
    is_surrogate=True,
    allen_data_dir='allen_data',
    aggregation_methods='mean',
    save_scores=True
):
    """
    Fetchs gene expression data from Allen Human Brain Atlas based on the
    parcellation file, conducts correlation analysis between gene expression
    and markers. Exports r_score and p_value of each gene to csv file.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT. Extracts
        (and returns) measures of region-wise gray matter density (GMD).
        Aggregation method is mean.
    project_path : str or os.PathLike
        path to project directory, in which the input/output folders should be
        available
    correlation_method : str
        'spearman' or 'pearson', kind of correlation to use for mass univariate
        analysis
    alpha : float
        alpha threshold for correlation analysis
    partial_correlation : bool, default = False
        Applies partial correlation for PCA covariates. Works with PCA.
        https://pingouin-stats.org/generated/pingouin.partial_corr.html
    perform_pca : bool, default = False,
        Applies Principal component analysis (PCA) using scikit-learn.
    n_pca_comp : int, default = 0
        Number of components to keep.
    custom_covariates_df : dict, default = None
        If PCA is disabled and partial_correlation is enabled,
        custom_covariates should be defined.
    is_surrogate : bool, default = True
        If given parcellation file is surrogate map, this option should be
        True.
    allen_data_dir : str or os.PathLike
        Gene expression data directory for Allen Human Brain Atlas
    aggregation_methods : list, default = mean
        List with aggregation methods the corresponding array with
        the calculated GMD based on the provided atlas.
    save_scores : bool, default = True
        If True, correlation scores will be saved as a csv file.

    Returns
    -------
    dict : corr_scores, significant_genes
        Correlation score of all genes
    """

    if aggregation_methods in ["mean"]:
        aggregation_methods = ["mean"]

    parcellation_path, parcellation_file_tail = os.path.split(
        parcellation_file
    )
    parcellation_name, _ = os.path.splitext(parcellation_file_tail)

    allen_dir = os.path.join(project_path, allen_data_dir)

    marker_aggregated, agg_func_params = _aggregate_marker(
        parcellation_file, marker_file, aggregation=aggregation_methods
    )

    perform_pca = False
    if n_pca_comp is None:
        perform_pca = False
        pca_dict = None
    else:
        assert isinstance(n_pca_comp, int), "n_pca_comp should be int or None!"
        assert n_pca_comp > 0, "n_pca_comp should be None or greater than 0!"
        perform_pca = True
        pca_dict = {"n_components": n_pca_comp}

    # corr analysis
    corr_scores, significant_genes, pca_components = get_gene_expression(
        marker_aggregated['mean'],
        parcellation_file,
        allen_data_dir=allen_dir,
        save_expressions=True,
        correlation_method=correlation_method,
        alpha=alpha,
        partial_correlation=partial_correlation,
        perform_pca=perform_pca,
        pca_dict=pca_dict,
        custom_covariates_df=custom_covariates_df,
    )

    # save corr_scores for each surrogate map
    if is_surrogate:
        fname = (
            f'correlationmethod-{correlation_method}_alphalevel-{alpha}'
            f'pcacovariates-{n_pca_comp}.tsv'
        )
        corr_scores.to_csv(
            os.path.join(
                output_path,
                'smap_corr_scores',
                fname
            ), sep="\t"
        )

    # If PCA is enabled, save components and corr_scores for original
    # parcellation
    _, tail = os.path.split(parcellation_file)
    if not os.path.isfile(f"{output_path}/{tail}"):
        os.system(f"cp {parcellation_file} {output_path}/.")

    if perform_pca:
        pca_covariates_dir = os.path.join(
            output_path, 'pca_covariates', f'pca_{n_pca_comp}')
        if not os.path.isdir(pca_covariates_dir):
            os.makedirs(pca_covariates_dir)

        corr_scores.to_csv(
            os.path.join(pca_covariates_dir, f'corr_scores.tsv'), sep="\t"
        )
        significant_genes.to_csv(
            os.path.join(
                pca_covariates_dir,
                f'significant_genes.tsv'),
            sep="\t")
        for label, comp in pca_components.items():
            comp.to_filename(
                os.path.join(pca_covariates_dir, f'{label}.nii.gz')
            )

    return corr_scores, significant_genes, pca_components


def export_significant_genes(
    parcellation_file,
    marker_file,
    project_path, 
    correlation_method,
    alpha,
    partial_correlation,
):
    """
    Searches .csv files in output directory, concats all gene correlation
    scores of surrogate maps.
    Exports significant genes based on the all correlation scores for each
    surrogate map and reference map score.
    Applies FDR correction and returns significant genes.
    Creates genes.txt and genes.tsv files based on the significant genes.

    Parameters
    ----------
    parcellation_file : str or os.PathLike
        Nifti of atlas to use for parcellation.
    marker_file : str or os.PathLike
        Nifti of voxel based morphometry as e.g. outputted by CAT. Extracts
        (and returns) measures of region-wise
        gray matter density (GMD). Aggregation method is mean.
    project_path : str or os.PathLike
        path to project directory, in which the input/output folders should be
        available
    correlation_method : str
        'spearman' or 'pearson', kind of correlation to use for mass univariate
        analysis
    alpha : float
        alpha threshold for correlation analysis
    partial_correlation : bool, default = False
        Applies partial correlation for PCA covariates. Works with PCA.
        https://pingouin-stats.org/generated/pingouin.partial_corr.html
    perform_pca : bool, default = False,
        Applies Principal component analysis (PCA) using scikit-learn.
    n_pca_comp : int, default = 0
        Number of components to keep.
    custom_covariates_df : dict, default = None
        If PCA is disabled and partial_correlation is enabled,
        custom_covariates should be defined.
    
    Returns
    -------
    rejected : list
        FDR corrected significant gene list.
    """
    
    # find all correllation csv files
    logger.info(f"Import significant genes from surrogate maps...")
    smap_corr_files = []
    # for 

    # read, concat, delete p-val column from all csv
    smashed_data = []
    for f in smap_corr_files:
        smashed_data.append(pd.read_csv(f, sep="\t", index_col=0))

    smashed_concat = pd.concat(smashed_data, axis=1).drop(columns=["pval"])

    # get corr score for original atlas
    reference_data, _, _ = get_corr_scores(
        parcellation_file, marker_file, project_path,
        is_surrogate=False, save_scores=False
    )
    reference_data.drop(columns=["pval"], inplace=True)

    # apply empirical p-value formula
    stat = smashed_concat.T.values
    stat0 = reference_data.T.values
    empirical_pval = _empirical_pval(stat, stat0)
    df = pd.DataFrame(
        {
            "genes": reference_data.index,
            "pval": empirical_pval
        }
    ).set_index("genes")

    # FDR correction
    reject, corrected, *_ = multipletests(
        df["pval"], alpha=alpha, method="fdr_bh",
        is_sorted=False, returnsorted=False
    )
    df["fdr"] = corrected
    df.drop(columns=["pval"], inplace=True)
    rejected = df[reject]
    np.savetxt(
        os.path.join(
            output_path,
            "genes.txt"),
        rejected.index,
        fmt="%s")

    return rejected


def run_webgestalt(
        genes_file='genes.txt',
        r_path='/usr/bin/Rscript',
        r_arg='--vanilla',
        r_exec='./../r_script/webgestalt.r'):
    """
    Runs Webgestalt R package to conduct enrichment analysis.

    Parameters
    ----------
    genes_file : str or os.PathLike
        Significant gene list file. Should be txt.
    r_path : str or os.PathLike
        Installed R script path.
    r_arg : str
        Argument for running script in command line.
    r_exec : str or os.PathLike
        R script file for gene set enrichment analysis.

    Returns
    -------
    Creates gene enrichment analysis report.
    """

    logger.info(f"Gene enrichment analysis [webgestalt]..")

    # override - get R file from the library folder
    r_exec = r_script_file()

    if os.stat(genes_file).st_size == 0:
        logger.info(f"genes.txt file is empty.")
        return False

    pc1 = perf_counter()
    p = subprocess.Popen(
        [r_path, r_arg, r_exec, genes_file],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1,
        universal_newlines=True,
    )

    for line in p.stdout:
        print(line, end="")  # process line here
    p.wait()

    if p.returncode != 0:
        # raise ValueError(p.returncode, p.args)
        print(p.returncode, p.args)

    if p.returncode == 0:
        print(f"command {p.args} succeeded")
    elif p.returncode <= 125:
        print(
            f"command failed, exit-code={p.returncode} error: {str(p.stderr)}")
    elif p.returncode == 127:
        print(f"program not found  {str(p.stderr)}")
    else:
        pass

    pc2 = perf_counter()
    elapsed_time = ["webgestalt time elapsed", (pc2 - pc1) / 60]
    print(elapsed_time)


def _aggregate_marker(atlas, vbm, aggregation=None, limits=None):
    """
    Constructs a masker based on the input atlas_nifti, applies resampling of
    the atlas if necessary and applies the masker to
    the vbm_nifti to extract brain-imaging based vbm markers.
    So far the aggregation methods "winsorized mean", "mean" and
    "std" are supported.

    Parameters
    ----------
    atlas_nifti : niimg-like object
        Nifti of atlas to use for parcellation.
    vbm_nifti: niimg-like object
        Nifti of voxel based morphometry as e.g. outputted by CAT.
    aggregation: list
        List with strings of aggregation methods to apply. Defaults to
        aggregation = ['winsorized_mean', 'mean', 'std'].
    limits: array
        Array with lower and upper limit for the calculation of the winsorized
        mean. Only needed when 'winsorized_mean' was specified
        in aggregation. If wasn't specified defaults to [0.1, 0.1].

    Returns
    -------
    marker_aggregated : dict
        Dictionary with keys being each of the chosen aggregation methods
        and values the corresponding array with the calculated marker based on 
        the provided atlas. The array therefore as the shape of the chosen 
        number of ROIs (granularity).
    marker_func_params: dict
        Dictionary with parameters used for the aggregation function. Keys:
        respective aggregation function, values: dict with responding
        parameters
    """

    atlas_nifti = image.load_img(atlas)
    vbm_nifti = image.load_img(vbm)

    # defaults (validity is checked in _get_funcbyname())
    if aggregation is None:  # Don't put mutables as defaults, use None instead
        aggregation = ["winsorized_mean", "mean", "std", "median"]
    if limits is None:
        limits = [0.1, 0.1]

    # aggregation function parameters (validity is checked in
    # _get_funcbyname())
    agg_func_params = {"winsorized_mean": {"limits": limits}}

    # definitions
    # sort rois to be related to the order of i_roi (and get rid of 0 entry)
    rois = sorted(np.unique(image.get_data(atlas_nifti)))[1:]  # roi numbering
    n_rois = len(rois)  # granularity
    marker_aggregated = {
        x: np.ones(shape=(n_rois)) * np.nan for x in aggregation
    }

    # resample atlas if needed
    if not atlas_nifti.shape == vbm_nifti.shape:
        atlas_nifti = image.resample_to_img(
            atlas_nifti, vbm_nifti, interpolation="nearest"
        )
        logger.info("Atlas nifti was resampled to resolution of VBM nifti.")

    logger.info("make masker and apply")

    # make masker and apply
    for i_roi, roi in enumerate(rois):
        mask = image.math_img(f"img=={roi}", img=atlas_nifti)
        marker = masking.apply_mask(imgs=vbm_nifti, mask_img=mask)  # gmd per roi
        # logger.info(f'Mask applied for roi {roi}.')
        # aggregate (for all aggregation options in list)
        for agg_name in aggregation:
            # logger.info(f'Aggregate GMD in roi {roi} using {agg_name}.')
            agg_func = _get_funcbyname(
                agg_name, agg_func_params.get(
                    agg_name, None))
            marker_aggregated[agg_name][i_roi] = agg_func(gmd)
    logger.info(f"{aggregation} was computed for all {n_rois} ROIs.\n")

    return marker_aggregated, agg_func_params


def _get_funcbyname(name, func_params):
    """
    Helper function to generically apply any function. Here used to apply
    different aggregation functions for extraction of gray matter density
    (GMD).

    Parameters
    ----------
    name : str
        Name to identify the function. Currently supported names and
        corresponding functions are:
        'winsorized_mean' -> scipy.stats.mstats.winsorize
        'mean' -> np.mean
        'std' -> np.std

    func_params : dict
        Dictionary containing functions that need further parameter
        specifications. Keys are the function and values are dictionaries
        with the parameter specifications.
        E.g. 'winsorized_mean': func_params = {'limits': [0.1, 0.1]}

    Returns
    -------
    respective function with inputted (partial) parameters.
    """

    # check validity of names
    _valid_func_names = {"winsorized_mean", "mean", "std"}

    # apply functions
    if name == "winsorized_mean":
        # check validity of func_params
        limits = func_params.get("limits")
        if all((lim >= 0.0 and lim <= 1) for lim in limits):
            logger.info(f"Limits for winsorized mean are set to {limits}.")
        else:
            raise_error(
                "Limits for the winsorized mean must be between 0 and 1.")
        # partially interpret func_params
        return partial(winsorized_mean, **func_params)
    if name == "mean":
        return np.mean  # No func_params
    if name == "std":
        return np.std
    if name == "median":
        return np.median

    else:
        raise_error(
            f"Function {name} unknown. Please provide any of "
            f"{_valid_func_names}")


def winsorized_mean(data, axis=None, **win_params):
    """
    Helper function to chain winsorization and mean to compute winsorized
    mean.

    Parameters
    ----------
    data : array
        Data to calculate winsorized mean on.
    win_params : dict
        Dictionary containing the keyword arguments for the winsorize function.
        E.g. {'limits': [0.1, 0.1]}

    Returns
    -------
    Winsorized mean of the inputted data with the winsorize settings applied
    as specified in win_params.
    """

    win_dat = winsorize(data, axis=axis, **win_params)
    win_mean = win_dat.mean(axis=axis)

    return win_mean


def raise_error(error):
    logger.error(error)
    raise
