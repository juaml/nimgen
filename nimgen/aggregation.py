from .expressions import *
from .utils import logger
from pathlib import Path
from slugify import slugify
import subprocess
from nilearn import image, masking
from functools import partial
from scipy.stats.mstats import winsorize
from time import perf_counter
import uuid
import pathlib
from joblib import Parallel, delayed
import multiprocessing

def create_sample_dir(project_dir, parcellation_file):
    output_dir = os.path.join(project_dir, f"output")
    sample_id = uuid.uuid4().hex.upper()[0:6]
    sample_name = "_".join([pathlib.Path(parcellation_file).stem, sample_id])
    sample_path = os.path.join(output_dir, sample_name)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    if not os.path.isdir(sample_path):
        os.mkdir(sample_path)
    return sample_path



def get_pvalues(
    atlas_file,
    marker_file,
    aggregation_methods,
    output_dir,
    allen_data_dir="allen_data",
    **kwargs,
):
    # aggregation methoduna gore ROI sayisi kadar weights dizisi olusturur, onu import eder marker dosyasi gibi
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    atlas = image.load_img(atlas_file)
    marker = image.load_img(marker_file)
    gmd_aggregated, agg_func_params = get_gmd(
        atlas, marker, aggregation=aggregation_methods
    )

    all_pvalues = {}
    for agg_name in gmd_aggregated:
        sample_name = f"{slugify(Path(atlas_file).stem)}_{slugify(Path(marker_file).stem)}_{slugify(agg_name)}"
        sample_dir = os.path.join(output_dir, sample_name)
        if not os.path.isdir(sample_dir):
            os.mkdir(sample_dir)
        logger.info(f"save files to: {sample_name}")
        pvalues, sign = get_gene_expression(
            gmd_aggregated[agg_name],
            atlas_file,
            allen_data_dir=allen_data_dir,
            save_expressions=True,
            multiple_correction=None,
        )
        pvalues.to_csv(os.path.join(sample_dir, "genes.csv"), sep="\t")
        all_pvalues[agg_name] = pvalues

    return all_pvalues, sign



def get_pvalues_parallel(atlases, num_cores, **kwargs):
    #kwargs['output_dir'] = os.path.join(kwargs['output_dir'], 'smashed_expressions')    
    output_dir = kwargs.get('output_dir', 1)
    pc1 = perf_counter()
    smashed_data = Parallel(n_jobs=num_cores)(
        delayed(get_pvalues)(atlas, **kwargs) for atlas in atlases
    )
    # joblib returns in dict
    smashed_data = [i[0] for i in smashed_data]
    smashed_data = [s['mean'] for s in smashed_data]
    pc2 = perf_counter()
    elapsed_time= ["parallel.get_pvalues", (pc2 - pc1) / 60]
    np.savetxt(os.path.join(output_dir,'elapsed_time2.txt'), elapsed_time, fmt='%s')
    print(elapsed_time)
    np.save(os.path.join(kwargs['output_dir'], "smashed_data.npy"), smashed_data)
    return smashed_data



def save_significant_genes(sign_genes, output_dir):
    for posneg in ["positive", "negative"]:
        if posneg == "negative":
            posneg_genes = sign_genes[sign_genes["r_score"] < 0]
        else:
            posneg_genes = sign_genes[sign_genes["r_score"] > 0]
        posneg_genes.to_csv(
            os.path.join(output_dir, f"{posneg}_sign_genes.tsv"), sep="\t"
        )

    sign_genes.to_csv(os.path.join(output_dir, "sign_genes.tsv"), sep="\t")


def run_webgestalt(
    genes_file,
    executable_r_file="r/webg.r",
    r_path="/usr/local/bin/Rscript",
    r_arg="--vanilla",
):
    
    pc1 = perf_counter()
    p = subprocess.Popen(
        [r_path, r_arg, executable_r_file, genes_file],
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
        raise print(p.returncode, p.args)

    pc2 = perf_counter()
    elapsed_time= ["webgestalt", (pc2 - pc1) / 60]
    np.savetxt(os.path.join(os.path.dirname(genes_file),'elapsed_time3.txt'), elapsed_time, fmt='%s')
    print(elapsed_time)


def get_gmd(atlas_nifti, vbm_nifti, aggregation=None, limits=None):
    """
    Constructs a masker based on the input atlas_nifti, applies resampling of
    the atlas if necessary and applies the masker to
    the vbm_nifti to extract (and return) measures of region-wise gray matter
    density (GMD). So far the aggregtaion methods "winsorized mean", "mean" and
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
    gmd_aggregated : dict
        Dictionary with keys being each of the chosen aggregation methods
        and values the corresponding array with the calculated GMD based on the
        provided atlas. The array therefore as the shape of the chosen number
        of ROIs (granularity).
    agg_func_params: dict
        Dictionary with parameters used for the aggregation function. Keys:
        respective aggregation function, values: dict with responding
        parameters
    """

    # defaults (validity is checked in _get_funcbyname())
    if aggregation is None:  # Don't put mutables as defaults, use None instead
        aggregation = ["winsorized_mean", "mean", "std", "median"]
    if limits is None:
        limits = [0.1, 0.1]

    # aggregation function parameters (validity is checked in _get_funcbyname())
    agg_func_params = {"winsorized_mean": {"limits": limits}}

    # definitions
    # sort rois to be related to the order of i_roi (and get rid of 0 entry)
    rois = sorted(np.unique(image.get_data(atlas_nifti)))[1:]  # roi numbering
    n_rois = len(rois)  # granularity
    gmd_aggregated = {x: np.ones(shape=(n_rois)) * np.nan for x in aggregation}

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
        gmd = masking.apply_mask(imgs=vbm_nifti, mask_img=mask)  # gmd per roi
        # logger.info(f'Mask applied for roi {roi}.')
        # aggregate (for all aggregation options in list)
        for agg_name in aggregation:
            # logger.info(f'Aggregate GMD in roi {roi} using {agg_name}.')
            agg_func = _get_funcbyname(agg_name, agg_func_params.get(agg_name, None))
            gmd_aggregated[agg_name][i_roi] = agg_func(gmd)
    logger.info(f"{aggregation} was computed for all {n_rois} ROIs.\n")

    return gmd_aggregated, agg_func_params


def _get_funcbyname(name, func_params):
    """
    Helper function to generically apply any function. Here used to apply
    different aggregation functions for extraction of gray matter density (GMD).

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
            raise_error("Limits for the winsorized mean must be between 0 and 1.")
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
            f"Function {name} unknown. Please provide any of " f"{_valid_func_names}"
        )


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
