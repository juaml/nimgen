"""Utility functions commonly used in a variety of modules."""

# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Leonard Sasse <l.sasse@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
# License: AGPL

import logging
import os
import subprocess
import sys
import warnings
from distutils.version import LooseVersion
from pathlib import Path


from neuromaps.parcellate import Parcellater
import numpy as np
import pandas as pd
from nilearn import image


logger = logging.getLogger("nimgen")

# logging.basicConfig(format=format)
stream_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
stream_handler.setFormatter(formatter)

logger.propagate = False
logger.setLevel(logging.INFO)
logger.addHandler(stream_handler)


def _get_git_head(path):
    """Aux function to read HEAD from git."""
    if not path.exists():
        raise_error("This path does not exist: {}".format(path))
    command = ("cd {gitpath}; " "git rev-parse --verify HEAD").format(
        gitpath=path
    )
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    del process
    return proc_stdout


def get_versions(sys):
    """Get versions for each module.

    If it's a git-installed package, get the git hash too.

    Parameters
    ----------
    sys : module
        The sys module object.
    Returns
    -------
    module_versions : dict
        The module names and corresponding versions.
    """
    module_versions = {}
    for name, module in sys.modules.items():
        if "." in name:
            continue
        if name in ["_curses", "_glmnet"]:
            continue
        module_version = LooseVersion(getattr(module, "__version__", None))
        module_version = getattr(module_version, "vstring", None)
        if module_version is None:
            module_version = None
        elif "git" in module_version:
            git_path = Path(module.__file__).resolve().parent
            head = _get_git_head(git_path)
            module_version += "-HEAD:{}".format(head)

        module_versions[name] = module_version
    return module_versions


def _safe_log(versions, name):
    if name in versions:
        logger.info(f"{name}: {versions[name]}")


def log_versions():
    """Log versions of the core libraries, for reproducibility purposes."""
    versions = get_versions(sys)
    logger.info("===== Lib Versions =====")
    _safe_log(versions, "numpy")
    _safe_log(versions, "scipy")
    _safe_log(versions, "sklearn")
    _safe_log(versions, "pandas")
    _safe_log(versions, "abagen")
    _safe_log(versions, "neuromaps")
    _safe_log(versions, "brainsmash")
    logger.info("========================")


def raise_error(msg, klass=ValueError):
    """Raise an error, but first log it.

    Parameters
    ----------
    msg : str
        Error message
    klass : class of the error to raise. Defaults to ValueError
    """
    logger.error(msg)
    raise klass(msg)


def warn(msg, category=RuntimeWarning):
    """Warn, but first log it.

    Parameters
    ----------
    msg : str
        Warning message
    category : instance of Warning
        The warning class. Defaults to ``RuntimeWarning``.
    """
    logger.warning(msg)
    warnings.warn(msg, category=category)


def remove_nii_extensions(nii_file):
    """Remove file extension from .nii or nii.gz path."""
    nii_name = Path(nii_file)

    while nii_name.suffix in {".nii", ".gz"}:
        nii_name = nii_name.with_suffix("")

    return str(nii_name)


def read_csv_tsv(path):
    """Read both csv and tsv, file type inferred by extension."""
    _, ext = os.path.splitext(path)
    extensions = {".csv": ",", ".tsv": "\t"}

    return pd.read_csv(path, index_col=0, sep=extensions[ext])


def covariates_to_nifti(parcellation, covariates_df):
    """Create nifti images for given PCA covariates.

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
    if os.path.isfile(parcellation):
        parcellation = image.load_img(parcellation)

    parcellation_array = np.array(parcellation.dataobj)
    covariate_niftis = {}
    for covariate_label, covariate in covariates_df.items():
        null_mat = np.zeros(parcellation_array.shape)
        null_mat[parcellation_array == 0] = -20000
        marker_min = np.min(covariate)
        marker_max = np.max(covariate)
        for label, value in covariate.iteritems():
            null_mat[parcellation_array == label + 1] = value

        pc_nii = image.new_img_like(parcellation, null_mat)
        pc_nii.header["cal_min"] = marker_min
        pc_nii.header["cal_max"] = marker_max
        covariate_niftis[covariate_label] = pc_nii

    return covariate_niftis


def _cols_to_nifti(
    filename,
    parcellation_file,
    outfolder,
    n_cols=None
):
    """Convert columns of an array into nifti files.
    
    Parameters
    ----------
    filename : str or os.PathLike
        File name of numpy (.npy) file with array.
    parcellation_file : str or os.PathLike
        File name of parcellation to use.
    outfolder : str or os.PathLike
        Directory in which to save output.
    n_cols : int or None
        How many columns to convert.
        If None, it will do all.
    
    """
    # Create array specific outpath.
    filename = Path(filename)
    array_name = Path(filename.name).stem
    outfolder = Path(outfolder)
    assert outfolder.exists(), f"{outfolder} does not exist."
    outfolder = outfolder / array_name
    outfolder.mkdir(exist_ok=True)
    
    # load array and convert
    array = np.load(filename)
    _, cols = array.shape
    if n_cols is None:
        n_cols = cols

    assert isinstance(n_cols, int), "n_cols should be an integer."
    for idx in range(n_cols):
        column = array[:, idx]
        nm_parcellater = Parcellater(parcellation_file, space="MNI152")
        nifti = nm_parcellater.inverse_transform(column)
        outfile = outfolder / f"{array_name}_column-{idx}.nii.gz"
        nifti.to_filename(outfile)
    

def _read_sign_genes(sign_genes):
    if isinstance(sign_genes, pd.DataFrame):
        sign_genes = sign_genes.index
    elif isinstance(sign_genes, list):
        return sign_genes
    elif os.path.isfile(sign_genes):
        _, ext = os.path.splitext(sign_genes)
        if ext in [".csv", ".tsv"]:
            extensions = {".csv": ",", ".tsv": "\t"}
            sign_genes = pd.read_csv(
                sign_genes,
                header=None,
                index_col=0,
                dtype=str,
                sep=extensions[ext],
            ).index
        elif ext in [".txt"]:
            sign_genes = list(np.loadtxt(sign_genes, dtype=str))
    else:
        raise ValueError(
            "'sign_genes' should be a pd.DataFrame," " .csv/.tsv or .txt file!"
        )

    return list(sign_genes)
