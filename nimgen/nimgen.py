"""Console script for CLI."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import argparse
import shutil

import yaml

from nimgen.logutils import configure_logging
from nimgen.pipelines.base_steps import _step_1, _step_2, _step_3, _step_4
from nimgen.pipelines.htcondor import HTCondor
from nimgen.utils import _cols_to_nifti


def parse_args():
    """Initialise the CLI script by parsing arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "nimgen CLI to run pipelines for mass univariate analysis of "
            "brain imaging data and gene expression data obtained in the "
            "Allen Human Brain Atlas"
        )
    )

    parser.add_argument(
        "-v",
        "--verbosity",
        help="Level of logging verbosity.",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )

    subparsers = parser.add_subparsers(help="sub-command help")

    # create "create" command
    create_parser = subparsers.add_parser(
        "create",
        help=(
            "Create a pipeline using a yaml configuration file."
            " Input should be the path to a valid yaml file specifying "
            "pipeline configuration."
        ),
    )
    create_parser.add_argument(
        "config_yaml", help="Path to yaml file with pipeline configuration."
    )

    create_parser.add_argument(
        "-f",
        "--force_overwrite",
        action="store_true",
        help="If activated, previous jobs directories will be overwritten.",
    )
    create_parser.set_defaults(func=create)

    # create step 1 command
    step1_parser = subparsers.add_parser(
        "step_1",
        help="Run step 1 of the nimgen pipeline to create distance matrices.",
    )
    step1_parser.add_argument(
        "parcellation_file",
        help="Parcellation file for which to create the distance matrix.",
    )
    step1_parser.set_defaults(func=step_1)

    # create step 2 command
    step2_parser = subparsers.add_parser(
        "step_2",
        help=(
            "Run step 2 of the nimgen pipeline to create null maps"
            "for a given marker/parcellation combination."
        ),
    )
    step2_parser.add_argument(
        "parcellation_file",
        help="Parcellation file for which to create the null maps.",
    )
    step2_parser.add_argument(
        "marker_file", help="Marker file for which to create the null maps."
    )
    step2_parser.add_argument(
        "n_perm", type=int, help="Number of null maps to create."
    )
    step2_parser.add_argument(
        "seed",
        type=int,
        help="Seed to use for operations involving randomness.",
    )
    step2_parser.set_defaults(func=step_2)

    # create step 3 command
    step3_parser = subparsers.add_parser(
        "step_3",
        help=(
            "Run step 3 of the nimgen pipeline to perform"
            "correlation analysis between a given marker "
            "null map and gene expression profiles."
        ),
    )
    step3_parser.add_argument(
        "parcellation_file",
        help="Parcellation file for which to run the correlation analysis.",
    )
    step3_parser.add_argument(
        "marker_file",
        help=(
            "Marker file for which to choose a given"
            "null map and perform the correlation analysis."
        ),
    )
    step3_parser.add_argument(
        "n_perm",
        type=int,
        help=("Number of null maps that were created in the second step."),
    )
    step3_parser.add_argument(
        "smap_id",
        type=int,
        help=(
            "Unique index of the null map to use"
            "for this marker/parcellation combination."
        ),
    )
    step3_parser.add_argument(
        "seed",
        type=int,
        help="Seed to use for operations involving randomness.",
    )
    step3_parser.add_argument(
        "allen_data_dir", help="Directory in which AHBA data is cached."
    )
    step3_parser.add_argument(
        "correlation_method",
        type=str,
        help="{'spearman', 'pearson' or 'dcorr'}",
    )
    step3_parser.add_argument(
        "--n_pca_covariates",
        dest="n_pca_covariates",
        type=int,
        help="Number of gene expression components to remove as covariates",
    )

    step3_parser.set_defaults(func=step_3)

    # create step 4 command
    step4_parser = subparsers.add_parser(
        "step_4",
        help=(
            "Run step 4 of the nimgen pipeline to perform"
            "correlation analysis between a given marker "
            "and gene wprofiles."
        ),
    )
    step4_parser.add_argument(
        "parcellation_file",
        help="Parcellation file for which to run the correlation analysis.",
    )
    step4_parser.add_argument(
        "marker_file",
        help=("Marker file for which to" " perform the correlation analysis."),
    )
    step4_parser.add_argument(
        "n_perm",
        type=int,
        help=("Number of null maps that were created in the second step."),
    )
    step4_parser.add_argument(
        "seed",
        type=int,
        help="Seed to use for operations involving randomness.",
    )
    step4_parser.add_argument(
        "allen_data_dir", help="Directory in which AHBA data is cached."
    )
    step4_parser.add_argument(
        "correlation_method",
        type=str,
        help="{'spearman', 'pearson' or 'dcorr'}",
    )
    step4_parser.add_argument(
        "alpha", type=float, help="Alpha level of significance."
    )
    step4_parser.add_argument(
        "--n_pca_covariates",
        dest="n_pca_covariates",
        type=int,
        help="Number of gene expression components to remove as covariates",
    )
    step4_parser.add_argument(
        "--r_path",
        type=str,
        help="Path to Rscript executable.",
    )

    step4_parser.set_defaults(func=step_4)

    # create cols_to_nifti command
    col_nii_parser = subparsers.add_parser(
        "cols_to_nifti",
        help="Convert columns of (numpy) array to nifti given a parcellation.",
    )
    col_nii_parser.add_argument(
        "array_file",
        help="Path to file with (numpy) array.",
    )
    col_nii_parser.add_argument(
        "parcellation_file",
        help="Path to parcellation nifti file",
    )
    col_nii_parser.add_argument(
        "outfolder",
        help="Folder in which to store nifti outputs.",
    )
    col_nii_parser.add_argument(
        "--n_cols",
        "-n",
        dest="n_cols",
        type=int,
        help="Number of columns to convert to nifti.",
    )
    col_nii_parser.set_defaults(func=cols_to_nifti)
    return parser.parse_args()


def yaml_to_dict(path_to_file):
    """Read yaml file with pipeline specifications.

    Parameters
    ----------
    path_to_file : str, os.PathLike
        path to yaml file

    Returns
    --------
    config_dict : dict
        dictionary with pipeline configurations

    """
    with open(path_to_file, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def create(args):
    """Take configuration dict and construct and return pipeline object.

    Parameters
    ----------
    config_dict : dict
        valid pipeline specification dictionary

    Returns
    --------
    pipeline object

    """

    config_dict = yaml_to_dict(args.config_yaml)
    valid_pipelines = {"HTCondor": HTCondor}
    assert (
        config_dict["pipeline"]["type"] in valid_pipelines.keys()
    ), f"Only pipelines implemented are {valid_pipelines.keys()}!"
    pipeline = valid_pipelines[config_dict["pipeline"]["type"]](**config_dict)
    print("Creating nimgen pipeline directory...")
    pipeline.create(args.force_overwrite)
    shutil.copy(args.config_yaml, pipeline.jobs_dir)

    return pipeline


def step_1(args):
    """Wrap first step for command line."""
    _step_1(args.parcellation_file)


def step_2(args):
    """Wrap second step for command line."""
    _step_2(
        parcellation_file=args.parcellation_file,
        marker_file=args.marker_file,
        n_perm=args.n_perm,
        seed=args.seed,
    )


def step_3(args):
    """Wrap third step for command line."""
    _step_3(
        parcellation_file=args.parcellation_file,
        marker_file=args.marker_file,
        n_perm=args.n_perm,
        smap_id=args.smap_id,
        seed=args.seed,
        allen_data_dir=args.allen_data_dir,
        correlation_method=args.correlation_method,
        n_pca_covariates=args.n_pca_covariates,
    )


def step_4(args):
    """Wrap fourth step for command line."""
    if args.r_path is not None:
        _step_4(
            parcellation_file=args.parcellation_file,
            marker_file=args.marker_file,
            n_perm=args.n_perm,
            seed=args.seed,
            allen_data_dir=args.allen_data_dir,
            correlation_method=args.correlation_method,
            alpha=args.alpha,
            r_path=args.r_path,
        )
    else:
        _step_4(
            parcellation_file=args.parcellation_file,
            marker_file=args.marker_file,
            n_perm=args.n_perm,
            seed=args.seed,
            allen_data_dir=args.allen_data_dir,
            correlation_method=args.correlation_method,
            alpha=args.alpha,
        )


def cols_to_nifti(args):
    """Wrap _cols_to_nii function for command line."""
    _cols_to_nifti(
        filename=args.array_file,
        parcellation_file=args.parcellation_file,
        outfolder=args.outfolder,
        n_cols=args.n_cols,
    )


def main():
    """Run nimgen CLI."""
    args = parse_args()
    configure_logging(args.verbosity)
    args.func(args)
