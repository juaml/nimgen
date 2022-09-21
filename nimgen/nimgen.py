"""Console script for CLI."""

# Authors: Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import argparse
import os

import yaml

from nimgen.pipelines.htcondor import HTCondor


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
        "--create",
        "-c",
        dest="create",
        help=(
            "create a pipeline using a yaml configuration file."
            "Input should be the path to a valid yaml file specifying "
            "pipeline configuration."
        ),
    )
    parser.add_argument(
        "--run",
        "-r",
        dest="run",
        help=(
            "Create (if it has not been created yet) and run a pipeline"
            " using a yaml configuration file."
            "Input should be the path to a valid yaml file specifying "
            "pipeline configuration."
        ),
    )

    return parser.parse_args()


def validate_args(args):
    """Check that values for keyword arguments are valid.

    Parameters
    ----------
    args : args
        arguments parsed by argparse.ArgumentParser

    Returns
    --------
    path_to_yaml : str
        path to a yaml file determining pipeline configuration

    """
    if args.create is not None:
        if not os.path.isfile(args.create):
            raise FileNotFoundError(f"{args.create} not found!")
        return args.create
    elif args.run is not None:
        if not os.path.isfile(args.create):
            raise FileNotFoundError(f"{args.run} not found!")
        return args.run


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


def create_pipeline(**config_dict):
    """Take configuration dict and construct and return pipeline object.

    Parameters
    ----------
    config_dict : dict
        valid pipeline specification dictionary

    Returns
    --------
    pipeline object

    """
    valid_pipelines = {"HTCondor": HTCondor}
    assert (
        config_dict["pipeline"]["type"] in valid_pipelines.keys()
    ), f"Only pipelines implemented are {valid_pipelines.keys()}!"
    pipeline = valid_pipelines[config_dict["pipeline"]["type"]](**config_dict)
    pipeline.create()
    return pipeline


def main():
    """Run nimgen CLI."""
    args = parse_args()
    print("You are running the nimgen CLI!")
    yaml_file = validate_args(args)
    config_dict = yaml_to_dict(yaml_file)
    create_pipeline(**config_dict)
