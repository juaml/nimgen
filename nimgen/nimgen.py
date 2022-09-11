import argparse
import os
from ptpython.ipython import embed


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "nimgen CLI to run pipelines for mass univariate analysis of "
            "brain imaging data and gene expression data obtained in the "
            "Allen Human Brain Atlas"
        )
    )
    parser.add_argument(
        "--create", "-c",
        dest="create",
        help=(
            "create a pipeline using a yaml configuration file."
            "Input should be the path to a valid yaml file specifying "
            "pipeline configuration."
        )
    )
    parser.add_argument(
        "--run", "-r",
        dest="run",
        help=(
            "Create (if it has not been created yet) and run a pipeline"
            " using a yaml configuration file."
            "Input should be the path to a valid yaml file specifying "
            "pipeline configuration."
        )
    )

    return parser.parse_args()


def validate_args(args):

    if (args.create is not None) and (args.run is not None):
        assert args.create == args.run, (
            "It is recommended you use either --create or --run, not both."
            " By default run will also create a pipeline directory if it "
            "doesn't exist yet!"
        )
    elif args.create is not None:
        if not os.path.isfile(args.create):
            raise FileNotFoundError(f"{args.create} not found!")
    elif args.run is not None:
        if not os.path.isfile(args.create):
            raise FileNotFoundError(f"{args.run} not found!")


def yaml_to_dict(path_to_file):
    pass


def main():

    args = parse_args()
    print("You are running the nimgen CLI!")
    # validate_args(args)
    embed()
