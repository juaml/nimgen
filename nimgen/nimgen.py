import argparse
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
    return parser.parse_args()


def yaml_to_dict(path_to_file):
    pass

def create(configs):
    pass

def main():

    args = parse_args()
    configs = yaml_to_dict(args.create)
    pipeline = create(configs)