"""Module to interface with R and Webgestalt for Gene Enrichment Analysis."""

# Authors: Yasir Demirtaş <tyasird@gmail.com>
#          Leonard Sasse <l.sasse@fz-juelich.de>
# License: AGPL

import os
import subprocess
from time import perf_counter

from .utils import logger


def r_script_file():
    """Return webgestalt.r abspath.

    It is generally inside in the r_script folder in the package.
    """
    r_file = (
        os.path.dirname(os.path.abspath(__file__))
        + "/../r_script/webgestalt.r"
    )
    return r_file


def run_webgestalt(
    genes_file="genes.txt",
    r_path="/usr/bin/Rscript",
    r_arg="--vanilla",
    r_exec="./../r_script/webgestalt.r",
):
    """Run Webgestalt R package to conduct enrichment analysis.

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

    logger.info("Gene enrichment analysis [webgestalt]..")

    # override - get R file from the library folder
    r_exec = r_script_file()

    if os.stat(genes_file).st_size == 0:
        logger.info("genes.txt file is empty.")
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
            f"command failed, exit-code={p.returncode} error: {str(p.stderr)}"
        )
    elif p.returncode == 127:
        print(f"program not found  {str(p.stderr)}")
    else:
        pass

    pc2 = perf_counter()
    elapsed_time = ["webgestalt time elapsed", (pc2 - pc1) / 60]
    print(elapsed_time)
