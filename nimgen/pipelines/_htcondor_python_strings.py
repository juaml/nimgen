RUN_IN_VENV = """#!/usr/bin/bash
source {venv}/bin/activate
$@
deactivate
"""

RUN_CONDA = """#!/usr/bin/bash
eval "$(conda shell.bash hook)"
conda activate {conda_env}
$@
conda deactivate
"""
""
