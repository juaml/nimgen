#!/usr/bin/bash
source /data/project/nimgen/nimgen_validation_study/venv/nimgen_validation/bin/activate
python3 $@
deactivate
