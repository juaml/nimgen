
STEP_ONE_FSTRING = """
#!/usr/bin/env python3

import sys
from nimgen.pipelines import base_steps

# {}
if __name__ == "__main__":
    parcellation_file = sys.argv[1]
    base_steps.step_1(parcellation_file)
"""

STEP_TWO_FSTRING = """
#!/usr/bin/env python3

import sys
from nimgen.pipelines import base_steps


if __name__ == "__main__":
    name_markers_dir = "{}"
    name_output_dir = "{}"
    allen_data_dir = "{}"

    parcellation_file = sys.argv[1]
    marker_file = sys.argv[2]
    smap_id = sys.argv[3]
    correlation_method = sys.argv[4]
    alpha = float(sys.argv[5])
    n_pca_covariates = literal_eval(sys.argv[6])
    partial_correlation = literal_eval(sys.argv[7])

    base_steps.step_2(
        parcellation_file=parcellation_file,
        marker_file=marker_file,
        name_markers_dir=name_markers_dir,
        name_output_dir=name_output_dir,
        smap_id=smap_id,
        allen_data_dir=allen_data_dir,
        correlation_method=correlation_method,
        alpha=alpha,
        n_pca_covariates=n_pca_covariates,
        partial_correlation=partial_correlation,
    )
"""

STEP_THREE_FSTRING = """
#!/usr/bin/env python3

import sys
from nimgen.pipelines import base_steps


if __name__ == "__main__":
    base_steps.step_3(
        parcellation_file=sys.argv[1],
        marker_file=sys.argv[2],
        name_markers_dir="{}",
        name_output_dir="{}",
        allen_data_dir="{}",
        r_path="{}",
        correlation_method=sys.argv[3],
        alpha=float(sys.argv[4]),
        n_pca_covariates=literal_eval(sys.argv[5]),
        partial_correlation=literal_eval(sys.argv[6]),
    )

"""
