
STEP_ONE_FSTRING = """
#!/usr/bin/env python3
from nimgen import (
    export_voxel_coordinates, generate_distance_matrices, create_sample_path
)
from nimgen.pipelines.htcondor import remove_nii_extensions
import sys
import os
import glob


PARCELLATION_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]

N_SURROGATE_MAPS = {}
project_path="{}"

if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
    raise ValueError('Input file not found.')

# create sample path specific for marker and atlas
_, parc_name = os.path.split(remove_nii_extensions(PARCELLATION_FILE))

parcellation_path = os.path.join({}, parc_name)

# check surrogate maps for given atlas
smaps_dir = os.path.join(parcellation_path, 'smaps', '*.nii')
smaps = glob.glob(smaps_dir)

# if number of surrogate maps are enough, do not create distance matrix
if len(smaps) >= int(N_SURROGATE_MAPS):
    print(
        f'matrix file creation skipped.'
        f' In the surrogate map folder there are already enough'
        f' surrogate maps.'
    )
    pass
else:
    coord_file, parcel_file = export_voxel_coordinates(
        PARCEL_FILE, MARKER_FILE)
    matrix_files = generate_distance_matrices(PARCEL_FILE, MARKER_FILE)

"""

STEP_TWO_FSTRING = """
#!/usr/bin/env python3

from nimgen import generate_surrogate_map, get_corr_scores
import sys
import os


parcel_file = sys.argv[1]
marker_file = sys.argv[2]
smap_id = sys.argv[3]
correlation_method = sys.argv[4]
alpha = float(sys.argv[5])
n_pca_covariates = int(sys.argv[6])

project_path="{}"

if not os.path.isfile(parcel_file) or not os.path.isfile(parcel_file):
    raise ValueError('Input file not found.')

# generate surrogate map for given atlas
surrogate_map = generate_surrogate_map(
    parcel_file, marker_file, smap_id, project_path=project_path
)

# perform correlation analysis for given surrogate map and marker
corr_scores, sign_genes, pca_components = get_corr_scores(
    parcellation_file=surrogate_map,
    marker_file=marker_file,
    project_path=project_path,
    correlation_method=correlation_method,
    alpha=alpha,
)

# save corr_scores for each surrogate map
if is_surrogate:
    fname = (
        f'correlationmethod-{{correlation_method}}_alphalevel-{{alpha}}'
        f'pcacovariates-{{n_pca_comp}}-smapid-{{smap-id}}.tsv'
    )
    corr_scores.to_csv(
        os.path.join(
            output_path,
            'smap_corr_scores',
            fname
        ), sep="\t"
    )


"""

STEP_THREE_FSTRING = """
#!/usr/bin/env python3
from nimgen import (
    export_significant_genes, run_webgestalt, create_sample_path,
    correlated_gene_expression, gene_coexpression
)
import os
import sys

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]
project_path = "{}"
markers_path = "{}"

marker_head, marker_tail = os.path.split(MARKER_FILE)
marker_name, _ = os.path.splitext(marker_tail)

parcel_head, parcel_tail = os.path.split(PARCEL_FILE)
parcel_name, _ = os.path.split(parcel_tail)
output_path = os.path.join(
    "{}", parcel_name, marker_head, marker_name,
)
MARKER_FILE_ABS = os.path.join(markers_path, MARKER_FILE)

genes_file = os.path.join(output_path, 'genes.txt')

# Searches .csv files in output directory
# concats all gene correlation scores of surrogate maps.
# Exports significant genes based on the all
# correlation scores for each surrogate map.
# Apply FDR and create genes.txt

export_significant_genes = export_significant_genes(
    PARCEL_FILE, MARKER_FILE, project_path
)

# Run gene set enrichment analysis using genes.txt via webgestalt R package
run_webgestalt(
    r_path="{}",
    genes_file=genes_file
)

# calculate gene co-expression based on significant genes
corr_gene_exp_matrix = correlated_gene_expression(
    PARCEL_FILE, "all", metric='spearman'
)
coexp_matrix = gene_coexpression(PARCEL_FILE, genes_file, metric='spearman')

# calculate gene co-expression based on all genes
corr_all_gene_exp_matrix = correlated_gene_expression(
    PARCEL_FILE, "all", metric='spearman'
)
coexp_all_matrix = gene_coexpression(
    PARCEL_FILE, genes_file, metric='spearman'
)

# save as csv in the output path
coexp_matrix.to_csv(
    os.path.join(
        output_path,
        f'significant-genes_gene_by_gene_matrix.tsv'
    ), sep="\t"
)
corr_gene_exp_matrix.to_csv(
    os.path.join(
        output_path,
        f'significant-genes_region_by_region_matrix.tsv'
    ), sep="\t"
)

# save matrices based on all genes
coexp_all_matrix.to_csv(
    os.path.join(
        output_path,
        f'all-genes_gene_by_gene_matrix.tsv'
    ), sep="\t"
)
corr_all_gene_exp_matrix.to_csv(
    os.path.join(
        output_path,
        f'all-genes_region_by_region_matrix.tsv'
    ), sep="\t"
)
"""

STEP_FOUR_FSTRING = """
#!/usr/bin/env python3
from nimgen import create_sample_path, get_corr_scores
from nimgen.expressions import correlated_gene_expression, gene_coexpression
import os
import sys

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]

PARTIAL_CORR = bool(sys.argv[3].split('=')[1])
PERFORM_PCA = bool(sys.argv[4].split('=')[1])
N_PCA_COMP = int(sys.argv[5].split('=')[1])

project_path="{}"

_, output_path, _ = create_sample_path(PARCEL_FILE, MARKER_FILE, project_path)

corr_scores, significant_genes, pca_components = get_corr_scores(
    PARCEL_FILE,
    MARKER_FILE,
    project_path,
    partial_correlation=PARTIAL_CORR,
    perform_pca=PERFORM_PCA,
    n_pca_comp=N_PCA_COMP,
    custom_covariates_df=None,
    is_surrogate=False,
)

corr_gene_exp_matrix = correlated_gene_expression(
    PARCEL_FILE, significant_genes, metric='spearman'
)
coexp_matrix = gene_coexpression(
    PARCEL_FILE, significant_genes, metric='spearman'
)

# save as csv in the output path
coexp_matrix.to_csv(
    os.path.join(
        output_path, "pca_covariates", f"pca_{{N_PCA_COMP}}",
        f'significant-genes_pca_{{N_PCA_COMP}}_gene_by_gene_matrix.tsv'
    ), sep="\t"
)
corr_gene_exp_matrix.to_csv(
    os.path.join(
        output_path, "pca_covariates", f"pca_{{N_PCA_COMP}}",
        f'significant-genes_pca_{{N_PCA_COMP}}_region_by_region_matrix.tsv'
    ), sep="\t"
)

"""

