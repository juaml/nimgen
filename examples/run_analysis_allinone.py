#!/usr/bin/env python3
import nimgen
import sys
import os
import glob

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]
N_SURROGATE_MAPS = int(sys.argv[3].split('=')[1])
PARTIAL_CORR = bool(sys.argv[4].split('=')[1])
PERFORM_PCA = bool(sys.argv[5].split('=')[1])
N_PCA_COMP = sys.argv[6].split('=')[1]


if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
    raise ValueError('Input file not found.')

# create sample path specific for marker and atlas
_, output_path, parcellation_path = nimgen.create_sample_path(PARCEL_FILE, MARKER_FILE, True)

# check surrogate maps for given atlas
smaps_dir = os.path.join(parcellation_path, 'smaps', '*.nii')
smaps = glob.glob(smaps_dir)

# if number of surrogate maps are enough, do not create distance matrix
if len(smaps) >= int(N_SURROGATE_MAPS):
    print(f'matrix file creation skipped. In the surrogate map folder there are allready {len(smaps)} surrogate maps.')
    pass
else:
    coord_file, parcel_file = nimgen.export_voxel_coordinates(PARCEL_FILE, MARKER_FILE)
    matrix_files = nimgen.generate_distance_matrices(PARCEL_FILE, MARKER_FILE)


if not os.path.isfile(PARCEL_FILE) or not os.path.isfile(MARKER_FILE):
    raise ValueError('Input file not found.')


for SMAP_ID in range(N_SURROGATE_MAPS):

    # generate surrogate map for given atlas
    surrogate_map = nimgen.generate_surrogate_map(PARCEL_FILE, MARKER_FILE, SMAP_ID)

    # perform correlation analysis for given surrogate map and marker
    nimgen.get_corr_scores(
        parcellation_file=surrogate_map,
        marker_file=MARKER_FILE,
    )


genes_file = os.path.join(output_path, 'genes.txt')

# Searches .csv files in output directory, concats all gene correlation scores of surrogate maps.
# Exports significant genes based on the all correlation scores for each surrogate map.
# Apply FDR and create genes.txt
export_significant_genes = nimgen.export_significant_genes(PARCEL_FILE, MARKER_FILE)

# Run gene set enrichment analysis using genes.txt via webgestalt R package
nimgen.run_webgestalt(r_path='/home/ydemirtas/miniconda3/envs/r_env/bin/Rscript', genes_file=genes_file)

# calculate gene co-expression based on significant genes
corr_gene_exp_matrix = nimgen.correlated_gene_expression( PARCEL_FILE, genes_file, metric='spearman')
coexp_matrix = nimgen.gene_coexpression(PARCEL_FILE, genes_file, metric='spearman')


# save as csv in the output path
coexp_matrix.to_csv(os.path.join(output_path, f'gene_by_gene_matrix.csv'), sep="\t")   
corr_gene_exp_matrix.to_csv(os.path.join(output_path, f'region_by_region_matrix.csv'), sep="\t")   


# for only original atlas perform PCA 
corr_scores, significant_genes, pca_components = nimgen.get_corr_scores(
        PARCEL_FILE,
        MARKER_FILE,
        partial_correlation=PARTIAL_CORR,
        perform_pca=PERFORM_PCA,
        n_pca_comp=N_PCA_COMP,
        custom_covariates_df=None,
        is_surrogate=False,
)
