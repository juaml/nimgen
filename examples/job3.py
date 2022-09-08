#!/usr/bin/env python3
from nimgen import export_significant_genes, run_webgestalt, create_sample_path, correlated_gene_expression, gene_coexpression
import os
import sys

PARCEL_FILE = sys.argv[1]
MARKER_FILE = sys.argv[2]

_, output_path, _ = create_sample_path(PARCEL_FILE, MARKER_FILE)
genes_file = os.path.join(output_path, 'genes.txt')

# Searches .csv files in output directory, concats all gene correlation scores of surrogate maps.
# Exports significant genes based on the all correlation scores for each surrogate map.
# Apply FDR and create genes.txt
export_significant_genes = export_significant_genes(PARCEL_FILE, MARKER_FILE)

# Run gene set enrichment analysis using genes.txt via webgestalt R package
run_webgestalt(r_path='/home/ydemirtas/miniconda3/envs/r_env/bin/Rscript', genes_file=genes_file)

# calculate gene co-expression based on significant genes
corr_gene_exp_matrix = correlated_gene_expression( PARCEL_FILE, genes_file, metric='spearman')
coexp_matrix = gene_coexpression(PARCEL_FILE, genes_file, metric='spearman')


# save as csv in the output path
coexp_matrix.to_csv(os.path.join(output_path, f'gene_by_gene_matrix.csv'), sep="\t")   
corr_gene_exp_matrix.to_csv(os.path.join(output_path, f'region_by_region_matrix.csv'), sep="\t")   


