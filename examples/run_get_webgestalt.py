# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
# License: AGPL
import nimgen
import pandas as pd

atlas = './data/Power_5mm.nii'
weights = pd.read_csv('./data/rehoavg.txt', header=None)
allen_data_dir = './data/allen'

genes = nimgen.get_gene_expression(
    weights, atlas, allen_data_dir=allen_data_dir, save_expressions=True)

nimgen.webgestalt.webgestalt(genes, enrich_db_category='geneontology',
                             enrich_db_name='Cellular_Component')
