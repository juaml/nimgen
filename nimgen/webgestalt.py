# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
# License: AGPL
import tempfile
from pathlib import Path

_URL = "http://www.webgestalt.org/option.php"


def _post_selenium(driver, url, data):
    input_template = ('{k} <input type="text" name="{k}"'
                      ' id="{k}" value="{v}"><BR />\n')
    multi_template = ('{k} <textarea name="{k}"'
                      'id="{k}">{v}</textarea><BR />\n')
    inputs = ""
    if data:
        for k, v in data.items():
            if isinstance(v, list):
                inputs += multi_template.format(k=k, v='\n'.join(v))
            else:
                inputs += input_template.format(k=k, v=v)
    html = (f'<html><body>\n<form action="{url}" method="post" id="formid">\n'
            f'{inputs}<input type="submit" id="inputbox">\n'
            f'</form></body></html>')

    with tempfile.TemporaryDirectory() as tmpdir:
        html_file = Path(tmpdir) / 'temp.html'
        with open(html_file, "w") as text_file:
            text_file.write(html)

        driver.get(f"file://{html_file}")
        driver.find_element_by_id('inputbox').click()


def webgestalt(genes, enrich_db_category, enrich_db_name, enrich_method='ORA',
               id_type='genesymbol', refset='agilent_sureprint_g3_ge_8x60k',
               organism='hsapiens', sig_method='fdr', sig_value=0.05,
               browser='chrome'):
    """Set-up an analysis to be performed in webgestalt.
    See http://www.webgestalt.org

    Requires selenium and the browser driver installed. 
    See https://juaml.github.io/nimgen/installation.html

    Parameters
    ----------
    genes : list(str)
        A list will all the genes that will take part in the webgestalt
        analysis.
    enrich_db_category : str
        Functional database category. Options are geneontology, pathway,
        network, disease, drug, phenotype, chromosomalLocation,
        community-contributed, others.
    enrich_db_name : str
        Depends on enrich_db_category. Check http://www.webgestalt.org for
        options.
    enrich_method : str
        Options are: ORA, GSEA, NTA. Defaults to 'ORA'.
    id_type: str
        Gene ID Type of the gene list. Check http://www.webgestalt.org for
        options. Defaulkts to 'genesymbol'.
    refset : str
        Built-in reference set or microarray platform (only when
        `enrich_method='ORA'`). Defaults to 'agilent_sureprint_g3_ge_8x60k'.
    organism : str
        Organism of interest. Options are athaliana, btaurus, celegans,
        cfamiliaris, dmelanogaster, drerio, ggallus, hsapiens, mmusculus,
        rnorvegicus, scerevisiae, sscrofa, others. Defaults to 'hsapiens'.
    sig_method : str
        Methods to use for identifying the enriched categories. 'fdr' means the
        enriched categories will be identified based on FDR threshold and
        'top' means the categories will be first ranked based on the FDR and
        then the top most significant categories will be selected. Defaults to
        'fdr'.
    sig_value : numerical
        Threshold for the significance method. Defaults to 0.05.
    browser : str
        Browser to use. Supported options are 'chrome' (default) and 'firefox'.

    Raises
    ------
    ValueError
        If there is a problem with the browser and/or selenium

    """
    from selenium import webdriver
    driver = None
    if browser == 'chrome':
        driver = webdriver.Chrome()
    elif browser == 'firefox':
        driver = webdriver.Firefox()
    else:
        raise ValueError('Supported browsers are Chrome (browser=chrome) '
                         'and Firefox (browser=firefox)')
    params = {
        'organism': organism,
        'enrich_method': enrich_method,
        'enriched_database_category': enrich_db_category,
        'enriched_database_name': enrich_db_name,
        'id_type': id_type,
        'ref_set': refset,
        'gene_list': genes,
        'sig_method': sig_method,
        'sig_value': sig_value
    }
    _post_selenium(driver, _URL, params)
