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


def webgestalt(genes, enrich_db_category, enrich_db_name,
               id_type='genesymbol', refset='agilent_sureprint_g3_ge_8x60k',
               organism='hsapiens', enrich_method='ORA', browser='chrome'):
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
        'gene_list': genes
    }
    _post_selenium(driver, _URL, params)
