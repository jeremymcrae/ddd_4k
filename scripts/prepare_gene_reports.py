"""
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import getpass

import psycopg2
from tabulate import tabulate
import markdown
from weasyprint import HTML

from ddd_4k.constants import PHENOTYPES, SANGER_IDS
from ddd_4k.load_files import open_phenotypes

def get_variant_details():
    """
    """
    
    pwd = getpass.getpass("DDD database password: ")
    conn = psycopg2.connect(database="ddd_prod", user="ddd_login_ro", \
        host="ddd-lims-db", port="5444", password=pwd)
    
    conn.execute("select p.stable_id,p.decipher_id,clin.* from fe.person p join clinician clin using(id_clinician) where p.proband = true and p.failed = false")
    
    # we need decipher ID, sex, variant location, ref and alt alleles, vep
    # consequence, gene symbol, protein consequence (e.g. "p.Arg222*")referring
    # clinician, validation status. HGVS for transcript and protein,
    #
    
    return md

def get_anthropometric(sample_ids, pheno):
    """ prepare a markdown formatted table of anthropometric data for the probands
    """
    
    samples = pheno[pheno["person_stable_id"].isin(sample_ids)]
    
    table = samples[["patient_id", "last_clinical_assessment",
        "decimal_age_at_assessment", "birthweight", "height", "weight", "ofc"]]
    
    # format the table as markdown
    header = ["ID", "date", "age", "birthweight", "height", "weight", "OFC"]
    table = table.values.tolist()
    tab = tabulate(table, header, tablefmt="pipe")
    
    tab = "\n" + tab + "\n"
    
    return tab

def get_hpo_phenotypes(sample_ids, pheno):
    """ prepare markdown-formatted tables of HPO terms per proband
    """
    
    samples = pheno[pheno["person_stable_id"].isin(sample_ids)]
    
    md = ""
    for (x, row )in samples.iterrows():
        codes = row["child_hpo"].split("|")
        terms = row["child_terms"].split("|")
        
        table = zip(codes, terms)
        tab = tabulate(table, ["HPO code", "term"], tablefmt="pipe")
        child_id = row["person_stable_id"]
        
        md += "\n#### {}\n".format(child_id)
        md += tab + "\n"
    
    return md

phenotypes = open_phenotypes(PHENOTYPES, SANGER_IDS)

sample_ids = ["DDDP111164", "DDDP102163"]

md_text = ""
md_text += get_anthropometric(sample_ids, pheno)
md_text += get_hpo_phenotypes(sample_ids, pheno)

html = markdown.markdown(md_text, ['markdown.extensions.extra'])
html = HTML(string=html)

html.write_pdf('temp.pdf')
