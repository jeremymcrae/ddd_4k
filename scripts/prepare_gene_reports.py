

from tabulate import tabulate
import markdown
from weasyprint import HTML

from ddd_4k.constants import PHENOTYPES, SANGER_IDS
from ddd_4k.load_files import open_phenotypes

def get_variant_details():
    """
    """
    
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
