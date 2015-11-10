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

import argparse
import getpass
import json

import pandas
import psycopg2
from tabulate import tabulate
import markdown
from weasyprint import HTML

from ddd_4k.constants import PHENOTYPES, SANGER_IDS
from ddd_4k.load_files import open_phenotypes

from denovonear.ensembl_requester import EnsemblRequest

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--variants", \
        help="Path to table of variants in novel geens.")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to table of phenotypes.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS, \
        help="Path to table of alternate IDs for participants.")
    parser.add_argument("--output", default="temp.txt", \
        help="Path to send output to.")
    
    args = parser.parse_args()
    
    return args

class EnsemblVariant(EnsemblRequest):
    
    def parse_alternate_consequences(self, variant, gene, vep_cq):
        """ at the moment we're looking for the HGVS nomenclature for the
        canonical transcript, where that transcript needs to match our most
        severe consequence.
        
        TODO: I might need to shift this to using HGVS nomenclature for the
        longest protein_coding transcript, where that transcript consequence
        matches our most severe consrquence prediction.
        
        Args:
            variant: dictionary of possible consequences or variant, obtained
                from the Ensembl REST API.
            gene: HGNC symbol for the variant that we need to match against.
            vep_cq: VEP consequence that we need to match against.
        
        Returns:
            dictionary entry for the canonical transcript for the correct gene
            and VEP consequence. Returns None if not found.
        """
        
        for cq in variant["transcript_consequences"]:
            if "gene_symbol" not in cq:
                continue
            if cq["gene_symbol"] != gene:
                continue
            if cq["biotype"] not in ["protein_coding", "polymorphic_pseudogene"]:
                continue
            if vep_cq not in cq["consequence_terms"]:
                continue
            
            if "canonical" in cq:
                return cq
    
    def get_variant_hgvs(self, chrom, pos, ref, alt, gene, vep_cq):
        """ extract the HGVS predictions for a variant
        """
        
        self.request_attempts = 0
        end_pos = pos + len(ref) - 1
        
        headers = {"Content-Type": "application/json"}
        coordinate = "{}:{}:{}".format(chrom, pos, end_pos)
        ext = "/vep/human/region/{}/{}?hgvs=1;canonical=1".format(coordinate, alt)
        
        r = self.ensembl_request(ext, coordinate, headers)
        variant = json.loads(r)[0]
        
        return self.parse_alternate_consequences(variant, gene, vep_cq)

def get_clinician_details(sample_ids, cur):
    """ get details of clinicians who recruited each proband into the DDD
    
    Args:
        sample_ids: list of DDD sample IDs
        cur: sql cursor, to access the DDD database
    
    Returns:
        dictionary of clinician details, indexed by DDD sample ID.
    """
    
    cur.execute("select p.stable_id,p.decipher_id,clin.* " \
        "from fe.person p join clinician clin using(id_clinician) " \
            "where p.proband = true and p.failed = false and p.stable_id = any(%s);", (sample_ids, ))
    
    query = cur.fetchall()
    
    # reformat the table into a dictionary, so we can pull out the details given
    # a sample ID
    clinicians = {}
    for x in query:
        sample_id = x[0]
        clinicians[sample_id] = {}
        clinicians[sample_id]["decipher_id"] = x[1]
        clinicians[sample_id]["email"] = x[3]
        clinicians[sample_id]["clinician"] = x[4]
    
    return clinicians

def get_variant_details(variants, ensembl, cur):
    """
    we need decipher ID, sex, variant location, ref and alt alleles, vep
    consequence, gene symbol, protein consequence (e.g. "p.Arg222*"), referring
    clinician, validation status. HGVS for transcript and protein,
    
    Args:
        variants: pandas DataFrame of variants, with sample IDs, coordiantes,
            alleles, consequence and validation statuses.
        ensembl: EnsemblVariant object, for pulling HGVS nomenclature out.
        cur: sql cursor, to access the DDD database
    """
    
    sample_ids = list(variants["person_stable_id"])
    clinicians = get_clinician_details(sample_ids, cur)
    
    variants["transcript_hgvs"] = None
    variants["protein_hgvs"] = None
    variants["clinician"] = None
    variants["clinician_email"] = None
    variants["decipher_id"] = None
    for idx, var in variants.iterrows():
        hgvs = ensembl.get_variant_hgvs(var["chrom"], var["pos"], var["ref"],
            var["alt"], var["symbol"], var["consequence"])
        
        if hgvs is not None:
            if "hgvsp" in hgvs:
                variants["protein_hgvs"].loc[idx] = hgvs["hgvsp"]
            if "hgvsc" in hgvs:
                variants["transcript_hgvs"].loc[idx] = hgvs["hgvsc"]
        
        sample_id = var["person_stable_id"]
        variants["clinician"][idx] = clinicians[sample_id]["clinician"]
        variants["clinician_email"][idx] = clinicians[sample_id]["email"]
        variants["decipher_id"][idx] = clinicians[sample_id]["decipher_id"]
    
    return variants

def get_anthropometric(sample_ids, pheno):
    """ prepare a markdown formatted table of anthropometric data for the probands
    
    Args:
        sample_ids: list or pandas Series of DDD person IDs for probands
        pheno: pandas DataFrame of phenotypic data
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
    
    Args:
        sample_ids: list or pandas Series of DDD person IDs for probands
        pheno: pandas DataFrame of phenotypic data
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


def main():
    
    args = get_options()
    variants = pandas.read_table(args.variants, sep="\t")
    
    ensembl = EnsemblVariant(cache_folder="cache", genome_build="grch37")
    
    pwd = getpass.getpass("DDD database password: ")
    with psycopg2.connect(database="ddd_prod", user="ddd_login_ro", \
            host="ddd-lims-db", port="5444", password=pwd) as conn:
        with conn.cursor() as cur:
            var_table = get_variant_details(variants, ensembl, cur)
    
    phenotypes = open_phenotypes(PHENOTYPES, SANGER_IDS)
    
    md_text = ""
    md_text += get_anthropometric(variants["person_id"], pheno)
    md_text += get_hpo_phenotypes(variants["person_id"], pheno)
    
    html = markdown.markdown(md_text, ['markdown.extensions.extra'])
    html = HTML(string=html)
    
    html.write_pdf('temp.pdf')

if __name__ == '__main__':
    main()
