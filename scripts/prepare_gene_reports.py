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
import sys
import os

import pandas
import psycopg2

from ddd_4k.constants import PHENOTYPES, SANGER_IDS
from ddd_4k.load_files import open_phenotypes
from ddd_4k.ensembl_variant import EnsemblVariant

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--de-novos", \
        help="Path to table of variants in novel genes.")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to table of phenotypes.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS, \
        help="Path to table of alternate IDs for participants.")
    parser.add_argument("--output-dir", default="gene_reports", \
        help="Folder to send output tables to.")
    
    args = parser.parse_args()
    
    return args

def get_clinician_details(sample_ids, cur):
    """ get details of clinicians who recruited each proband into the DDD
    
    Args:
        sample_ids: list of DDD sample IDs
        cur: sql cursor, to access the DDD database
    
    Returns:
        dictionary of clinician details, indexed by DDD sample ID.
    """
    
    cur.execute("select p.stable_id,p.decipher_id,clin.email,clin.name " \
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
        clinicians[sample_id]["email"] = x[2]
        clinicians[sample_id]["clinician"] = x[3]
    
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
    
    columns = ["decipher_id", "sex", "chrom", "pos", "ref", "alt",
        "consequence", "symbol", "status", "transcript_hgvs", "protein_hgvs",
        "clinician", "clinician_email"]
    table = variants[columns]
    
    return table

def get_clinical_details(variants, pheno):
    """ prepare a table of anthropometric data for the probands
    
    Args:
        sample_ids: list or pandas Series of DDD person IDs for probands
        pheno: pandas DataFrame of phenotypic data
        ensembl: EnsemblVariant object, for pulling HGVS nomenclature out.
    """
    
    sample_ids = variants["decipher_id"]
    
    pheno["decipher_id"] = pheno["patient_id"].astype(str)
    
    samples = pheno[pheno["decipher_id"].isin(sample_ids)]
    samples = samples.merge(variants, how="left", on="decipher_id")
    
    columns = ["patient_id", "gender", "chrom", "pos", "ref", "alt",
        "consequence", "transcript_hgvs", "protein_hgvs", "symbol",
        "mothers_age", "fathers_age", "birthweight_sd",
        "gestation", "birth_ofc_sd", "decimal_age_at_assessment", "height_sd",
        "weight_sd", "ofc_sd", "sat_independently", "walked_independently",
        "first_words", "child_terms",  "syndrome", "additional_comments"]
    table = samples[columns]
    
    return table

def get_tables(ensembl, cur, variants, phenotypes):
    """ get tables of variant and clinical details for variants/probands for a gene
    
    Args:
        ensembl: EnsemblVariant object, for pulling HGVS nomenclature out.
        cur: sql cursor, to access the DDD database
        variants: pandas DataFrame of variants foir a single gene, with
            sample IDs, coordinates, alleles, consequence and validation statuses.
        phenotypes: pandas DataFrame of phenotypic data
    
    Returns:
        tuple of (variant table, clinical table) dataframes
    """
    
    var_table = get_variant_details(variants, ensembl, cur)
    
    clinical_table = get_clinical_details(var_table, phenotypes)
    
    return (var_table, clinical_table)

def main():
    
    args = get_options()
    variants = pandas.read_table(args.de_novos, sep="\t")
    phenotypes = open_phenotypes(PHENOTYPES, SANGER_IDS)
    
    ensembl = EnsemblVariant(cache_folder="cache", genome_build="grch37")
    
    tries = 0
    while tries < 3:
        try:
            pwd = getpass.getpass("DDD database password: ")
            conn = psycopg2.connect(database="ddd_prod", user="ddd_login_ro", \
                    host="ddd-lims-db", port="5444", password=pwd)
            break
        except psycopg2.OperationalError:
            tries += 1
    
    if tries == 3:
        sys.exit("too many wrong password attempts")
    
    with conn.cursor() as cur:
        
        var_table, clinical_table = get_tables(ensembl, cur, variants, phenotypes)
        
        var_path = os.path.join(args.output_dir, "ddd_4k_variants.txt")
        clin_path = os.path.join(args.output_dir, "ddd_4k_clinical.txt")
        
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        
        var_table.to_csv(var_path, sep="\t", na_rep="NA", index=False)
        clinical_table.to_csv(clin_path, sep="\t", na_rep="NA", index=False)
    
    conn.close()

if __name__ == '__main__':
    main()
