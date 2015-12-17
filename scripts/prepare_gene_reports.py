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

from ddd_4k.constants import PHENOTYPES, SANGER_IDS, DIAGNOSED
from ddd_4k.load_files import open_phenotypes
from ddd_4k.clinical_phenotypes import get_clinical_details, get_variant_details
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
    parser.add_argument("--diagnosed", default=DIAGNOSED, \
        help="Path to table of participants with diagnoses.")
    parser.add_argument("--output-dir", default="gene_reports", \
        help="Folder to send output tables to.")
    
    args = parser.parse_args()
    
    return args

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
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    
    # remove the vairnats in diagnosed individuals, since they didn't contribute
    # to identify the gene as being genomewide significant. These "diagnosed"
    # individuals still need to be inspected, to make sure they don't conflict
    # with the evidence from the other probands.
    variants = variants[~variants["person_stable_id"].isin(diagnosed["person_id"])]
    
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
