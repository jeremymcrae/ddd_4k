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

import pandas
from statsmodels.stats.multitest import fdrcorrection

from ddd_4k.constants import DENOVO_PATH, VALIDATIONS, SANGER_IDS, THRESHOLD

RESULTS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.without_diagnosed.all.2015-11-24.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--de-novos", default=DENOVO_PATH, \
        help="Path to table of de novo variants.")
    parser.add_argument("--validations", default=VALIDATIONS, \
        help="Path to table of de novo validation results.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS, \
        help="Path to table of tbale mapping sanger IDs to DDD IDs.")
    parser.add_argument("--results", default=RESULTS_PATH, \
        help="Path to table of association results.")
    parser.add_argument("--output", default="novel_gene_variants.txt", \
        help="Path to send output to.")
    
    args = parser.parse_args()
    
    return args

def open_de_novos(de_novos_path, validations_path, sanger_ids_path):
    """ load the de novos, and their validation data
    """
    
    validations = pandas.read_table(validations_path, sep="\t")
    de_novos = pandas.read_table(de_novos_path, sep="\t")
    de_novos = de_novos.merge(validations, how="left",
        left_on=["person_stable_id", "chrom", "pos", "ref", "alt", "symbol", "consequence"],
        right_on=["person_id", "chrom", "start_pos", "ref_allele", "alt_allele", "hgnc", "consequence"])
    
    
    sanger_ids = pandas.read_table(sanger_ids_path, sep="\t")
    sanger_ids = sanger_ids[~sanger_ids["person_stable_id"].duplicated() ]
    de_novos = de_novos.merge(sanger_ids, how="left", on="person_stable_id")
    
    de_novos = de_novos[["person_stable_id", "decipher_id", "sex", "chrom", "pos",
        "ref", "alt", "symbol", "consequence", "status"]]
    
    return de_novos

def get_sites_for_validation(de_novos, new_genes):
    """ find which sites still require validation
    """
    
    validations = pandas.DataFrame({"hgnc": new_genes, "tested": None, "untested": None, "invalid": None})
    for pos in validations.index:
        # find the sites in each significant gene
        valids = de_novos[de_novos["symbol"] == validations["hgnc"][pos]]
        validations["tested"][pos] = sum(~valids["status"].isnull())
        validations["invalid"][pos] = sum(valids["status"].isin(["inherited", "false_positive"]))
        validations["untested"][pos] = sum(valids["status"].isnull())
    
    # figure out which sites in the strongly associated genes have not yet been
    # tested and double check the sites in the strongly associated genes that
    # failed validation
    untested = de_novos[de_novos["symbol"].isin(new_genes) & de_novos["status"].isnull() ]
    invalid = de_novos[de_novos["symbol"].isin(new_genes) & de_novos["status"].isin(["inherited", "false_positive"]) ]
    
    untested[["person_stable_id", "chrom", "pos", "ref", "alt", "symbol", \
        "consequence"]].to_csv("ddd_4k_validations.additional_genes.2015-10-12.txt", sep="\t", index=False)
    invalid[["person_stable_id", "chrom", "pos", "ref", "alt", "symbol",  \
        "consequence"]].to_csv("ddd_4k_validations.repeat_validations.2015-10-12.txt", sep="\t", index=False)

def main():
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations, args.sanger_ids)
    de_novos = de_novos[~de_novos["consequence"].isin(["synonymous_variant"]) ]
    
    results = pandas.read_table(args.results, sep="\t")
    
    # find the genes that exceed a multiple testing corrected genonmewide threshold
    new_genes = results["hgnc"][results["p_min"] < THRESHOLD]
    new_genes = new_genes[~new_genes.isnull()]
    
    get_sites_for_validation(de_novos, new_genes)
    
    sites = de_novos[de_novos["symbol"].isin(new_genes) & ~de_novos["status"].isin(["inherited", "false_positive"]) ]
    sites = sites.sort(["symbol", "pos"])
    
    sites[["person_stable_id", "sex", "chrom", "pos", "ref", "alt", "symbol", \
        "consequence", "status"]].to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
