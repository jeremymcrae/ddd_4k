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
import math

import pandas
from statsmodels.stats.multitest import fdrcorrection

from ddd_4k.load_files import open_de_novos, add_decipher_ids, \
    get_significant_results
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
    parser.add_argument("--output-association-table", \
        default="association_table.tsv", \
        help="Path to send output to.")
    
    args = parser.parse_args()
    
    return args

def get_sites_for_validation(de_novos_path, validations_path, new_genes):
    """ find which sites still require validation
    """
    
    # reload the de novos file, but this time retain the sites which failed
    # validation, so that we can count the sites which failed, as well as the
    # sites which passed.
    de_novos = open_de_novos(de_novos_path, validations_path, exclude_invalid=False)
    
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
    
    return (untested, invalid)

def prepare_association_table(table):
    
    # count the number of DDD and external de novos in the different functional
    # categories
    for pop in ["ddd", "meta"]:
        table[pop + ".lof"] = table[pop + ".lof_indel"] + table[pop + ".lof_snv"]
        table[pop + ".missense"] = table[pop + ".missense_indel"] + table[pop + ".missense_snv"]
    
    table["meta.lof"] = table["meta.lof"] - table["ddd.lof"]
    table["meta.missense"] = table["meta.missense"] - table["ddd.missense"]
    
    # check whether the p-value came from the DDD only, or the meta-analysis
    table["Test"] = table["p_min"] == table["ddd.p_min"]
    table["Test"] = table["Test"].map({True: "DDD", False: "Meta"})
    table = table.sort("p_min")
    
    table["Gene"] = table["hgnc"]
    table["P value"] = table["p_min"].map('{:.1e}'.format)
    
    # format the number of DDD and external de novos as "DDD n (external n)"
    meta_nonzero = table["meta.missense"] > 0
    table["Missense"][~meta_nonzero] = table["ddd.missense"][~meta_nonzero].astype(int).apply(str)
    table["Missense"][meta_nonzero] = \
        table["ddd.missense"][meta_nonzero].astype(int).apply(str) + " (" + \
        table["meta.missense"][meta_nonzero].astype(int).apply(str) + ")"
    meta_nonzero = table["meta.lof"] > 0
    table["PTV"][~meta_nonzero] = table["ddd.lof"][~meta_nonzero].astype(int).apply(str)
    table["PTV"][meta_nonzero] = \
        table["ddd.lof"][meta_nonzero].astype(int).apply(str) + " (" + \
        table["meta.lof"][meta_nonzero].astype(int).apply(str) + ")"
    
    # figure out whether the tested subset had clustering
    table["Clustering"] = [ x["ddd.p_missense_clust"] if x["test"] == "DDD" \
        else x["meta.p_missense_clust"] for row, x in table.iterrows() ]
    table["Clustering"][table["Clustering"].isnull()] = 1
    
    # define the gene as having clustered mutations if the clustering p-value is
    # less than 0.01 (the standard alpha for this manuscript)
    table["Clustering"] = [ "Yes" if x < 0.01 else "No" for x in table["Clustering"] ]
    
    formatted = table[["Gene", "Missense", "PTV", "P value", "Test", "Clustering"]]
    
    return formatted

def main():
    args = get_options()
    
    de_novos = open_de_novos(de_novos_path, validations_path)
    de_novos = add_decipher_ids(de_novos, args.sanger_ids)
    
    de_novos = de_novos[["person_stable_id", "decipher_id", "sex", "chrom", "pos",
        "ref", "alt", "symbol", "consequence", "status"]]
    
    table = get_significant_results(args.results, THRESHOLD)
    new_genes = table["hgnc"]
    
    sites = de_novos[de_novos["symbol"].isin(new_genes)]
    sites = sites.sort(["symbol", "pos"])
    
    sites[["person_stable_id", "sex", "chrom", "pos", "ref", "alt", "symbol", \
        "consequence", "status"]].to_csv(args.output, sep="\t", index=False)
    
    formatted = prepare_association_table(table)
    formatted.to_csv(args.output_association_table, sep="\t", index=False)
    
    # find sites for further validation experiments
    (untested, invalid) = get_sites_for_validation(args.de_novos,
        args.validations, de_novos, new_genes)
    untested_path = "ddd_4k_validations.additional_genes.2015-10-12.txt"
    retest_path = "ddd_4k_validations.repeat_validations.2015-10-12.txt"
    untested[["person_stable_id", "chrom", "pos", "ref", "alt", "symbol", \
        "consequence"]].to_csv(untested_path, sep="\t", index=False)
    invalid[["person_stable_id", "chrom", "pos", "ref", "alt", "symbol",  \
        "consequence"]].to_csv(retest_path, sep="\t", index=False)
    

if __name__ == '__main__':
    main()
