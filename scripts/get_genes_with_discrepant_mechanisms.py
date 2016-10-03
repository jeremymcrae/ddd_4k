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

from ddd_4k.constants import KNOWN_GENES, THRESHOLD
from mupit.open_ddd_data import open_known_genes

RESULTS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.all.2015-10-12.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="identify known genes where we"
        "have a different mechanism from what is known for the gene")
    parser.add_argument("--known-genes", default=KNOWN_GENES, \
        help="Path to table of known genes.")
    parser.add_argument("--results", default=RESULTS_PATH, \
        help="Path to table of association results.")
    parser.add_argument("--output", default="missing_mechanism_genes.txt", \
        help="Path to send output to.")
    
    args = parser.parse_args()
    
    return args

def get_lof_genes(results, known):
    """ find the genes where we have evidence of loss-of-function mutations,
    but where the literature suggests mutations shouldn't be loss-of-function.
    
    Args:
        results: pandas dataframe of association results
        known: pandas dataframe of known genes, including their mechanism
    
    Returns:
        pandas dataframe of genes with discrepant mechanisms
    """
    
    # identify the genes in DDG2P that don't have loss-of-function mechanisms
    non_lof = known["gene"][known["mech"] != "Loss of function"]
    non_lof = sorted(non_lof.unique())
    
    # some of the genes have multiple mechanisms, we only want the set without
    # any loss-of-function mechanism
    non_lof = [ x for x in non_lof if "Loss of function" not in list(known["mech"][known["gene"] == x]) ]
    
    # find the genes which have a significant loss-of-function enrichment
    lof_significant = results[(results["meta.p_lof"] < THRESHOLD) | (results["ddd.p_lof"] < THRESHOLD)]
    
    # find the DDG2P genes with a significant loss-of-function enrichment, but
    # without a known loss-of-function mechanism
    return lof_significant[lof_significant["hgnc"].isin(non_lof) ]

def get_gof_genes(results, known):
    """ find the genes where we have evidence of gain-of-function mutations,
    but where the literature suggests mutations shouldn't be gain-of-function.
    
    Args:
        results: pandas dataframe of association results
        known: pandas dataframe of known genes, including their mechanism
    
    Returns:
        pandas dataframe of genes with discrepant mechanisms
    """
    
    # find the known geens with only loss-of-function mechanism
    lof = []
    lof_set = set(["Loss of function"])
    for (symbol, rows) in known.groupby("gene"):
        if set(rows["mech"]) == lof_set:
            lof.append(symbol)
    
    # Now the genes which have gain-of-function mutations, but where the known
    # genes indicates a loss-of-function mechanism. This has a different
    # p-value threshold, since weonly want to know about the subset of genes
    # with multiple missense variants. This might not be a good threshold.
    gof_threshold = 0.05/sum(~results["meta.p_missense_clust"].isnull())
    gof = results[~results["meta.p_missense_clust"].isnull() & ~results["ddd.p_missense_clust"].isnull()]
    gof_significant = gof[(gof["meta.p_missense_clust"] < gof_threshold) | (gof["ddd.p_missense_clust"] < gof_threshold)]
    
    return gof_significant[gof_significant["hgnc"].isin(lof) ]

def main():
    args = get_options()
    
    results = pandas.read_table(args.results, sep="\t")
    
    results = results[~results["p_min"].isnull()]
    
    ddg2p = open_known_genes(args.known_genes)
    
    missing_lof = get_lof_genes(results, ddg2p)
    missing_gof = get_gof_genes(results, ddg2p)
    
    missing_lof["reason"] = "missing loss-of-function mechanism"
    missing_gof["reason"] = "missing gain-of-function mechanism"
    
    missing = missing_lof.append(missing_gof, ignore_index=True)
    
    missing.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
