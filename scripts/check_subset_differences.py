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

from numpy import log10, median, mean, logical_not
import pandas

import matplotlib
matplotlib.use("Agg")
import seaborn

from ddd_4k.constants import KNOWN_GENES
from ddd_4k.load_files import open_known_genes

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check the impact"
        "of including different subsets on the ability to detect genes involved"
        "in developmental disorders.")
    parser.add_argument("--baseline", \
        help="Path to file containing enrichment results for the baseline.")
    parser.add_argument("--modified", \
        help="Path to file containing enrichment results for the modified subset.")
    parser.add_argument("--known-genes", default=KNOWN_GENES,
        help="Path to table of genes known to be involved in developmental disorders.")
    parser.add_argument("--output", default="temp.pdf",
        help="Path to plot output to.")
    
    args = parser.parse_args()
    
    return args

def load_results(path, known_genes):
    """
    """
    
    dominant = known_genes["gencode_gene_name"][known_genes["mode"].isin(["Monoallelic", "X-linked dominant"])].unique()
    
    table = pandas.read_table(path, sep="\t")
    
    table["dominant"] = table["hgnc"].isin(dominant)
    table["log_p"] = -log10(table["enrichment_p_value"])
    table = table.dropna()
    
    return table


def main():
    
    args = get_options()
    
    known_genes = open_known_genes(args.known_genes)
    
    baseline = load_results(args.baseline, known_genes)
    modified = load_results(args.modified, known_genes)
    
    baseline["baseline"] = baseline["log_p"]
    modified["modified"] = modified["log_p"]
    
    merged = baseline[["hgnc", "dominant", "baseline"]].merge(modified[["hgnc", "dominant", "modified"]], how="outer", \
        on=["hgnc", "dominant"])
    
    merged = merged[merged["baseline"] > -log10(0.05/18500)]
    merged = merged.dropna()
    
    merged["delta"] = merged["modified"] - merged["baseline"]
    print(median(merged["delta"][merged["dominant"]]))
    print(median(merged["delta"][merged["dominant"] == False]))
    
    fig = seaborn.factorplot(x="dominant", y="delta", data=merged, kind="violin")
    fig.savefig(args.output, format="pdf")
    
    
    
if __name__ == '__main__':
    main()
