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
import os

from numpy import log10, median, mean, logical_not
from statsmodels.stats.multitest import fdrcorrection
import pandas

import matplotlib
matplotlib.use("Agg")
import seaborn

from ddd_4k.constants import KNOWN_GENES, THRESHOLD
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
    parser.add_argument("--modified", nargs="+", \
        help="Paths to files containing enrichment results for the modified subsets.")
    parser.add_argument("--known-genes", default=KNOWN_GENES,
        help="Path to table of genes known to be involved in developmental disorders.")
    parser.add_argument("--output", default="temp.pdf",
        help="Path to plot output to.")
    
    args = parser.parse_args()
    
    return args

def load_results(path, known_genes):
    """
    """
    
    dominant_modes = ["Monoallelic", "X-linked dominant"]
    symbols = known_genes["gencode_gene_name"]
    dominant = symbols[known_genes["mode"].isin(dominant_modes)].unique()
    
    table = pandas.read_table(path, sep="\t")
    
    table["dominant"] = table["hgnc"].isin(dominant)
    table["log_p"] = -log10(table[["p_lof", "p_func"]].min(axis=1))
    table = table.dropna()
    
    table = table[["hgnc", "dominant", "log_p"]]
    
    return table

def main():
    
    args = get_options()
    
    known_genes = open_known_genes(args.known_genes)
    
    baseline = load_results(args.baseline, known_genes)
    modified = [ load_results(x, known_genes) for x in args.modified ]
    
    baseline["baseline"] = baseline["log_p"]
    baseline = baseline.drop("log_p", axis=1)
    
    # rename the significance column for the subsets with part of the filename
    modified_names = [ os.path.basename(x).split(".")[3] for x in args.modified ]
    for (x, table) in enumerate(modified):
        name = "modified.{}".format(modified_names[x])
        table[name] = table["log_p"]
    
    modified = [ x.drop("log_p", axis=1) for x in modified ]
    
    # merge all the dataframes from the various modified subsets together
    current = modified[0][["hgnc", "dominant"]]
    for frame in modified:
        current = current.merge(frame, how="outer", on=["hgnc", "dominant"])
    
    merged = baseline.merge(current, how="outer", on=["hgnc", "dominant"])
    merged = merged.dropna()
    
    # Determine the set of genes for comparison. For this, I've decided on the
    # genes that have FDR < 0.05, but which are not genomewide, since these are
    # the genes which might be expected to shift in significance when additional
    # cohorts are included.
    fdr = fdrcorrection(10**-merged["baseline"])[1]
    nearly = merged[(fdr < 0.05) & (merged["baseline"] < -log10(THRESHOLD))]
    
    # figure out the difference in p-values between the modified subsets versus
    # the baseline p-value.
    columns = [ x for x in nearly.columns if "modified" in x ]
    for x in columns:
        name = x.replace("modified.", "")
        nearly[name] = merged[x] - merged["baseline"]
    
    melted = pandas.melt(nearly, id_vars=["hgnc", "dominant"], \
        value_vars=modified_names, value_name="delta")
    
    fig = seaborn.factorplot(x="dominant", y="delta", hue="variable",
        data=melted, kind="violin", size=8, aspect=1.5, legend_out=False, width=1)
    fig.savefig(args.output, format="pdf")

if __name__ == '__main__':
    main()
