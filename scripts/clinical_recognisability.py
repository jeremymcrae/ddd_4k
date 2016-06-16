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
import sys
import tempfile
import math

import pandas
from numpy import median, mean, sqrt, log10
from scipy.stats import norm, poisson
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

import seaborn

from ddd_4k.constants import VALIDATIONS, DENOVO_PATH
from ddd_4k.load_files import open_de_novos

from mupit.mutation_rates import get_default_rates, get_expected_mutations
from mupit.open_ddd_data import get_ddd_rates

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

RATES_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.mutation_rates.2015-11-24.txt"
NEURODEV_PATH = "/nfs/users/nfs_j/jm33/neurodevelopmental.dominant_lof_DRF.xlsx"
RESULTS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.all.2015-11-24.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to examine differences " \
        "in observed numbers of de novo mutations to expected numbers, across " \
        "different levels of clinical recognisability. This is done for gene s" \
        "which are loss-of-function neurodevelopmental, and which have been " \
        "scored for recognisability by a clinical geneticist.")
    parser.add_argument("--rates", default=RATES_PATH,
        help="Path to table of consequence specific mutation rates per gene.")
    parser.add_argument("--de-novos", default=DENOVO_PATH,
        help="Path to table of candidate de novo mutations.")
    parser.add_argument("--validations", default=VALIDATIONS,
        help="Path to table of validation data.")
    parser.add_argument("--neurodevelopmental-genes", default=NEURODEV_PATH,
        help="Path to table of neurodevelopmental genes with loss-of-function mechanisms.")
    parser.add_argument("--results", default=RESULTS_PATH,
        help="Path to table of results from testing for significance.")
    parser.add_argument("--output",
        default="results/clinical_recognisability.pdf",
        help="Path to plot graph to.")
    
    args = parser.parse_args()
    
    return args

def load_neurodevelopmental(path):
    """ load the dominant LoF neurodevelopmental genes, with recognizability
    """
    
    neurodev = pandas.read_excel(path, sheetname="neurodevelopmental.dominant_lof")
    
    return neurodev

def count_de_novos(de_novo_path, validations_path, lof_only=False):
    """ count the loss-of-function mutations for each gene.
    """
    
    de_novos = open_de_novos(de_novo_path, validations_path)
    
    if lof_only:
        de_novos = de_novos[de_novos["category"] == "truncating"]
    
    counts = {}
    for (gene, rows) in de_novos.groupby("hgnc"):
        counts[gene] = len(rows)
    
    return counts

def plot_recognisability(neurodev, output):
    """ plot the observed expected ratios at different recognisability indexes
    """
    
    neurodev["ratio"] = neurodev["observed"]/neurodev["expected"]
    
    # some genes lack an expected count, since we don't have mutation rates for
    # genes which lack any mutations in them. Their ratio must be 0.
    neurodev["ratio"][neurodev["ratio"].isnull()] = 0
    
    fig = seaborn.factorplot(x="Recognisable", y="ratio", data=neurodev,
        kind="box", order=["1+2", 3, 4, 5], size=6, fliersize=0)
    
    lim = fig.ax.set_ylim((0, 80))
    
    lab = fig.ax.set_ylabel("observed/expected")
    lab = fig.ax.set_xlabel("Clinical recognisability")
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0)
    matplotlib.pyplot.close()

def estimate_missing_variants(neurodev):
    """ estimate the number of missing neurodevelopmental genes
    
    Args:
        neurodev: pandas dataframe of neurodevelopmental genes, with their
            clinical recognisability.
        baseline: baseline proportion of observed to expected mutations at
            lowest clinical recognisability.
        slope: slope of linear regression, so that each increment of
            recognisability reduces the baseline proportion by this much.
        
    Returns:
        number of missing genes
    """
    
    # concatenate the genes in ech recognisability category.
    joined = pandas.pivot_table(neurodev, rows="Recognisable",
        values=["observed", "expected"], aggfunc=sum)
    joined["ratio"] = joined["observed"]/joined["expected"]
    joined["Recognisable"] = joined.index
    
    baseline = joined[joined["Recognisable"] != "5"]
    baseline = sum(baseline["ratio"])/len(baseline)
    
    different = joined[joined["Recognisable"] == "5"]
    delta = baseline * different["expected"] - different["observed"]
    
    return delta[0]

def plot_concatenated_recognisability(neurodev, output):
    """
    """
    
    # concatenate the genes in ech recognisability category.
    joined = pandas.pivot_table(neurodev, rows="Recognisable",
        values=["observed", "expected"], aggfunc=sum)
    joined["ratio"] = joined["observed"]/joined["expected"]
    
    # get the upper and lower confidence intervals for the observed counts
    lower, upper = poisson.interval(0.95, joined["observed"])
    joined["lower"] = abs(lower/joined["expected"] - joined["ratio"])
    joined["upper"] = abs(upper/joined["expected"] - joined["ratio"])
    
    joined["Recognisable"] = joined.index
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.bar(range(len(joined)), joined["ratio"], align="center",
        yerr=[joined["lower"], joined["upper"]],
        ecolor="black", capsize=10, error_kw={'capthick': 2})
    
    # fix the axis limits and ticks
    e = ax.set_xlim((-0.5, len(joined) - 0.5))
    e = ax.set_xticks(range(len(joined)))
    e = ax.set_xticklabels(joined["Recognisable"])
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    e = ax.set_xlabel("recognisability class")
    e = ax.set_ylabel("observed/expected")
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)
    pyplot.close()

def plot_rate_by_p_value(neurodev, results):
    """ plot the mutation rate by p-value for neurodevelopmental genes. Shade by
    clinical recognisability.
    """
    
    results = dict(zip(results["hgnc"], results["p_min"]))
    neurodev["p_value"] = -log10(neurodev["hgnc"].map(results))
    neurodev = neurodev[~neurodev["p_value"].isnull()]
    
    # data = neurodev[["expected", "p_value", "Recognisable"]].copy()
    
    colors = ["green", "blue", "red", "black"]
    groups = ["1+2", "3", "4", "5"]
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    for group, color in zip(groups, colors):
        data = neurodev[neurodev["Recognisable"] == group]
        e = ax.plot(data["expected"], data["p_value"], linestyle='None',
            color=color, marker=".", markersize=10, label=group)
    
    e = ax.legend(fontsize="small")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    e = ax.set_xlabel("Expected protein-truncating mutation rate")
    e = ax.set_ylabel("-log10(P)")
    
    fig.savefig("results/rate_by_p_value.pdf", format="pdf", bbox_inches='tight',
        pad_inches=0, transparent=True)
    pyplot.close()

def get_enrichment_factor(neurodev):
    """ get the enrichment factor for neurodevelopmental genes with low recognizability
    
    Args:
        neurodev: pandas DataFrame of counts per gene, containing columns for
            clinical recognizability ("Recognisable"), the
    
    Returns:
        float value for enrichment factor.
    """
    
    # get the known dominant haploinsufficient genes with low clinical
    # recognisability
    low_recog = neurodev[neurodev["Recognisable"] != "5"]
    
    return sum(low_recog["observed"])/sum(low_recog["expected"])

def main():
    args = get_options()
    
    # rates = get_ddd_rates(args.rates)
    rates = get_default_rates()
    
    results = pandas.read_table(args.results)
    
    # determine the number of mutations we expect per gene, given consequence
    # specific mutation rates for each gene.
    expected = get_expected_mutations(rates, male=2407, female=1887)
    # columns = ["lof_indel", "lof_snv", "missense_indel", "missense_snv"]
    columns = ["lof_indel", "lof_snv"]
    expected["expected"] = expected[columns].sum(axis=1)
    
    lof_only = False
    if columns == ["lof_indel", "lof_snv"]:
        lof_only = True
    
    counts = count_de_novos(args.de_novos, args.validations, lof_only)
    neurodev = load_neurodevelopmental(args.neurodevelopmental_genes)
    
    neurodev["observed"] = 0
    for hgnc in counts:
        if hgnc in neurodev["hgnc"].values:
            neurodev["observed"][neurodev["hgnc"] == hgnc] = counts[hgnc]
    
    # calculate the observed/expected ratio for individual genes
    expected = dict(zip(expected["hgnc"], expected["expected"]))
    neurodev["expected"] = neurodev["hgnc"].map(expected)
    neurodev = neurodev[~neurodev["expected"].isnull()]
    
    # Join the two least recognisable groups, since they are smaller than the
    # other groups.
    recode = {1: "1+2", 2: "1+2", 3: "3", 4: "4", 5: "5"}
    neurodev["Recognisable"] = neurodev["Recognisable"].map(recode)
    
    # plot_recognisability(neurodev, args.output)
    plot_concatenated_recognisability(neurodev, args.output)
    
    missing = estimate_missing_variants(neurodev)
    print("missing lof variants: {}".format(missing))
    
    enrich_factor = get_enrichment_factor(neurodev)
    print("enrichment: {}".format(enrich_factor))
    
    plot_rate_by_p_value(neurodev, results)

if __name__ == '__main__':
    main()
