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

import pandas
from numpy import median, mean
import matplotlib
matplotlib.use("Agg")

import seaborn

from ddd_4k.constants import VALIDATIONS, DENOVO_PATH
from ddd_4k.load_files import open_de_novos

from mupit.mutation_rates import get_default_rates, get_expected_mutations
from mupit.open_ddd_data import get_ddd_rates

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

RATES_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.mutation_rates.2015-11-24.txt"

def load_neurodevelopmental(path):
    """ load the dominant LoF neurodevelopmental genes, with recognizability
    """
    
    neurodev = pandas.read_excel(path, sheetname="neurodevelopmental.dominant_lof")
    
    return neurodev

def count_de_novos(de_novo_path, validations_path):
    """ count the loss-of-function mutations for each gene.
    """
    
    de_novos = open_de_novos(de_novo_path, validations_path)
    
    counts = {}
    for (gene, rows) in de_novos.groupby("hgnc"):
        counts[gene] = len(rows)
    
    return counts

def plot_recognisability(neurodev):
    """ plot the observed expected ratios at different recognisability indexes
    """
    
    fig = seaborn.factorplot(x="Recognisable", y="ratio", data=neurodev,
        kind="box", order=["1+2", 3, 4, 5], size=6, fliersize=0)
    
    lim = fig.ax.set_ylim((0, 80))
    
    lab = fig.ax.set_ylabel("observed/expected")
    lab = fig.ax.set_xlabel("Clinical recognisability")
    
    fig.savefig("clinical_recognisability.pdf", format="pdf")
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
    
    neurodev = neurodev[~neurodev["ratio"].isnull()].copy()
    
    baseline = median(neurodev["ratio"][neurodev["Recognisable"] != 5])
    recognisable = median(neurodev["ratio"][neurodev["Recognisable"] == 5])
    
    delta_ratio = (baseline - recognisable)
    
    # figure out how many we should have expected in the group of most clinically
    # recognisable genes.
    expected_sum = neurodev["expected"][neurodev["Recognisable"] == 5].sum()
    
    # count how many were missing from those genes, given the difference in
    # obs/exp ratio for the most recognisable genes to the less recognisable
    # genes
    return delta_ratio * expected_sum

def main():
    neurodev_path = "/nfs/users/nfs_j/jm33/neurodevelopmental.dominant_lof_DRF.xlsx"
    
    rates = get_ddd_rates(RATES_PATH)
    
    # TODO: check sex counts
    male = 2407
    female = 1887
    expected = get_expected_mutations(rates, male, female)
    expected["expected"] = expected[["lof_indel", "lof_snv", "missense_indel",
        "missense_snv"]].sum(axis=1)
    
    expected = expected[["hgnc", "expected"]]
    
    counts = count_de_novos(DENOVO_PATH, VALIDATIONS)
    neurodev = load_neurodevelopmental(neurodev_path)
    
    neurodev["observed"] = 0
    for hgnc in counts:
        if hgnc in neurodev["hgnc"].values:
            neurodev["observed"][neurodev["hgnc"] == hgnc] = counts[hgnc]
    
    # plot LoF counts by mutation rates, shade by recognizability
    neurodev = neurodev.merge(expected, how="left", on=["hgnc"])
    
    # calculate the observed/expected ratio for individual genes
    neurodev["ratio"] = neurodev["observed"]/neurodev["expected"]
    
    neurodev["ratio"][neurodev["ratio"].isnull()] = 0
    
    # join the group "1" with group "2", since both groups have fewer genes than
    # the other groups.
    recode = {1: "1+2", 2: "1+2", 3:3, 4:4, 5:5}
    neurodev["Recognisable"] = neurodev["Recognisable"].map(recode)
    
    plot_recognisability(neurodev)
    
    missing = estimate_missing_variants(neurodev)
    print(missing)

if __name__ == '__main__':
    main()
