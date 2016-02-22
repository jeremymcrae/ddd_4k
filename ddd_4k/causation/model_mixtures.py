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

from __future__ import division, print_function

import numpy
import pandas

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, gridspec
import seaborn

from ddd_4k.causation.classify_known_genes import classify_monoallelic_genes
from ddd_4k.causation.pli_functions import include_constraints, get_constraint_bins
from ddd_4k.causation.plot_pli_bins import plot_by_hi_bin
from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.aggregate_pli_bins import aggregate

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def model_mixing(known, de_novos, expected, constraints):
    """ identify the optimal mixing proportion of haploinsufficient and
    nonhaploinsufficient across to match the observed mixture.
    
    Args:
        known: pandas DataFrame of known developmental disorder genes, including
            columns for hgnc symbols ("gencode_gene_name"), a column for the mode of
            inheritance ("mode"), and a column for the mechanism of action
            ("mech"). We select the dominant genes, and further stratify to
            the haploinsufficient and nonhaploinsufficient dominant genes.
        de_novos: pandas DataFrame of candidate de novos within the cohort.
        expected: counts of de novos expected per gene given our cohort size,
            across different functional categories, for all genes in the genome.
    
    Returns:
        proportion of de novos as haploinsufficient that best approximates the
        difference in observed to expected counts across the pLI bins.
    """
    
    # classify known dominant genes
    mono = classify_monoallelic_genes(known)
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    # identify which pLI quantile each gene falls into
    merged = include_constraints(merged, constraints)
    merged["pLI_bin"] = get_constraint_bins(merged, bins=20)
    
    hi_merged = merged[merged["hgnc"].isin(mono["haploinsufficient"])]
    non_hi_merged = merged[merged["hgnc"].isin(mono["nonhaploinsufficient"])]
    
    target = aggregate(merged, ["missense", "lof"])
    hi_start = aggregate(hi_merged, ["lof", "missense"])
    non_hi_start = aggregate(non_hi_merged, ["missense"])
    
    fig = pyplot.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(3, 1)
    
    all_ax = pyplot.subplot(gs[0])
    hi_ax = pyplot.subplot(gs[1])
    nonhi_ax = pyplot.subplot(gs[2])
    
    plot_by_hi_bin(target, "delta", title="all genes", ax=all_ax)
    plot_by_hi_bin(hi_start, "delta", title="HI genes", ax=hi_ax)
    plot_by_hi_bin(non_hi_start, "delta", title="non-HI genes", ax=nonhi_ax)
    
    fig.savefig("results/obs_to_exp_delta_by_hi_bin.pdf", format="pdf", bbox_inches='tight', pad_inches=0)
    pyplot.close()
    
    return optimise_mixing(hi_start, non_hi_start, target, "results/hi_bin_mixing_optimisation.pdf")

def optimise_mixing(hi_start, non_hi_start, target, output):
    """ identify optimal mixing proportion of HI and non-HI genes to reproduce
    observed frequencies across the pLI bins.
    
    Args:
        hi_start: pandas DataFrame of observed to expected differences across
            the pLI bins, for the known dominant haploinsufficient genes.
        non_hi_start: pandas DataFrame of observed to expected differences across
            the pLI bins, for the known dominant nonhaploinsufficient genes.
        target: pandas DataFrame of observed to expected differences across
            the pLI bins, for all genes with observed candidate de novos.
        output: path to save pdf plot to
    
    Returns:
        proportion of loss-of-function variants required to best capture the
        observed frequencies at pLI bins.
    """
    
    increments = 200.0
    lof_freqs = [ x/increments for x in range(int(increments) + 1) ]
    difference = []
    for lof_frequency in lof_freqs:
        mis_frequency = 1 - lof_frequency
        
        mixed = hi_start["delta"] * lof_frequency + \
            non_hi_start["delta"] * mis_frequency
        
        difference.append(sum((mixed - target["delta"])**2))
    
    mixtures = pandas.DataFrame({"HI_frequency": lof_freqs, "goodness_of_fit": difference})
    plot_optimisation(mixtures, output)
    
    return list(mixtures["HI_frequency"])[numpy.argmin(mixtures["goodness_of_fit"])]

def plot_optimisation(mixtures, output):
    """ plot the goodness of fit for the differences between the observed
    mixture to simulated mixtures of haploinsufficient and nonhaploinsufficient.
    
    Args:
        mixtures: pandas DataFrame with a column for the simulated HI mixture
            frequency (0.0-1.0), and a column for goodness of fit.
        output: path to save a pdf to.
    """
    
    # plot the optimisation efforts
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(mixtures["HI_frequency"], mixtures["goodness_of_fit"],
        marker="None")
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("sum of squares (observed - simulated)")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0)
