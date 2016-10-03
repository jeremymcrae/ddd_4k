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

import pandas
import numpy
from scipy.interpolate import UnivariateSpline

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, gridspec

from ddd_4k.causation.classify_known_genes import classify_monoallelic_genes
from ddd_4k.causation.pli_functions import include_constraints, get_constraint_bins
from ddd_4k.causation.plot_pli_bins import plot_by_hi_bin
from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.goodness_of_fit import get_goodness_of_fit
from ddd_4k.causation.model_permutation import permute_fits
from ddd_4k.causation.mixture_optimum_variance import variance_around_optimum

from mupit.constants import LOF_CQ, MISSENSE_CQ

import seaborn
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def model_mixing(known, de_novos, expected, constraints, check_modelling=False, check_variance=False):
    """ identify the optimal mixing proportion of haploinsufficient and
    nonhaploinsufficient across to match the observed mixture.
    
    Args:
        known: pandas DataFrame of known developmental disorder genes, including
            columns for hgnc symbols ("gene"), a column for the mode of
            inheritance ("mode"), and a column for the mechanism of action
            ("mech"). We select the dominant genes, and further stratify to
            the haploinsufficient and nonhaploinsufficient dominant genes.
        de_novos: pandas DataFrame of candidate de novos within the cohort.
        expected: counts of de novos expected per gene given our cohort size,
            across different functional categories, for all genes in the genome.
        constraints: pandas DataFrame of constraint scores for all genes in
            the genome with columns for 'gene' and 'pLI'.
        check_modelling: whether to check how well the method can identify
            the correct optimal proportion.
        check_variance: whether to check the variance around the optimal mixing
            proportion by resampling variants according to the optimal ratio.
    
    Returns:
        proportion of de novos as haploinsufficient that best approximates the
        difference in observed to expected counts across the pLI bins.
    """
    
    de_novos['person_id'] = de_novos['person_stable_id']
    
    # classify known dominant genes
    mono = classify_monoallelic_genes(known)
    
    # identify which pLI quantile each gene falls into
    expected['synonymous_expected'] = expected['synonymous_snv']
    expected = include_constraints(expected, constraints)
    bins = [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]
    expected["pLI_bin"] = get_constraint_bins(expected, bins=bins, rate_correct=True)
    # expected["pLI_bin"] = get_constraint_bins(expected, bins=10, rate_correct=True)
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    hi_merged = merged[merged["hgnc"].isin(mono["haploinsufficient"])]
    non_hi_merged = merged[merged["hgnc"].isin(mono["nonhaploinsufficient"])]
    
    missense_excess = aggregate(merged, ["missense"], normalise=True, bins=bins)
    lof_excess = aggregate(hi_merged, ["lof", "missense"], normalise=True, bins=bins)
    gof_excess = aggregate(non_hi_merged, ["missense"], normalise=True, bins=bins)
    
    plot_default_distributions(missense_excess, lof_excess, gof_excess,
        output="results/obs_to_exp_delta_by_hi_bin.pdf")
    
    fits = get_goodness_of_fit(missense_excess, lof_excess, gof_excess, increments=200)
    plot_optimisation(fits, output="results/hi_bin_mixing_optimisation.pdf")
    
    optimal = list(fits["proportion"])[numpy.argmin(fits["goodness_of_fit"])]
    plot_optimisation(fits, output="results/hi_bin_mixing_optimisation.pdf")
    
    # # uncomment the line below if you want to check how a permuted dataset will
    # # behave.
    if check_modelling:
        (slope, intercept) = permute_fits(de_novos, expected, mono, bins, missense_excess, lof_excess, gof_excess, optimal, increments=100, permutations=100)
    
    if check_variance:
        (below, above) = variance_around_optimum(de_novos, expected, mono, optimal, bins, slope, intercept, permutations=1000)
        print('95% CI for proportion from resampling: {0}-{1}'.format(below, above))
    
    return optimal

def plot_default_distributions(missense_excess, lof_excess, gof_excess, output):
    """ plot the unmodified distributions of excess variants by constraint
    quantile for missense DNMs in all genes, DNMs with loss-of-function
    mechanisms, and DNMs with gain-of-function mechanisms.
    
    Args:
        missense_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for all genes with observed candidate de novos.
        lof_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for LoF DNMs in dominant haploinsufficient genes.
        gof_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for missense DNMs in dominant nonhaploinsufficient genes
    """
    
    fig = pyplot.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(3, 1)
    
    all_ax = pyplot.subplot(gs[0])
    hi_ax = pyplot.subplot(gs[1])
    nonhi_ax = pyplot.subplot(gs[2])
    
    plot_by_hi_bin(missense_excess, "delta", title="all genes", ax=all_ax)
    plot_by_hi_bin(lof_excess, "delta", title="HI genes", ax=hi_ax)
    plot_by_hi_bin(gof_excess, "delta", title="non-HI genes", ax=nonhi_ax)
    
    fig.savefig("results/obs_to_exp_delta_by_hi_bin.pdf", format="pdf",
        bbox_inches='tight', pad_inches=0, transparent=True)
    pyplot.close()

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
    
    e = ax.plot(mixtures["proportion"], mixtures["goodness_of_fit"], marker=".",
        linestyle='none')
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("sum of squares (observed - simulated)")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)
