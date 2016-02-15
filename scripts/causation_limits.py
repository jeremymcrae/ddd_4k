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

from __future__ import division

import numpy
import pandas

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.load_files import open_de_novos, open_known_genes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, CONSTRAINTS_URL
from ddd_4k.count_mutations_per_person import get_count_by_person

from mupit.mutation_rates import get_default_rates, get_expected_mutations
from mupit.count_de_novos import get_de_novo_counts
from mupit.gene_enrichment import gene_enrichment

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def classify_monoallelic_genes(known):
    """ classify monoalleleic genes into haploinsufficient and
        nonhaploinsufficient sets
    
    Args:
        known: pandas DataFrame of known developmental dfisorder genes
    
    Returns:
        dictionary of gene sets, one entry for haploinsufficient genes, one
        entry for nonhaploinsufficient genes
    """
    
    # identify known haploinsufficient dominant genes
    monoallelic = known[known["mode"].isin(["Monoallelic", "X-linked dominant"])]
    haploinsufficient = set(monoallelic["gencode_gene_name"][
        monoallelic["mech"].isin(["Loss of function"])])
    
    # identify known nonhaploinsufficient dominant genes
    nonhaploinsufficient = set(monoallelic["gencode_gene_name"][
        monoallelic["mech"].isin(["Activating", "Dominant negative"])])
    
    # remove genes which fall into both categories, since de novos in those will
    # be a mixture of the two models, and we need to cleanly separate the
    # underlying models.
    overlap = haploinsufficient & nonhaploinsufficient
    haploinsufficient -= overlap
    nonhaploinsufficient -= overlap
    
    return {"haploinsufficient": haploinsufficient,
        "nonhaploinsufficient": nonhaploinsufficient}

def count_synonymous_per_gene(de_novos):
    """ count the number of synonymous mutations per gene
    
    Args:
        de_novos: pandas DataFrame containing one row per candidate de novo
            mutation.
    
    Returns:
        pandas DataFrame of counts of synonymous variants per gene
    """
    
    synonymous = de_novos[de_novos["consequence"] == "synonymous_variant"]
    counts = synonymous.pivot_table(rows="symbol", values="alt", aggfunc=len)
    counts = pandas.DataFrame({"hgnc": counts.index, "observed": counts})
    
    counts.index = range(len(counts))
    
    return counts

def get_pp_dnm_threshold(de_novos, expected):
    """ get the pp_dnm threshold where we expect as many synonymous de novos as we observed
    
    Args:
        de_novos: pandas DataFrame of all candidate coding de novos, including
            synonymous variants
        rates: pandas DataFrame of numbers of expected mutations within the
            cohort, acorss different categories. This includes a "synonymous_snv"
            column. The table includes all genes in the genome (or as close as
            possible), since the number of expected
    
    Returns:
        pp_dnm value, which if we exclude candidates below this value, then the
        ratio of observed synonymous variants to expected equals one. Applying
        the threshold to all candidate de novos threshold ensures we do not have
        more candidates due to errors than expected by chance in missense and
        loss-of-function categories.
    """
    
    rates = dict(zip(expected["hgnc"], expected["synonymous_snv"]))
    
    epsilon = 0.0009
    low = 0
    high = 1
    while True:
        # define a mid point to bisect at
        mid = (low + high)/2
        # select the candidates above the midpoint threshold
        candidates = de_novos[de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > mid) ]
        synonymous = count_synonymous_per_gene(candidates)
        
        # make sure we have expected numbers available. Some genes lack
        # expected numbers, so we exclude those genes, rather than counting
        # extra genes, for which we can never know if we are near expectations.
        synonymous["expected"] = synonymous["hgnc"].map(rates)
        synonymous = synonymous[~synonymous["expected"].isnull()]
        
        # we divide by the total expected number of synonymous variants across
        # all genes in the genome, rather than just the geens we have observed a
        # syonymous variant in.
        ratio = sum(synonymous["observed"])/sum(expected["synonymous_snv"])
        
        if ratio > 1:
            low = mid
        elif ratio < 1:
            high = mid
        
        if abs(ratio - 1) < epsilon:
            break
    
    return mid

def merge_observed_and_expected(de_novos, expected):
    """ get per gene counts of observed and expected de novo mutations
    
    Args:
        expected: pandas DataFrame of counts of expected de novo mutations for
            all genes (or nearly all) in the genome.
        de_novos: pandas DataFrame of candidate de novo mutations observed in
            the cohort.
    
    Returns:
        pandas DataFrame of rows per gene of observed and expected counts for
        loss-of-function and missense candidates.
    """
    
    # sum the observed loss-of-function and missense candidates per gene
    observed = get_de_novo_counts(de_novos)
    observed["lof_observed"] = observed[["lof_snv", "lof_indel"]].sum(axis=1)
    observed["missense_observed"] = observed[["missense_snv", "missense_indel"]].sum(axis=1)
    
    # restrict the observed counts to the columns for merging.
    observed = observed[["hgnc", "lof_observed", "missense_observed"]]
    
    # merge the observed counts with the expected counts, for genes without
    # counts to merge, replace NA values with 0.
    expected = expected[["hgnc", "lof_indel", "lof_snv", "missense_indel",
        "missense_snv", "synonymous_snv"]].copy()
    merged = expected.merge(observed, how="left", on="hgnc")
    merged["lof_observed"][merged["lof_observed"].isnull()] = 0
    merged["missense_observed"][merged["missense_observed"].isnull()] = 0
    
    # get the numbers of expected loss-of-function, missense and
    # synonymous mutations.
    merged["lof_expected"] = merged[["lof_snv", "lof_indel"]].sum(axis=1)
    merged["missense_expected"] = merged[["missense_snv", "missense_indel"]].sum(axis=1)
    merged = merged.rename(columns={'synonymous_snv': 'synonymous_expected'})
    
    # restrict the table to the minimum required columns
    merged = merged[["hgnc", "lof_observed", "missense_observed",
        "lof_expected", "missense_expected", "synonymous_expected"]].copy()
    
    return merged

def include_constraints(table, constraints):
    """ add a column for pLI constraints scores
    
    Args:
        table: pandas DataFrame, with one gene per row
        constraints: table of constraints scores, including columns for HGNC
            symbols ("gene") and constraints score ("pLI").
    
    Returns:
        same pandas DataFrame, with an extra pLI column (and minus genes without
        pLI scores.)
    """
    
    # load the pLI data, and append column to the table
    recode = dict(zip(constraints["gene"], constraints["pLI"]))
    table["pLI"] = table["hgnc"].map(recode)
    
    # drop genes which lacka constraint score
    table = table[~table["pLI"].isnull()]
    
    return table

def get_constraint_bins(table, bins=20, rate_correct=False):
    """ identify the quantile bins each gene falls into.
    
    Args:
        table: pandas DataFrame, with one gene per row.
        bins: number of quantiles to divide into.
        rate_correct: whether to adjust the bins to account for different
            mutation rates of genes depending on the pLI score. If used, each
            bin will expect equal numbers of synonymous mutations.
    
    Returns:
        pandas Series, indicating the bin each gene falls into.
    """
    
    # identify which pLI quantile each gene falls into
    quantiles = [ x/float(bins) for x in range(bins + 1) ]
    
    if not rate_correct:
        pLI_bin, bins = pandas.qcut(table["pLI"], q=quantiles,
            labels=quantiles[:-1], retbins=True)
    else:
        summed_synonymous = sum(table["synonymous_expected"])
        
        pLI_bin = []
        current_bin = 0
        current_sum = 0
        pli_sorted = table.sort("pLI")
        for key, row in pli_sorted.iterrows():
            if current_sum > summed_synonymous/float(bins):
                current_bin += 1
                current_sum = 0
            
            current_sum += row["synonymous_expected"]
            pLI_bin.append(quantiles[current_bin])
        
        pli_sorted["pLI_bin"] = pLI_bin
        pli_sorted = pli_sorted.sort_index()
        pLI_bin = pli_sorted["pLI_bin"]
    
    return pLI_bin

def get_overall_consequence_ratios(expected, de_novos):
    """ determine the ratio of observed to expected for synonymous, missense
    and loss-of-function canddiate de novos
    
    Args:
        expected: pandas DataFrame of counts of expected de novo mutations for
            all genes (or nearly all) in the genome.
        de_novos: pandas DataFrame of candidate de novo mutations observed in
            the cohort.
    
    Returns:
        dictionary of counts and ratios for synonymous, missense and
        loss-of-function consequences.
    """
    
    synonymous_count = len(de_novos[de_novos["consequence"] == "synonymous_variant"])
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    lof_ratio = sum(merged["lof_observed"])/sum(merged["lof_expected"])
    missense_ratio = sum(merged["missense_observed"])/sum(merged["missense_expected"])
    
    return {"synonymous": {"ratio": 1.0, "count": synonymous_count},
        "loss-of-function": {"ratio": lof_ratio, "count": sum(merged["lof_observed"])},
        "missense": {"ratio": missense_ratio, "count": sum(merged["missense_observed"])}}

def plot_overall_ratios(ratios, output):
    """ plot the ratio of observed to expected de novo counts for overall consequences
    
    Args:
        ratios: dictionary of de novo counts and observed/expected ratios for
            "synonymous", "loss-of-function" and "missense" candidate de novos.
        output: path to save plot as pdf to.
    """
    
    # format the ratios and counts
    temp = pandas.DataFrame(ratios).transpose()
    temp["consequence"] = temp.index
    temp = temp.reindex(["synonymous", "missense", "loss-of-function"])
    temp.index = range(len(temp))
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.bar(range(len(temp)), temp["ratio"], align="center")
    
    # annotate the plot, to show the baseline, and the numbers of candidate de
    # novos in each category
    e = ax.axhline(1.0, color="black", linestyle="dashed")
    for key, row in temp.iterrows():
        excess = row["count"] * ((row["ratio"] - 1) / row["ratio"])
        e = ax.text(key, row["ratio"]+0.01,
            "n={0:.0f}\nexcess={1:.0f}".format(row["count"], excess),
            horizontalalignment='center')
    
    # fix the axis limits and ticks
    e = ax.set_xlim((-0.5, len(temp) - 0.5))
    e = ax.set_xticks(range(len(temp)))
    e = ax.set_xticklabels(temp["consequence"])
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    e = ax.set_xlabel("consequence class")
    e = ax.set_ylabel("observed/expected")
    
    fig.savefig(output, format="pdf")

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
    
    monoallelic = known[known["mode"].isin(["Monoallelic", "X-linked dominant"])]
    
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
    
    plot_by_hi_bin(target, "delta", "results/obs_to_exp_delta_by_hi_bin.all_genes.pdf")
    plot_by_hi_bin(hi_start, "delta", "results/obs_to_exp_delta_by_hi_bin.hi_genes.pdf")
    plot_by_hi_bin(non_hi_start, "delta", "results/obs_to_exp_delta_by_hi_bin.non_hi_genes.pdf")
    
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

def aggregate(merged, consequences=None, normalise=True):
    """ get the difference between observed and expected counts across pLI bins.
    
    This returns the difference in the number of observed de novo mutations to
    expected numbers of de novo mutations. The difference is given as a delta
    and as a ratio.
    
    Args:
        merged: pandas DataFrame containing one row per gene, which includes
            columns for observed and expected counts, by lof and missense
            categories. This also includes a column indicating which pLI
            vigicile each genes falls into.
        consequences: list of consequence types to include, such as ["missense"],
            or ["missense", "lof"] or ["lof"]
        normalise: boolean for whether to normalise the differences so that the
            total sums to one.
    
    Returns:
        pandas DataFrame with columns for pLI bin, and observed to expected
        differences as delta and ratio.
    """
    
    if consequences is None:
        consequences = ["lof", "missense"]
        
    observed_columns = [ x + "_observed" for x in consequences ]
    expected_columns = [ x + "_expected" for x in consequences ]
    
    aggregated = merged.pivot_table(rows="pLI_bin",
        values=observed_columns + expected_columns, aggfunc=sum)
    
    aggregated["pLI_bin"] = aggregated.index
    
    aggregated["observed"] = aggregated[observed_columns].sum(axis=1)
    aggregated["expected"] = aggregated[expected_columns].sum(axis=1)
    
    aggregated["delta"] = aggregated["observed"] - aggregated["expected"]
    aggregated["ratio"] = aggregated["observed"]/aggregated["expected"]
    
    if normalise:
        aggregated["delta"] = aggregated["delta"]/sum(aggregated["delta"])
        aggregated["ratio"] = aggregated["ratio"]/sum(aggregated["ratio"])
    
    # make sure we have all the quantile bins, even for tables that might lack
    # genes in a given pLI bin
    for x in [ x/20.0 for x in range(20) ]:
        if x not in list(aggregated["pLI_bin"]):
            row = pandas.DataFrame({"pLI_bin": [x], "delta": [0], "ratio": [0]})
            aggregated = aggregated.append(row, ignore_index=True)
    
    aggregated = aggregated[["pLI_bin", "observed", "expected", "delta",
        "ratio"]].copy()
    aggregated = aggregated.sort("pLI_bin")
    aggregated.index = range(len(aggregated))
    
    return aggregated

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
    
    fig.savefig(output, format="pdf")

def plot_by_hi_bin(aggregated, diff_type, output, expected=None, count_halves=False):
    """ plot difference in observed to expected across pLI vigiciles
    
    Args:
        aggregated: pandas DataFrame of observed to expected differences across
            the pLI bins.
        diff_type: list of consequence types to include, for example
            ["missense"] or ["missense", "lof"] or ["lof"]
        output: path to save the plot to.
        expected: position to plot horizontal line of expectation at, for
            example, at a ratio of 1.0.
        count_halves: whether to count the number of observed de novos in each
            half of the plot.
    """
    
    fig = pyplot.figure(figsize=(12, 6))
    ax = fig.gca()
    
    barwidth = 0.04
    e = ax.bar(aggregated["pLI_bin"], aggregated[diff_type], width=barwidth,
        align="center")
    
    # fix the axis limits and ticks
    quantiles = [ x/20.0 for x in range(20) ]
    e = ax.set_xlim((-0.05, max(quantiles) + 0.05))
    e = ax.set_xticks(quantiles)
    e = ax.set_xticklabels(quantiles, fontsize="medium")
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    
    if expected is not None:
        e = ax.axhline(expected, linestyle="dashed", color="black")
    
    if count_halves:
        e = ax.axvline(0.475, linestyle="dashed", color="black")
        lower = sum(aggregated["observed"][aggregated["pLI_bin"] <= 0.45])
        upper = sum(aggregated["observed"][aggregated["pLI_bin"] > 0.45])
        
        y_pos = ax.get_ylim()[1] * 0.9
        
        e = ax.text(0.35, y_pos, "n={0:.0f}".format(lower), horizontalalignment='center')
        e = ax.text(0.60, y_pos, "n={0:.0f}".format(upper), horizontalalignment='center')
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("normalised observed vs expected")
    
    fig.savefig(output, format="pdf")
    pyplot.close()

def excess_de_novos_from_pLI(de_novos, expected, constraints):
    """ identify the excess de novo mutations in more constrained genes
    
    Estimate excess of mutations in top half of lof and missense pLI bins. The
    pLI score for each gene is a good indicator of whether a gene can tolerate
    loss-of-function mutations. We can stratify genes by pLI scores, and
    determine the ratio of observed mutations to expected mutations for each pLI
    bin. Under a null model, we expect mutations to occur randomly in genes
    independent of pLI score. We can determine the excess as the difference in
    number of mutations between the upper half of genes versus the lower half of
    genes. When we break the variants down by consequence type, we expect the
    biggest exess to ocur in VEP predicted loss-of-function mutations. Variants
    which VEP classified as missense, but which have a loss-of-function
    mechanism will also be enriched in genes in the upper half of pLI scores.
    
    In this way we can estimate the number of VEP-classified missense variants
    which have a loss-of-function mechanism.
    
    Args:
        de_novos: pandas DataFrame of coding candidate de novo mutations
        expected: pandas DataFrame of numbers of expected mutations given our
            cohort size for different consequence types.
        constraints: pandas DataFrame of constraint scores, including columns
            for hgnc symbols ("gene") and loss-of-function constraints ("pLI").
    """
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    # include the synonymous observed counts in the combined counts table, to
    # show that variants with no consequence are dispersed randomly by pLI score.
    synonymous = count_synonymous_per_gene(de_novos)
    recode = dict(zip(synonymous["hgnc"], synonymous["observed"]))
    merged["synonymous_observed"] = merged["hgnc"].map(recode)
    merged["synonymous_observed"][merged["synonymous_observed"].isnull()] = 0
    
    # identify which pLI quantile each gene falls into
    merged = include_constraints(merged, constraints)
    merged["pLI_bin"] = get_constraint_bins(merged, bins=20)
    
    # get the ratio of observed to expected within the pLI bins for different
    # consequence types
    missense_by_hi = aggregate(merged, consequences=["missense"], normalise=False)
    lof_by_hi = aggregate(merged, consequences=["lof"], normalise=False)
    synonymous_by_hi = aggregate(merged, consequences=["synonymous"], normalise=False)
    
    # plot the differences in observed vs expected for the different pLI bins
    plot_by_hi_bin(missense_by_hi, "ratio", "missense_excess.pdf", expected=1, count_halves=True)
    plot_by_hi_bin(lof_by_hi, "ratio", "lof_excess.pdf", expected=1, count_halves=True)
    plot_by_hi_bin(synonymous_by_hi, "ratio", "synonymous_excess.pdf", expected=1, count_halves=True)

def plot_proportion_known_by_pLI(de_novos, expected,  constraints, known):
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    # identify known haploinsufficient dominant genes
    mono = classify_monoallelic_genes(known)
    
    # identify which pLI quantile each gene falls into
    merged = include_constraints(merged, constraints)
    merged["pLI_bin"] = get_constraint_bins(merged, bins=20)
    
    merged["haploinsufficient"] = merged["hgnc"].isin(mono["haploinsufficient"])
    merged["nonhaploinsufficient"] = merged["hgnc"].isin(mono["nonhaploinsufficient"])
    
    hi = merged.pivot_table(rows="pLI_bin", values="haploinsufficient", aggfunc=[sum, len])
    non_hi = merged.pivot_table(rows="pLI_bin", values="nonhaploinsufficient", aggfunc=[sum, len])
    
    hi["ratio"] = hi["sum"]/hi["len"]
    non_hi["ratio"] = non_hi["sum"]/non_hi["len"]
    
    hi["pLI_bin"] = hi.index
    non_hi["pLI_bin"] = non_hi.index
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    e = ax.plot(hi["pLI_bin"], hi["ratio"])
    e = ax.plot(non_hi["pLI_bin"], non_hi["ratio"])
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("proportion known pathogenic")
    
    fig.savefig("temp.pdf", format="pdf")
    
    hi_lof = merged[merged["haploinsufficient"]].pivot_table(rows="pLI_bin", values="lof_observed", aggfunc=sum)
    non_hi_mis = merged[merged["nonhaploinsufficient"] | merged["haploinsufficient"]].pivot_table(rows="pLI_bin", values="missense_observed", aggfunc=sum)
    
    lof = merged.pivot_table(rows="pLI_bin", values="lof_observed", aggfunc=sum)
    mis = merged.pivot_table(rows="pLI_bin", values="missense_observed", aggfunc=sum)
    
    known_proportions = pandas.DataFrame({"pLI_bin": lof.index, "lof_ratio": hi_lof/lof, "mis_ratio": non_hi_mis/mis})
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    e = ax.plot(known_proportions["pLI_bin"], known_proportions["lof_ratio"])
    e = ax.plot(known_proportions["pLI_bin"], known_proportions["mis_ratio"])
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("proportion known pathogenic")
    
    fig.savefig("temp2.pdf", format="pdf")

def main():
    de_novos = open_de_novos(DENOVO_PATH, VALIDATIONS, exclude_synonymous=False)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    de_novos["start_pos"] = de_novos["pos"]
    
    male = 2407
    female = 1887
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    constraints = pandas.read_table(CONSTRAINTS_URL)
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    # filtered = de_novos[(de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    filtered = de_novos[(~de_novos["pp_dnm"].isnull() & (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    
    ratios = get_overall_consequence_ratios(expected, filtered)
    plot_overall_ratios(ratios, "results/excess_by_consequence.pdf")
    
    proportions = model_mixing(known, filtered, expected, constraints)
    
    print(proportions)
    
    excess_de_novos_from_pLI(filtered, expected, constraints)
    
    plot_proportion_known_by_pLI(filtered, expected,  constraints, known)

if __name__ == '__main__':
    main()
