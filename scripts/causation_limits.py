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

from ddd_4k.load_files import open_de_novos, open_known_genes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, CONSTRAINTS_URL
from ddd_4k.count_mutations_per_person import get_count_by_person

from mupit.mutation_rates import get_default_rates, get_expected_mutations
from mupit.count_de_novos import get_de_novo_counts
from mupit.gene_enrichment import gene_enrichment

def open_haploinsufficiency(url):
    """ load the haploinsufficiency dataset
    
    Args:
        url: URL for haploinsufficiency dataset.
    
    Returns:
        pandas DataFrame of haploinsufficiency scores by gene symbol.
    """
    
    # get haploinsufficiency predictions genomewide. The data is from a report
    # in PLOS Genetics, formatted as a bed file:
    # PLOS Genet 6:e1001154 - doi:10.1371/journal.pgen.1001154
    # TODO: look for updated list of haploinsufficiency, with newer gene symbols
    # and get better threshold to classify haploinsufficiency.
    hi = pandas.read_table(url, skiprows=1, header=None)
    hi.columns = ["chrom", "start", "end", "name", "score", "strand",
        "thick_start", "thick_end", "color"]
    
    # reformat the gene and chromosome from the table
    hi["gene"] = [ x.split("|")[0] for x in hi["name"] ]
    hi["chrom"] = [ x.strip("chr") for x in hi["chrom"] ]
    
    return hi

def get_enrichment_ratios(constraints, haploinsufficiency):
    """ figure out the haploinsufficiency enrichment by constraint bin
    
    Args:
        constraints: pandas DataFrame of constraint scores by gene
        haploinsufficiency: pandas DataFrame of haploinsufficiency scores by gene
    
    Returns:
        pandas DataFrame of enrichment ratios for each pLI bin.
    """
    
    # classify each gene as belonging to one of 20 evenly spaced bins
    constraints["bin"], bins = pandas.cut(constraints["pLI"], bins=20, retbins=True)
    
    merged = constraints.merge(haploinsufficiency, on="gene")
    groups = merged.groupby("bin")
    
    # determine the enrichment ratio in each pLI bin
    threshold = 0.9
    values = [ sum(x["score"] > threshold)/len(x.index) for i, x in groups ]
    ratios = [ x/values[0] for x in values ]
    enrichment = pandas.DataFrame({"pLI": bins[:-1], "enrichment": ratios})
    
    return enrichment

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
        
        # select the synonymous candidates above the midpoint threshold
        syn = de_novos[(de_novos["consequence"] == "synonymous_variant") &
            (de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > mid)) ]
        syn = syn.pivot_table(rows="symbol", values="alt", aggfunc=len)
        synonymous = pandas.DataFrame({"hgnc": syn.index, "synonymous_snv": syn})
        
        # make sure we have expected numbers available. Some genes lack
        # expected numbers, so we exclude those genes, rather than counting
        # extra genes, for which we can never know if we are near expectations
        synonymous["expected"] = synonymous["hgnc"].map(rates)
        synonymous = synonymous[~synonymous["expected"].isnull()]
        
        ratio = sum(synonymous["synonymous_snv"])/sum(expected["synonymous_snv"])
        
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
    ratios = expected[["hgnc", "lof_indel", "lof_snv", "missense_indel", "missense_snv"]].copy()
    ratios = ratios.merge(observed, how="left", on="hgnc")
    ratios["lof_observed"][ratios["lof_observed"].isnull()] = 0
    ratios["missense_observed"][ratios["missense_observed"].isnull()] = 0
    
    # get the numbers of expected loss-of-function and missense mutations.
    ratios["lof_expected"] = ratios[["lof_snv", "lof_indel"]].sum(axis=1)
    ratios["missense_expected"] = ratios[["missense_snv", "missense_indel"]].sum(axis=1)
    
    # restrict the table to the minimum required columns
    ratios = ratios[["hgnc", "lof_observed", "missense_observed",
        "lof_expected", "missense_expected"]].copy()
    
    return ratios

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
    
    ratios = merge_observed_and_expected(de_novos, expected)
    
    lof_ratio = sum(ratios["lof_observed"])/sum(ratios["lof_expected"])
    missense_ratio = sum(ratios["missense_observed"])/sum(ratios["missense_expected"])
    
    return {"synonymous": {"ratio": 1.0, "count": synonymous_count},
        "loss-of-function": {"ratio": lof_ratio, "count": sum(ratios["lof_observed"])},
        "missense": {"ratio": missense_ratio, "count": sum(ratios["missense_observed"])}}

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

def model_mixing(known, de_novos, expected):
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
    
    # identify known haploinsufficient dominant genes
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
    
    ratios = merge_observed_and_expected(de_novos, expected)
    
    # load the pLI data, and append column to the table
    constraints = pandas.read_table(CONSTRAINTS_URL)
    recode = dict(zip(constraints["gene"], constraints["pLI"]))
    ratios["pLI"] = ratios["hgnc"].map(recode)
    ratios = ratios[~ratios["pLI"].isnull()]
    
    # identify which pLI quantile each gene falls into
    quantiles = [ x/20.0 for x in range(21) ]
    ratios["pLI_bin"], bins = pandas.qcut(ratios["pLI"], q=quantiles, labels=quantiles[:-1], retbins=True)
    
    hi_ratios = ratios[ratios["hgnc"].isin(haploinsufficient)]
    non_hi_ratios = ratios[ratios["hgnc"].isin(nonhaploinsufficient)]
    
    plot_by_hi_bin(ratios, ["missense", "lof"], "results/obs_to_exp_delta_by_hi_bin.all_genes.pdf")
    plot_by_hi_bin(hi_ratios, ["lof"], "results/obs_to_exp_delta_by_hi_bin.hi_genes.pdf")
    plot_by_hi_bin(non_hi_ratios, ["missense"], "results/obs_to_exp_delta_by_hi_bin.non_hi_genes.pdf")
    
    target = aggregate(ratios, ["missense", "lof"])
    hi_start = aggregate(hi_ratios, ["lof"])
    non_hi_start = aggregate(non_hi_ratios, ["missense"])
    
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

def aggregate(ratios, consequences=None, normalise=True):
    """ get the difference between observed and expected counts across pLI bins.
    
    This returns the difference in the number of observed de novo mutations to
    expected numbers of de novo mutations. The difference is given as a delta
    and as a ratio.
    
    Args:
        ratios: pandas DataFrame containing one row per gene, which includes
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
    
    if consequences is not None:
        consequences = ["lof", "missense"]
        
    observed_columns = [ x + "_observed" for x in consequences ]
    expected_columns = [ x + "_expected" for x in consequences ]
    
    aggregated = ratios.pivot_table(rows="pLI_bin",
        values=observed_columns + expected_columns, aggfunc=sum)
    
    aggregated["pLI_bin"] = aggregated.index
    
    observed = aggregated[observed_columns].sum(axis=1)
    expected = aggregated[expected_columns].sum(axis=1)
    
    aggregated["delta"] = observed - expected
    aggregated["ratio"] = observed/expected
    
    if normalise:
        aggregated["delta"] = aggregated["delta"]/sum(aggregated["delta"])
        aggregated["ratio"] = aggregated["ratio"]/sum(aggregated["ratio"])
    
    # make sure we have all the quantile bins, even for tables that might lack
    # genes in a given pLI bin
    for x in [ x/20.0 for x in range(20) ]:
        if x not in list(aggregated["pLI_bin"]):
            row = pandas.DataFrame({"pLI_bin": [x], "delta": [0], "ratio": [0]})
            aggregated = aggregated.append(row, ignore_index=True)
    
    aggregated = aggregated[["pLI_bin", "delta", "ratio"]].copy()
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

def plot_by_hi_bin(ratios, consequences, output):
    """ plot difference in observed to expected across pLI vigiciles
    
    Args:
        ratios: pandas DataFrame containing observed and expected counts per
            gene, along with pLI categories for each gene.
        consequences: list of consequence types to include, for example
            ["missense"] or ["missense", "lof"] or ["lof"]
        output: path to save the plot to.
    """
    
    aggregated = aggregate(ratios, consequences)
    
    fig = pyplot.figure(figsize=(12, 6))
    ax = fig.gca()
    
    barwidth = 0.04
    e = ax.bar(aggregated["pLI_bin"], aggregated["delta"], width=barwidth,
        align="center")
    
    # fix the axis limits and ticks
    quantiles = [ x/20.0 for x in range(20) ]
    e = ax.set_xlim((-0.05, max(quantiles) + 0.05))
    e = ax.set_xticks(quantiles)
    e = ax.set_xticklabels(quantiles, fontsize="medium")
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("normalised observed - expected")
    
    fig.savefig(output, format="pdf")
    pyplot.close()

def main():
    de_novos = open_de_novos(DENOVO_PATH, VALIDATIONS, exclude_synonymous=False)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    de_novos["start_pos"] = de_novos["pos"]
    
    male = 2407
    female = 1887
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    de_novos = de_novos[(de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    
    ratios = get_overall_consequence_ratios(expected, de_novos)
    plot_overall_ratios(ratios, "results/excess_by_consequence.pdf")
    
    proportions = model_mixing(known, de_novos, expected)
    
    print(proportions)

if __name__ == '__main__':
    main()
