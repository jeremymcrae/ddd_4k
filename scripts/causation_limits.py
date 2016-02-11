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

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

HAPLOINSUFFICIENCY_URL = "http://files.figshare.com/410746/Dataset_S1.txt"

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
        de_novos: pandas dataframe of all candidate coding de novos, including
            synonymous variants
        rates: pandas dataframe of numbers of expected mutations within the
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
        synonymous = de_novos[(de_novos["consequence"] == "synonymous_variant") &
            (de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > mid)) ]
        synonymous = pandas.pivot_table(synonymous, rows="symbol",
            values="alt", aggfunc=len)
        synonymous = pandas.DataFrame({"hgnc": synonymous.index,
            "synonymous_snv": synonymous})
        
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

def get_overall_consequence_ratios(expected, de_novos):
    """ determine the ratio of obeserved to expected for synonymous, missense
    and loss-of-function canddiate de novos
    
    Args:
        expected: pandas dataframe of counts of expected de novo mutations for
            all genes (or nearly all) in the genome.
        de_novos: pandas dataframe of candidate de novo mutations observed in
            the cohort.
    
    Returns:
        dictionary of counts and ratios for synonymous, missense and
        loss-of-function consequences.
    """
    
    observed = get_de_novo_counts(de_novos)
    synonymous_count = len(de_novos[de_novos["consequence"] == "synonymous_variant"])
    
    # sum the observed loss-of-function and missense candidates per gene
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
    
    # get the numbers of expecetd loss-of-function and missense mutations.
    ratios["lof_expected"] = ratios[["lof_snv", "lof_indel"]].sum(axis=1)
    ratios["missense_expected"] = ratios[["missense_snv", "missense_indel"]].sum(axis=1)
    
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
        excess = row["count"] - row["count"]/row["ratio"]
        e = ax.text(key, row["ratio"]+0.01, "n={0:.0f}\nexcess={1:.0f}".format(row["count"], excess), horizontalalignment='center')
    
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

def main():
    de_novos = open_de_novos(DENOVO_PATH, VALIDATIONS, exclude_synonymous=False)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    de_novos["start_pos"] = de_novos["pos"]
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    de_novos = de_novos[(de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    
    male = 2407
    female = 1887
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    ratios = get_overall_consequence_ratios(expected, de_novos)
    plot_overall_ratios(ratios, "results/excess_by_consequence.pdf")
    
    constraints = pandas.read_table(CONSTRAINTS_URL)
    recode = dict(zip(constraints["gene"], constraints["pLI"]))
    
    ratios["pLI"] = ratios["hgnc"].map(recode)
    
    quantiles = [ x/20 for x in range(21) ]
    ratios["bin"], bins = pandas.qcut(ratios["pLI"], q=quantiles, labels=quantiles[:-1], retbins=True)
    
    lof_enrich = pandas.pivot_table(ratios, values="lof_ratio", rows=["bin"], aggfunc=numpy.mean, fill_value=0)
    mis_enrich = pandas.pivot_table(ratios, values="mis_ratio", rows=["bin"], aggfunc=numpy.mean, fill_value=0)
    
    # haploinsufficiency = open_haploinsufficiency(HAPLOINSUFFICIENCY_URL)
    # enrichment = get_enrichment_ratios(constraints, haploinsufficiency)
    
    # # plot the enrichment ratio for each bin
    # fig = seaborn.lmplot(x="pLI", y="enrichment", data=enrichment, size=6, lowess=True)
    # fig.savefig("results/enrichment_ratio.pdf", type="pdf")

if __name__ == '__main__':
    main()
