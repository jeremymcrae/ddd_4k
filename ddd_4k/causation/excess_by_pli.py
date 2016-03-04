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

from matplotlib import pyplot, gridspec

from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.count_synonymous import count_synonymous_per_gene
from ddd_4k.causation.pli_functions import include_constraints, get_constraint_bins
from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.plot_pli_bins import plot_by_hi_bin

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
    merged["pLI_bin"] = get_constraint_bins(merged, bins=20, rate_correct=True)
    
    # get the ratio of observed to expected within the pLI bins for different
    # consequence types
    lof_by_hi = aggregate(merged, consequences=["lof"], normalise=False)
    missense_by_hi = aggregate(merged, consequences=["missense"], normalise=False)
    synonymous_by_hi = aggregate(merged, consequences=["synonymous"], normalise=False)
    
    fig = pyplot.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(3, 1 )
    
    lof_ax = pyplot.subplot(gs[0])
    mis_ax = pyplot.subplot(gs[1])
    syn_ax = pyplot.subplot(gs[2])
    
    # plot the differences in observed vs expected for the different pLI bins
    plot_by_hi_bin(lof_by_hi, "ratio", title="loss-of-function", expected=1, count_halves=True, ax=lof_ax)
    plot_by_hi_bin(missense_by_hi, "ratio", title="missense", expected=1, count_halves=True, ax=mis_ax)
    plot_by_hi_bin(synonymous_by_hi, "ratio", title="synonymous", expected=1, count_halves=True, ax=syn_ax)
    
    fig.savefig("results/excess_consequence_by_pli.pdf", format="pdf",
        bbox_inches='tight', pad_inches=0, transparent=True)
    pyplot.close()
