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

from __future__ import print_function, division

from scipy.stats import poisson

from ddd_4k.causation.classify_known_genes import classify_monoallelic_genes

def check_prevalence_from_baseline_lof(rates, known, mis_to_lof=2.0,
        missing=0.5, cnv_adjust=1.2105):
    """ estimate the prevalence of developmental disorders caused by de novo
    mutations from known mutation rates for the genes.
    
    Args:
        rates: pandas DataFrame of expected mutation rates for all genes in the
            genome (or near all).
        known: pandas DataFrame of genes known to be involved in developmental
            disorders, including the mode of inheritance and mechanism of action.
        mis_to_lof: ratio of  the number of missense mutations to
            loss-of-function mutations.
        missing: proportion of genes yet to be discovered (relative to the
            currently known genes).
        cnv_adjust: the ratio of predicted diagnostic yield of DDD probands with
            any dominant de novo to the predicted yield without CNV de novos.
    
    Returns:
        estimate for the prevalance of live births with developmental disorders
        caused by de novo mutations.
    """
    
    mono = classify_monoallelic_genes(known, remove_overlap=False)
    
    # restrict the rates to the known haploinsufficient developmental disorder
    # genes. Six of 238 genes lack mutation rates, which shouldn't substantially
    # affect downstream estimates.
    dominant = rates[rates["hgnc"].isin(mono["haploinsufficient"])]
    
    # get the summed loss-of-function mutation rate in known dominant genes
    lof_mu = sum(dominant[["splice_site", "frameshift", "non"]].sum(axis=1, skipna=True))
    
    # get the proportion of the population with a loss-of-function mutation in a
    # known dominant gene.
    lof_proportion = poisson.sf(0, lof_mu, loc=0)
    
    # get the corresponding proportion of the population with a missense
    # mutation in a known dominant gene.
    mis_proportion = lof_proportion * mis_to_lof
    
    # get the baseline proportion of the population with a mutation in a known
    # dominant gene
    prevalence_baseline = lof_proportion + mis_proportion
    
    # adjust the prevalance for the proportion of genes yet to be discovered.
    prevalence_adjusted = prevalence_baseline * (1.0 + missing) * cnv_adjust
    
    return prevalence_adjusted