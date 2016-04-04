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

import math

from scipy.stats import poisson, norm

from ddd_4k.causation.classify_known_genes import classify_monoallelic_genes

def prevalence_from_cohort_excess(cohort_n, excess, cnv_yield=0.1, missing_variants=119,
    enrichment=120.0):
    """ calculate the prevalence of live births with developmental disorders
    caused by de novo mutations.
    
    This prevalence estimate is based off the proportion of the cohort with an
    excess de novo mutation, which is scaled to the general population by an
    enrichment factor, determined from the excess enrichment of loss-of-function
    mutations within dominant haploinsufficient neurodevelopmental genes with
    low clinical recognisability.
    
    Args:
        cohort_n: size of the studied cohort
        excess: dictionary with numbers of excess mutations for each of
            synonymous, missense and loss-of-function mutations
        cnv_yield: estimated diagnostic yield from de novo CNVs in probands with
            developmental disorders, obtained from published estimates
            (doi:10.1038/ng.909, doi:10.1097/GIM.0b013e318194ee8f and
            doi:10.1016/j.ajhg.2010.04.006)
        missing_variants: number of variants missing from the cohort due to ease
             of clinical recognisability.
        enrichment: factor by which the cohort is enriched for loss-of-function
            de novo mutations in known haploinsufficient genes with low clinical
            recognisability.
    
    Returns:
        estimate of the prevalence of live births with developmental disorders
        caused by de novo mutations
    """
    
    cohort_yield = (excess["loss-of-function"]["excess"] + excess["missense"]["excess"])/cohort_n
    missing_yield = missing_variants/cohort_n + cnv_yield
    
    # adjust the yield in the DDD for the proportions that exist outside the
    # DDD, since these expand the DDD cohort size.
    cohort_yield = cohort_yield/(1 + missing_yield)
    
    return (cohort_yield + missing_yield)/enrichment

def prevalence_from_baseline_lof(rates, known, excess_to_lof, snv_yield, unobserved_snvs, unobserved_cnvs, ci=None):
    """ estimate the prevalence of developmental disorders caused by de novo
    mutations from known mutation rates for the genes.
    
    Args:
        rates: pandas DataFrame of expected mutation rates for all genes in the
            genome (or near all).
        known: pandas DataFrame of genes known to be involved in developmental
            disorders, including the mode of inheritance and mechanism of action.
        excess_to_lof: dictionary of counts of excess de novos and
            loss-of-function de novos in dominant DD-associated genes.
        snv_yield: yield of excess SNVs within the cohort_yield
        unobserved_snvs: number of SNVs missing from the cohort due to ease of
             clinical recognisability.
        unobserved_cnvs: dictionary of counts of excess
        ci: confidence interval (for example, 0.95 for 95% CI) or None if CI not
            required.
    
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
    lof_prevalence = poisson.sf(0, lof_mu, loc=0)
    
    # scale up by the ratio of excess de novos to loss-of-function de novos in
    # DD-associated genes.
    ratio = (excess_to_lof['excess'] + unobserved_snvs)/(excess_to_lof['dominant_lof'] + unobserved_snvs)
    prevalence_baseline = lof_prevalence * ratio
    
    cnv_yield = unobserved_cnvs['diagnosed']/float(unobserved_cnvs['cohort_n'])
    cnv_boost = ((snv_yield + cnv_yield)/(1 + cnv_yield))/snv_yield
    prevalence_adjusted = prevalence_baseline * cnv_boost
    
    if ci is not None:
        ci_interval = prevalence_CI(prevalence_adjusted, excess_to_lof, snv_yield,
            unobserved_snvs, unobserved_cnvs, ci)
        return (prevalence_adjusted, ci_interval)
    else:
        return prevalence_adjusted

def prevalence_CI(prevalence, excess_to_lof, snv_yield, unobserved_snvs, unobserved_cnvs, ci=0.95):
    ''' estimate the confidence interval for the prevalance of developmental
    disorders caused by de novo mutations.
    
    Args:
        prevalence: estimate of prevalance of DD due to dominant de novo mutations
        snv_yield: yield of excess SNVs within the cohort_yield
        unobserved_snvs: number of SNVs missing from the cohort due to ease of
             clinical recognisability.
        unobserved_cnvs: dictionary of counts of excess
        ci: confidence interval (for example, 0.95 for 95% CI)
    
    Returns:
        tuple of estimates for lower and upper confidence interval for the
        prevalance of live births with developmental disorders caused by de novo
        mutations.
    '''
    
    z = norm.ppf(1 - ((1 - ci)/2))
    
    # get the normal approximation to the CI for the ratio of excess to dominant
    # LoF in dominant DD-associated genes
    ratio = (excess_to_lof['excess'] + unobserved_snvs)/(excess_to_lof['dominant_lof'] + unobserved_snvs)
    inverse = 1/ratio
    delta = z * math.sqrt((inverse * (1 - inverse))/(excess_to_lof['excess'] + unobserved_snvs))
    
    ratio_lower = 1/(inverse + delta)
    ratio_upper = 1/(inverse - delta)
    
    # get the normal approximation to the CI for the boost in prevalence due to
    # unobserved de novo CNVs.
    cnv_yield = unobserved_cnvs['diagnosed']/float(unobserved_cnvs['cohort_n'])
    delta = z * math.sqrt((cnv_yield * (1 - cnv_yield))/unobserved_cnvs['cohort_n'])
    
    cnv_yield_lower = cnv_yield - delta
    cnv_yield_upper = cnv_yield + delta
    
    cnv_boost = ((snv_yield + cnv_yield)/(1 + cnv_yield))/snv_yield
    cnv_boost_lower = ((snv_yield + cnv_yield_lower)/(1 + cnv_yield_lower))/snv_yield
    cnv_boost_upper = ((snv_yield + cnv_yield_upper)/(1 + cnv_yield_upper))/snv_yield
    
    # scale the prevalence up and down to get the upper an lower confidence inter
    prevalence_lower = prevalence * (ratio_lower/ratio) * (cnv_boost_lower/cnv_boost)
    prevalence_upper = prevalence * (ratio_upper/ratio) * (cnv_boost_upper/cnv_boost)
    
    return (prevalence_lower, prevalence_upper)
