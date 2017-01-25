'''
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
'''

from __future__ import absolute_import, print_function, division

import argparse
import math

from scipy.stats import norm
import pandas

from mupit.open_ddd_data import open_known_genes
from ddd_4k.load_files import open_de_novos, open_phenotypes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, \
    CONSTRAINTS_URL, PHENOTYPES, TRIOS, SANGER_IDS, PREVIOUS_VALIDATIONS
from ddd_4k.causation.excess_by_consequence import get_consequence_excess, \
    plot_consequence_excess
from ddd_4k.causation.model_mixtures import model_mixing
from ddd_4k.causation.prevalence import plot_prevalence_by_age
from ddd_4k.causation.de_novo_threshold import get_pp_dnm_threshold
from ddd_4k.validation_rates import open_previous_validations, get_rates, plot_excess_by_pp_dnm_threshold
from ddd_4k.causation.excess_by_pli import excess_de_novos_from_pLI
from ddd_4k.causation.proportion_known_by_pli import plot_proportion_known_by_pLI
from ddd_4k.causation.open_uk_ages import open_uk_parent_ages
from ddd_4k.causation.prevalence_estimates import prevalence_from_baseline_lof

from mupit.mutation_rates import get_default_rates, get_expected_mutations

UK_AGES = '/nfs/users/nfs_j/jm33/apps/ddd_4k/data/Father_Mother_AgesAtBirth.xlsx'

def get_options():
    ''' parse the command line arguments
    '''
    
    parser = argparse.ArgumentParser(description='script to probe the limits '
        'of de novo causation in children with developmental disorders.')
    parser.add_argument('--de-novos', default=DENOVO_PATH,
        help='Path to table of candidate de novo mutations.')
    parser.add_argument('--validations', default=VALIDATIONS,
        help='Path to table of validation data.')
    parser.add_argument('--trios', default=TRIOS,
        help='Path to table of complete trios in the cohort.')
    parser.add_argument('--sanger-ids', default=SANGER_IDS,
        help='Path to table mapping between different ID formats.')
    parser.add_argument('--phenotypes', default=PHENOTYPES,
        help='Path to table of phenotypic data for probands.')
    parser.add_argument('--known-genes', default=KNOWN_GENES,
        help='Path to table of known developmental disorder genes.')
    parser.add_argument('--previous-validations', default=PREVIOUS_VALIDATIONS,
        help='Path to table of probands with known diagnostic variants.')
    parser.add_argument('--constraints', default=CONSTRAINTS_URL,
        help='Path or URL to table of constraint scores (pLI-based) for all'
            'genes in the genome.')
    parser.add_argument('--uk-ages', default=UK_AGES,
        help='Path to table of ages of parents at child\'s birth for mothers '
            'and fathers in the UK.')
    
    args = parser.parse_args()
    
    return args

def count_known_excess(filtered, known):
    ''' count the filtered de novos in DD-associated dominant genes.
    
    Protein altering de novos in DD-associated dominant genes are likely
    pathogenicand so must be part of the excess of de novo burden within the
    cohort.
    
    Args:
        filtered: pandas DataFrame of de novo variants.
        known: pandas DataFrame of DD-associated genes.
    
    Returns:
        dictionary of counts for truncating and 'missense' de novos.
    '''
    
    dominant = known['gene'][known['mode'].isin(['Monoallelic', 'X-linked dominant'])]
    variants = filtered[filtered['symbol'].isin(dominant)].copy()
    
    counts = variants['category'].value_counts()
    counts['missense'] = counts['functional']
    
    return dict(counts)

def print_known_in_excess(in_dominant, excess):
    lof_in_dominant = in_dominant['truncating']
    mis_in_dominant = in_dominant['missense']
    lof_excess = excess['truncating']['excess']
    mis_excess = excess['missense']['excess']
    
    print('excess lof in known: {:.0f}%'.format(100 * lof_in_dominant/lof_excess))
    print('excess mis in known: {:.0f}%'.format(100 * mis_in_dominant/mis_excess))
    
    print('excess func in known: {:.0f}%'.format(
        100 * (lof_in_dominant + mis_in_dominant)/(lof_excess + mis_excess)))

def proportion_in_dominant_hi_genes(filtered, known):
    ''' count the consequence types in dominant haploinsufficient DD genes
    
    Within the set of dominant haploinsufficient DD-associated genes, examine
    the breakdown of loss-of-function and missense consequences.
    
    Returns:
        dictionary of counts for truncating and 'missense' de novos
    '''
    
    dominant = known[known['mode'].isin(['Monoallelic', 'X-linked dominant'])]
    dominant_haploinsufficient = dominant['gene'][dominant['mech'] == 'Loss of function']
    dominant_haploinsufficient = filtered['symbol'].isin(dominant_haploinsufficient)
    
    # count the variants in genes with a dominant but not HI mode
    dom_nonhi_genes = dominant['gene'][dominant['mech'].isin(['Activating', 'Dominant negative'])]
    dom_nonhi = filtered['symbol'].isin(dom_nonhi_genes)
    nonhi_counts = filtered['category'][dom_nonhi].value_counts()
    print(nonhi_counts)
    
    variants = filtered[dominant_haploinsufficient].copy()
    
    counts = variants['category'].value_counts()
    counts['missense'] = counts['functional']
    
    return dict(counts)

def print_proportions_from_dominant_hi(dominant_hi, excess, alpha=0.95):
    ''' show the proportions of loss-of-function/altered function
    
    We can assume the protein truncating variants have a loss-of-function
    mechanism. In order to figure out how many of the missense excess has a
    loss-of-function mechanism, we use the numbers of  truncating and missense
    de novos in known DD-associated genes.
    
    Args:
        dominant_hi: dictionary of counts of de novos in dominant DD-associated
            genes. The entries are for 'missense' and 'truncating'
            consequences e.g. {'missense': 50, 'truncating': 60}.
        excess: dictionary of counts of excess de novos within the cohort,
            broken down by whether the excess is for missense, truncating or
            synonymous e.g. {'missense': {'excess': 100}, 'truncating': {'excess': 200}}.
        alpha: alpha for confidence intervals.
    '''
    
    lof_as_mis = dominant_hi['missense']/dominant_hi['truncating']
    
    lof_excess = excess['truncating']['excess']/ \
        (excess['missense']['excess'] + excess['truncating']['excess'])
    
    lof_proportion = lof_excess + lof_excess * lof_as_mis
    gof_proportion = 1 - lof_proportion
    
    # also estimate a 95% confidence interval
    z = norm.ppf(1 - ((1 - alpha)/2))
    hi_delta = z * math.sqrt((lof_as_mis * (1 - lof_as_mis))/(dominant_hi['missense'] + dominant_hi['truncating']))
    excess_delta = z * math.sqrt((lof_excess * (1 - lof_excess))/(excess['missense']['excess'] + excess['truncating']['excess']))
    lower = (lof_excess - excess_delta) + (lof_excess - excess_delta) * (lof_as_mis - hi_delta)
    upper = (lof_excess + excess_delta) + (lof_excess + excess_delta) * (lof_as_mis + hi_delta)
    
    print('LoF: {:.0f}%'.format(lof_proportion * 100))
    print('GoF: {:.0f}%'.format(gof_proportion * 100))
    print('LoF 95% CI: {0}-{1}%'.format(lower * 100, upper * 100))

def main():
    
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations, exclude_synonymous=False)
    known = open_known_genes(args.known_genes)
    de_novos['known'] = de_novos['symbol'].isin(known['gene'])
    de_novos['start_pos'] = de_novos['pos']
    constraints = pandas.read_table(args.constraints)
    uk_ages = open_uk_parent_ages(args.uk_ages)
    
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    trios = pandas.read_table(args.trios)
    phenotypes = phenotypes[phenotypes['patient_id'].isin(trios['decipher_id'])]
    
    male = 2408
    female = 1885
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    filtered = de_novos[(~de_novos['pp_dnm'].isnull() & (de_novos['pp_dnm'] > pp_dnm_threshold)) ]
    
    complete = open_de_novos(args.de_novos, exclude_synonymous=False)
    validations = open_previous_validations(args.previous_validations, complete)
    validation_rates = get_rates(validations, pp_dnm_threshold)
    
    excess = get_consequence_excess(expected, filtered, validation_rates.ppv, validation_rates.tpr)
    functional_excess = excess['truncating']['excess'] + excess['missense']['excess']
    snv_yield = functional_excess / (male + female)
    
    plot_consequence_excess(excess, 'results/excess_by_consequence.pdf')
    plot_excess_by_pp_dnm_threshold(de_novos, expected, validations, increments=100)
    
    dominant_hi = proportion_in_dominant_hi_genes(filtered, known)
    print_proportions_from_dominant_hi(dominant_hi, excess)
    
    in_dominant = count_known_excess(filtered, known)
    print_known_in_excess(in_dominant, excess)
    
    proportions = model_mixing(known, filtered, expected, constraints,
        check_modelling=False, check_variance=False)
    print('missense proportion as LoF: {}'.format(proportions))
    
    print((excess['truncating']['excess'] + excess['missense']['excess'] * proportions) / functional_excess)
    
    excess_de_novos_from_pLI(filtered, expected, constraints)
    
    # define the number of individuals with developmental disorders caused by
    # a de novo CNV. We do not observe these individuals in the DDD cohort.
    # These totals are for the whole genome-only cohorts listed in Table 2 of
    # Miller et al, ASHG 86:749-764, doi:10.1016/j.ajhg.2010.04.006
    unobserved_cnvs = {'cohort_n': 2159, 'diagnosed': 204}
    
    # define the number of SNVs that would be in genes with high clinical
    # recognisability, but which have not made it to our cohort. This is
    # estimated within scripts/clinical_recognisability.py
    unobserved_snvs = 119
    
    excess_to_lof = {'excess': functional_excess, 'dominant_lof': in_dominant['truncating']}
    prevalance_from_rates, ci = prevalence_from_baseline_lof(rates, known,
        excess_to_lof, snv_yield=snv_yield, unobserved_snvs=unobserved_snvs,
        unobserved_cnvs=unobserved_cnvs, ci=0.95)
    
    # define the number of de novo mutations per child from parents. Numbers
    # provided by Raheleh Rahbari, based on Rahbari et al, Nature Genetics 2016
    # 48: 126-133 doi:10.1038/ng.3469
    reference_mutations = {'dad_age': 29.53, 'mom_age': 29.87, 'mutations': 77}
    
    plot_prevalence_by_age(prevalance_from_rates, phenotypes, uk_ages,
        reference_mutations, dad_rate=1.53, mom_rate=0.86)
    print('prevalence: {}'.format(prevalance_from_rates))
    print('prevalence 95% CI: {}'.format(ci))
    
    print('SNV yield: {}'.format(snv_yield))
    

if __name__ == '__main__':
    main()
