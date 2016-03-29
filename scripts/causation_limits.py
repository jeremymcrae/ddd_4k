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

from __future__ import absolute_import, print_function, division

import argparse

import pandas

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, \
    CONSTRAINTS_URL, PHENOTYPES, TRIOS, SANGER_IDS, DIAGNOSED
from ddd_4k.causation.excess_by_consequence import get_consequence_excess, \
    plot_consequence_excess
from ddd_4k.causation.excess_by_pp_dnm_threshold import plot_excess_by_pp_dnm_threshold
from ddd_4k.causation.model_mixtures import model_mixing
from ddd_4k.causation.prevalence import get_birth_prevalence, plot_prevalence_by_age
from ddd_4k.causation.de_novo_threshold import get_pp_dnm_threshold
from ddd_4k.causation.excess_by_pli import excess_de_novos_from_pLI
from ddd_4k.causation.proportion_known_by_pli import plot_proportion_known_by_pLI
from ddd_4k.causation.open_uk_ages import open_uk_parent_ages
from ddd_4k.causation.prevalence_from_baseline import check_prevalence_from_baseline_lof

from mupit.mutation_rates import get_default_rates, get_expected_mutations

UK_AGES = "/nfs/users/nfs_j/jm33/apps/ddd_4k/data/Father_Mother_AgesAtBirth.xlsx"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to probe the limits " \
        "of de novo causation in children with developmental disorders.")
    parser.add_argument("--de-novos", default=DENOVO_PATH,
        help="Path to table of candidate de novo mutations.")
    parser.add_argument("--validations", default=VALIDATIONS,
        help="Path to table of validation data.")
    parser.add_argument("--trios", default=TRIOS,
        help="Path to table of complete trios in the cohort.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS,
        help="Path to table mapping between different ID formats.")
    parser.add_argument("--phenotypes", default=PHENOTYPES,
        help="Path to table of phenotypic data for probands.")
    parser.add_argument("--known-genes", default=KNOWN_GENES,
        help="Path to table of known developmental disorder genes.")
    parser.add_argument("--diagnosed", default=DIAGNOSED,
        help="Path to table of probands with known diagnostic variants.")
    parser.add_argument("--constraints", default=CONSTRAINTS_URL,
        help="Path or URL to table of constraint scores (pLI-based) for all" \
            "genes in the genome.")
    parser.add_argument("--uk-ages", default=UK_AGES,
        help="Path to table of ages of parents at child's birth for mothers " \
            "and fathers in the UK.")
    
    args = parser.parse_args()
    
    return args

def count_known_excess(filtered, known):
    """ count the filtered de novos in DD-associated dominant genes.
    
    Protein altering de novos in DD-associated dominant genes are likely
    pathogenicand so must be part of the excess of de novo burden within the
    cohort.
    
    Args:
        filtered: pandas DataFrame of de novo variants.
        known: pandas DataFrame of DD-associated genes.
    
    Returns:
        dictionary of counts for loss-of-function and 'missense' de novos.
    """
    
    dominant = known['gencode_gene_name'][known['mode'].isin(['Monoallelic', 'X-linked dominant'])]
    variants = filtered[filtered["symbol"].isin(dominant)].copy()
    
    counts = variants['category'].value_counts()
    counts['missense'] = counts['functional']
    
    return dict(counts)

def print_known_in_excess(in_dominant, excess):
    lof_in_dominant = in_dominant['loss-of-function']
    mis_in_dominant = in_dominant['missense']
    lof_excess = excess['loss-of-function']['excess']
    mis_excess = excess['missense']['excess']
    
    print('excess lof in known: {:.0f}%'.format(100 * lof_in_dominant/lof_excess))
    print('excess mis in known: {:.0f}%'.format(100 * mis_in_dominant/mis_excess))
    
    print('excess func in known: {:.0f}%'.format(
        100 * (lof_in_dominant + mis_in_dominant)/(lof_excess + mis_excess)))

def proportion_in_dominant_hi_genes(filtered, known):
    """ count the consequence types in dominant haploinsufficient DD genes
    
    Within the set of dominant haploinsufficient DD-associated genes, examine
    the breakdown of loss-of-function and missense consequences.
    
    Returns:
        dictionary of counts for loss-of-function and 'missense' de novos
    """
    
    dominant = known[known['mode'].isin(['Monoallelic', 'X-linked dominant'])]
    dominant_haploinsufficient = dominant['gencode_gene_name'][dominant['mech'] == 'Loss of function']
    dominant_haploinsufficient = filtered["symbol"].isin(dominant_haploinsufficient)
    
    variants = filtered[dominant_haploinsufficient].copy()
    
    counts = variants['category'].value_counts()
    counts['missense'] = counts['functional']
    
    return dict(counts)

def print_proportions_from_dominant_hi(dominant_hi, excess):
    ''' show the proportions of loss-of-function/gain of function
    
    These estimates use the split of
    '''
    
    lof_as_mis = dominant_hi['missense']/(dominant_hi['missense'] + dominant_hi['loss-of-function'])
    
    lof_excess = excess['loss-of-function']['excess']/ \
        (excess['missense']['excess'] + excess['loss-of-function']['excess'])
    
    lof_proportion = lof_excess + lof_excess * lof_as_mis
    gof_proportion = 1 - lof_proportion
    
    print('LoF: {:.0f}%'.format(lof_proportion * 100))
    print('GoF: {:.0f}%'.format(gof_proportion * 100))

def main():
    
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations, exclude_synonymous=False)
    known = open_known_genes(args.known_genes)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    de_novos["start_pos"] = de_novos["pos"]
    constraints = pandas.read_table(args.constraints)
    uk_ages = open_uk_parent_ages(args.uk_ages)
    
    # identify the probands with diagnostic de novo variants
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    diagnosed = diagnosed[(diagnosed["inheritance"] == "de_novo") &
        diagnosed["type"].isin(["snv", "indel"])]
    
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    trios = pandas.read_table(args.trios)
    phenotypes = phenotypes[phenotypes["patient_id"].isin(trios["decipher_id"])]
    
    male = 2407
    female = 1886
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    filtered = de_novos[(~de_novos["pp_dnm"].isnull() & (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    
    excess = get_consequence_excess(expected, filtered)
    plot_consequence_excess(excess, "results/excess_by_consequence.pdf")
    plot_excess_by_pp_dnm_threshold(de_novos, expected, increments=100)
    
    dominant_hi = proportion_in_dominant_hi_genes(filtered, known)
    print_proportions_from_dominant_hi(dominant_hi, excess)
    
    in_dominant = count_known_excess(filtered, known)
    print_known_in_excess(in_dominant, excess)
    
    proportions = model_mixing(known, filtered, expected, constraints)
    print('lof proportion'.format(proportions))
    
    excess_de_novos_from_pLI(filtered, expected, constraints)
    plot_proportion_known_by_pLI(filtered, expected,  constraints, known)
    
    prevalence = get_birth_prevalence(male + female, excess,
        cnv_yield=0.1, missing_variants=119.9, enrichment=118.8)
    
    excess_to_lof = (excess['loss-of-function']['excess'] + excess['missense']['excess'])/in_dominant['loss-of-function']
    print(excess_to_lof)
    prevalance_from_rates = check_prevalence_from_baseline_lof(rates, known, excess_to_lof)
    
    plot_prevalence_by_age(prevalance_from_rates, phenotypes, uk_ages, dad_rate=1.53, mom_rate=0.86)
    print(prevalence)
    print(prevalance_from_rates)
    

if __name__ == '__main__':
    main()
