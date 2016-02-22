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

from __future__ import absolute_import

import argparse

import pandas

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, \
    CONSTRAINTS_URL, PHENOTYPES, TRIOS, SANGER_IDS, DIAGNOSED
from ddd_4k.causation.excess_by_consequence import get_consequence_excess, \
    plot_consequence_excess
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

def main():
    
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations, exclude_synonymous=False)
    known = open_known_genes(args.known_genes)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    de_novos["start_pos"] = de_novos["pos"]
    constraints = pandas.read_table(args.constraints)
    uk_ages = open_uk_parent_ages(args.uk_ages)
    
    # identiry the probands with diagnostic de novo variants
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    diagnosed = diagnosed[diagnosed["inheritance"] == "de_novo"]
    diagnosed = diagnosed["person_id"].unique()
    
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    trios = pandas.read_table(args.trios)
    phenotypes = phenotypes[phenotypes["patient_id"].isin(trios["decipher_id"])]
    
    male = 2407
    female = 1887
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    pp_dnm_threshold = get_pp_dnm_threshold(de_novos, expected)
    filtered = de_novos[(~de_novos["pp_dnm"].isnull() & (de_novos["pp_dnm"] > pp_dnm_threshold)) ]
    
    excess = get_consequence_excess(expected, filtered)
    plot_consequence_excess(excess, "results/excess_by_consequence.pdf")
    
    proportions = model_mixing(known, filtered, expected, constraints)
    print(proportions)
    
    excess_de_novos_from_pLI(filtered, expected, constraints)
    plot_proportion_known_by_pLI(filtered, expected,  constraints, known)
    
    prevalence = get_birth_prevalence(male + female, excess,
        cnv_yield=0.1, missing_variants=119.9, enrichment=118.8)
    
    plot_prevalence_by_age(prevalence, phenotypes, diagnosed, uk_ages,
        mutations_per_year=2.5)
        
    prevalance_from_rates = check_prevalence_from_baseline_lof(rates, known,
        mis_to_lof=2.0, missing=0.5)
    
    print(prevalance_from_rates)
    

if __name__ == '__main__':
    main()
