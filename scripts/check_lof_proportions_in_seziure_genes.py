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

import os
import argparse

from scipy.stats import fisher_exact
import pandas

from ddd_4k.load_files import open_phenotypes, open_families, open_de_novos
from ddd_4k.constants import PHENOTYPES, FAMILIES, TRIOS, SANGER_IDS, DIAGNOSED, \
    VALIDATIONS, DENOVO_PATH, SEIZURE_GENES
from ddd_4k.rank_hpo import rank_terms
from ddd_4k.hpo_matches import find_hpo_matches

from hpo_similarity.ontology import Ontology

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to test for a" \
        "difference in the proportion of probands with different functional" \
        "categories of de novos. This compares the proportion in genes known" \
        "to be associated with seizures, and looks for differences in" \
        "probands who have been recorded as having seizures versus probands" \
        "without any record of seizures. The hypothesis is that probands with" \
        "a de novo in a seizure gene, but without any record of seizures" \
        "might have a difference in how severe their mutations are compared" \
        "to probands who do have a record of seizures.")
    parser.add_argument("--de-novos", default=DENOVO_PATH, \
        help="Path to table of variants in novel genes.")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to table of phenotypic data from probands.")
    parser.add_argument("--validations", default=VALIDATIONS, \
        help="Path to table of results from de novo validations.")
    parser.add_argument("--families", default=FAMILIES,
        help="Path to table of DDD families.")
    parser.add_argument("--trios", default=TRIOS,
        help="Path to table of DDD trios.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS,
        help="Path to table of mapping DDD IDs to decipher IDs.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    variants = open_de_novos(args.de_novos, args.validations)
    
    # diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    # variants = variants[~variants["person_stable_id"].isin(diagnosed["person_id"])]
    
    # open the phenotype data, and restrict it to the probands with complete trios
    pheno = open_phenotypes(args.phenotypes, args.sanger_ids)
    trios = pandas.read_table(args.trios, sep="\t")
    proband_ids = set(trios["proband_stable_id"])
    pheno = pheno[pheno["person_stable_id"].isin(proband_ids)]
    
    # open the HPO ontology, so we can get the set of terms which are relevant
    # to each disorder
    hpo_ontology = Ontology(None)
    graph = hpo_ontology.get_graph()
    
    seizure_root = "HP:0001250"
    seizure_terms = graph.get_descendants(seizure_root) | set([seizure_root])
    
    pheno["has_seizures"] = find_hpo_matches(pheno["child_hpo"], seizure_terms)
    
    variants = variants[variants["symbol"].isin(SEIZURE_GENES)]
    
    seizure_variants = variants[["person_stable_id", "sex", "chrom", "pos", "ref", "alt",
        "symbol", "category", "consequence"]].merge(
        pheno[["person_stable_id", "has_seizures"]], how="left",
        on=["person_stable_id"])
    
    has = seizure_variants[seizure_variants["has_seizures"]]["category"].value_counts()
    hasnt = seizure_variants[~seizure_variants["has_seizures"]]["category"].value_counts()
    
    print(fisher_exact([has, hasnt]))

if __name__ == '__main__':
    main()
