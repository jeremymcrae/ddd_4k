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

import pandas

from ddd_4k.load_files import open_phenotypes, open_families, open_de_novos
from ddd_4k.constants import PHENOTYPES, FAMILIES, TRIOS, SANGER_IDS, DIAGNOSED, \
    VALIDATIONS, DENOVO_PATH, THRESHOLD, SEIZURE_GENES
from ddd_4k.rank_hpo import rank_terms
from ddd_4k.hpo_matches import find_hpo_matches

from hpo_similarity.ontology import Ontology

SEIZURE_RESULTS = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.all.2015-11-24.txt"
# SEIZURE_RESULTS = "/nfs/users/nfs_j/jm33/apps/seizure/results/de_novos.ddd_4k.with_diagnosed.all.2015-10-12.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
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
    parser.add_argument("--results", default=SEIZURE_RESULTS,
        help="Path to table of of results from testing within seizure subset.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    results = pandas.read_table(args.results, sep="\t")
    results = results[results["p_min"] < THRESHOLD]
    
    variants = open_de_novos(args.de_novos, args.validations)
    
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
    
    for gene in results["hgnc"]:
        gene_vars = variants[variants["symbol"] == gene]
        gene_probands = gene_vars["person_stable_id"]
        
        n_probands = float(len(gene_probands))
        n_with_seizures = sum(pheno["has_seizures"][pheno["person_stable_id"].isin(gene_probands)])
        
        if n_probands > 0 and n_probands/n_with_seizures > 0.3 and gene not in SEIZURE_GENES:
            print(gene, n_with_seizures, len(gene_probands), gene in SEIZURE_GENES, list(results["in_ddg2p"][results["hgnc"] == gene]))
    

if __name__ == '__main__':
    main()
