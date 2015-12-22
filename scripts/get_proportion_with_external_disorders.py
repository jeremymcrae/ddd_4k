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

import argparse

import pandas

from ddd_4k.constants import PHENOTYPES, TRIOS
from ddd_4k.load_files import open_phenotypes
from ddd_4k.hpo_matches import find_hpo_matches

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import CalculateSimilarity

def get_options():
    
    parser = argparse.ArgumentParser(description="script to assess how many of" \
        "the DDD cohort have disorders from the external cohorts")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to file containing phenotypic information for the probands.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Paths to table of probands from complete exome-sequenced trios.")
    parser.add_argument("--output",
        default="ddd_proportion_with_external_disorders.tsv",
        help="file to save output to.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    # open the phenotype data, and restrict it to the probands with complete trios
    pheno = open_phenotypes(args.phenotypes)
    trios = pandas.read_table(args.trios, sep="\t")
    pheno = pheno[pheno["patient_id"].isin(trios["decipher_id"])]
    
    # open the HPO ontology, so we can get the set of terms which are relevant to
    # each disorder
    hpo_ontology = Ontology(None)
    hpo_graph = hpo_ontology.get_graph()
    graph = CalculateSimilarity({}, hpo_graph)
    
    # define the root nodes for each disorder
    roots = {"Autism spectrum disorder": ["HP:0000729"],
        "Congenital heart disorder": ["HP:0002564"],
        "Intellectual disability": ["HP:0001249", "HP:0012443", "HP:0100543"],
        "Seizures": ["HP:0001250"],
        "Schizophrenia": ["HP:0100753"]}
    
    columns = ["disorder", "root_term", "count", "proportion"]
    counts = pandas.DataFrame(columns=columns)
    for disorder in sorted(roots):
        # find the HPO terms that are relevant to each disorder
        terms = set([])
        for subterm in roots[disorder]:
            terms |= graph.get_descendants(subterm) | set([subterm])
        term = ",".join(roots[disorder])
        
        # figure out if each proband has a term that is in the disorder set
        matches = find_hpo_matches(pheno["child_hpo"], terms)
        temp = pandas.DataFrame({"disorder": [disorder], "root_term": [term], \
            "count": [sum(matches)], "proportion": [sum(matches)/len(pheno)]})
        counts = counts.append(temp, ignore_index=True)
    
    counts = counts[columns]
    counts.to_csv(args.output, sep="\t", index=False)
