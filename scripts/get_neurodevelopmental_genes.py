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

import argparse

import pandas

from ddd_4k.constants import KNOWN_GENES
from ddd_4k.load_files import open_known_genes
from ddd_4k.hpo_matches import find_hpo_matches

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import CalculateSimilarity

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to identify neuro" \
        "developmental dominant loss-of-function genes in the DDG2P. This " \
        "relies on the hpo_similarity code to identify genes for " \
        "neurodevelopmental HPO terms.")
    parser.add_argument("--known-genes", default=KNOWN_GENES, \
        help="Path to table of known genes, including a column for the target" \
            "organ.")
    parser.add_argument("--output", default="neurodevelopmental.dominant_lof.txt", \
        help="Path to write table of neuodevelopmental genes to.")
    
    args = parser.parse_args()
    
    return args

def get_organ_genes(dominant_lof):
    """ get a set of HGNC symbols for dominant LoF genes that affect the brain
    """
    
    brain_genes = dominant_lof[~dominant_lof["organs"].isnull() & dominant_lof["organs"].str.contains("Brain/Cognition")]
    organ_genes = set(brain_genes["gencode_gene_name"].unique())
    
    return organ_genes

def get_hpo_genes(dominant_lof):
    """ get a set of genes that have brain/cognition phenotypes
    """
    
    # load a graph of HPO terms
    hpo_ontology = Ontology(None)
    graph = hpo_ontology.get_graph()
    
    # get the terms which are descendents of brain/cognition terms
    cognition_root = "HP:0100543"
    brain_root = "HP:0012443"
    cognition_terms = graph.get_descendants(cognition_root) | set([cognition_root])
    brain_terms = graph.get_descendants(brain_root) | set([brain_root])
    
    neurodev_terms = cognition_terms | brain_terms
    
    # find the dominant LoF genes which have at least one brain/cognition term
    matches = find_hpo_matches(dominant_lof["hpo_codes"], neurodev_terms)
    term_genes = set(dominant_lof["gencode_gene_name"][matches].unique())
    
    return term_genes

def main():
    args = get_options()
    
    known = open_known_genes(args.known_genes)
    dominant_lof = known[known["mode"].isin(["Monoallelic", "X-linked dominant"]) & \
        (known["mech"] == "Loss of function")]
    
    organ_genes = get_organ_genes(dominant_lof)
    term_genes = get_hpo_genes(dominant_lof)
    
    # make a dataframe that includes all genes
    genes = pandas.DataFrame({"hgnc": sorted(term_genes | organ_genes)})
    genes["organ"] = genes["hgnc"].isin(organ_genes)
    genes["hpo"] = genes["hgnc"].isin(term_genes)
    
    genes.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
