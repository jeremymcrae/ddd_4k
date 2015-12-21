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

import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot
import pandas
import seaborn

from ddd_4k.constants import TRIOS

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import ICSimilarity
from hpo_similarity.load_files import load_participants_hpo_terms

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--de-novos", \
        help="Path to table of variants in novel genes.")
    parser.add_argument("--phenotypes", \
        default="/lustre/scratch113/projects/ddd/users/jm33/" \
            "de_novos.ddd_4k.phenotypes_by_proband.json", \
        help="JSON file of HPO terms per proband.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Path to table of alternate IDs for participants.")
    
    args = parser.parse_args()
    
    return args

def get_all_terms(graph, probands, proband_hpo):
    """ get the superset of terms used within the probands
    """
    
    all_terms = set()
    for proband in probands:
        terms = proband_hpo[proband]
        
        for term in terms:
            all_terms |= set([term]) | set(graph.get_ancestors(term))
    
    return all_terms

def drop_worse_ranked_terms(table, graph):
    """ drop out terms that are an ancestor of a more highly ranked term
    
    Run through the rows of terms in order of their score, and only include the
    row if that term has not been used before (either the term or its ancestors),
    and if the term immediately above it has not been used either.
    
    Args:
        table: pandas DataFrame of HPO terms, ranked by their "score"
            (calculated as the number of mutated probands with the term,
            multiplied by the information content for the term.)
        graph: ICSimilarity graph object for the HPO ontology.
    
    Returns:
        pandas DataFrame of
    """
    
    used_terms = set()
    to_include = []
    for num, row in table.iterrows():
        term = row["term"]
        no_direct_ancestor = len(set(graph.predecessors(term)) & used_terms) == 0
        not_in_terms = term not in used_terms
        
        if not_in_terms and no_direct_ancestor:
            used_terms |= set([term]) | graph.get_ancestors(term)
            
        to_include.append(not_in_terms & no_direct_ancestor)
    
    return table[to_include]

def rank_terms(graph, probands, proband_hpo):
    """ rank HPO terms used by probands who share a mutation in a gene.
    
    Ranks the HPO terms used probands who share a mutation. This only includes
    terms used by >1 proband, and ranks according to their rarity and how
    many of the probands with mutations share that term.
    
    Args:
        graph: an ICSimilarity object, which has details of how often each term
            was used in the DDD probands.
        probands: list of proband IDs for the probands who share a mutated gene
        proband_hpo: dictionary of HPO terms per proband, indexed by proband ID.
    
    Returns:
        pandas DataFrame of terms, including the term, the name for the term,
        the information content of the term, the count of how many probands with
        mutations share that term, and a score (count x information content).
    """
    
    all_terms = get_all_terms(graph, probands, proband_hpo)
    
    columns = ["term", "name", "IC", "count", "score"]
    table = pandas.DataFrame(columns=columns)
    # rank how many times each term is used, against the rarity
    for term in sorted(all_terms):
        ic = graph.calculate_information_content(term)
        name = graph.node[term]["name"]
        
        # identify all the subterms, and count how many probands have one.
        desc = set([term]) | graph.get_descendants(term)
        count = sum([len(set(proband_hpo[x]) & desc) > 0 for x in probands])
        
        temp = pandas.DataFrame({"term": [term], "name": [name], "IC": [ic],
            "count": [count], "score": [ic * count]})
        table = table.append(temp, ignore_index=True)
    
    # min_count = max(2, floor(len(probands/3)))
    min_count = 2
    # only include terms which are used for more than one proband
    table = table[table["count"] >= min_count]
    
    # sort the columns and rows
    table = table[columns]
    table = table.sort("score", ascending=False)
    table = drop_worse_ranked_terms(table, graph)
    
    return table

def plot_gene(graph, probands, proband_hpo, hgnc, trios, table):
    """ plot the terms relevant for a gene as a heatmap
    
    Args:
        graph: ICSimilarity object, used to look up and down the HPO graph.
        probands: list of proband IDs for the gene.
        proband_hpo: dictionary of HPO terms by proband for the full ddd 4k
            cohort, indexed by proband IDs.
        hgnc: HGNC symbol for the gene, used to determine a pdf filename.
        trios: table of probands in exome-sequenced, including decipher and DDD IDs.
        table: pandas DataFrame of relevant terms for the gene, one row per term.
    """
    
    decipher = trios[trios["proband_stable_id"].isin(probands)]
    decipher_ids = dict(zip(decipher["proband_stable_id"], decipher["decipher_id"]))
    
    columns = sorted(decipher_ids.values())
    data = pandas.DataFrame(columns=columns)
    
    for num, row in table.iterrows():
        term = row["term"]
        name = row["name"]
        count = row["count"]
        
        desc = set([term]) | graph.get_descendants(term)
        has_term = [ len(set(proband_hpo[x]) & desc) > 0 for x in probands ]
        has_term = [ [count] if x else [0] for x in has_term ]
        
        decipher = [ decipher_ids[x] for x in probands ]
        
        temp = pandas.DataFrame(dict(zip(decipher, has_term)), index=[name])
        data = data.append(temp)
    
    data = data[columns]
    data = data[data.columns].astype(float)
    
    # plot the heatmap
    ax = seaborn.heatmap(data, cmap="Blues", cbar=False, square=True)
    fig = ax.get_figure()
    fig.savefig("{}.pdf".format(hgnc), format="pdf")
    pyplot.close()

def main():
    args = get_options()
    
    args.de_novos = "/lustre/scratch113/projects/ddd/users/jm33/results/novel_gene_variants.ddd_4k.2015-11-24.txt"
    variants = pandas.read_table(args.de_novos, sep="\t")
    
    # open the phenotype data, and restrict it to the probands with complete trios
    trios = pandas.read_table(args.trios, sep="\t")
    proband_ids = set(trios["proband_stable_id"])
    
    # open the HPO ontology, so we can get the set of terms which are relevant
    # to each disorder
    hpo_ontology = Ontology(None)
    graph = hpo_ontology.get_graph()
    
    alt_node_ids = hpo_ontology.get_alt_ids()
    obsolete_ids = hpo_ontology.get_obsolete_ids()
    
    # load the HPO terms for the probands, and subset the probands to the ones
    # in the DDD 4K trios
    proband_hpo = load_participants_hpo_terms(args.phenotypes, \
        alt_node_ids, obsolete_ids)
    proband_hpo = { x: proband_hpo[x] for x in proband_hpo if x in proband_ids }
    
    graph.tally_hpo_terms(proband_hpo)
    
    for hgnc in sorted(variants["symbol"].unique()):
        gene = variants[variants["symbol"] == hgnc]
        probands = list(gene["person_stable_id"])
        table = rank_terms(graph, probands, proband_hpo)
        
        if len(table) > 0:
            plot_gene(graph, probands, proband_hpo, hgnc, trios, table)

if __name__ == '__main__':
    main()
