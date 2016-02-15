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

from __future__ import division

from math import ceil

import pandas

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
    
    min_count = max(2, ceil(len(probands)/4.0))
    # min_count = 2
    # only include terms which are used for more than one proband
    table = table[table["count"] >= min_count]
    
    # sort the columns and rows
    table = table[columns]
    table = table.sort("score", ascending=False)
    table = drop_worse_ranked_terms(table, graph)
    
    return table
