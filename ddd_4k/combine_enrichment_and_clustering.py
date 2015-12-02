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

import pandas
from statsmodels.stats.multitest import fdrcorrection

from ddd_4k.constants import NUM_GENES, ALPHA
from ddd_4k.combine_p_values import fishers_method

def get_gene_results(enrich_path, cluster_path):
    """ combine the enrichment and clustering results
    """
    
    enrich = open_enrichment(enrich_path)
    cluster = open_clustering(cluster_path)
    enrich = enrich.merge(cluster, how="left", on="hgnc")
    enrich["p_combined"] = enrich[["p_func", "p_clust"]].apply(fishers_method, axis=1)
    
    # select the best performing model for each gene
    enrich["p_min"] = enrich[["p_lof", "p_combined"]].min(axis=1)
    enrich = enrich[~enrich["p_min"].isnull()]
    enrich = enrich.sort("p_min")
    enrich["fdr"] = fdrcorrection(enrich["p_min"])[1]
    
    # set a genomewide threshold for the external studies, based on the number
    # genes in the genome.
    num_tests = 2 * NUM_GENES
    threshold = ALPHA/num_tests
    
    enrich["genomewide"] = enrich["p_min"] < threshold
    enrich["suggestive"] = enrich["fdr"] < ALPHA
    enrich["has_clustering"] = enrich["p_clust"] < ALPHA
    
    return enrich

def open_enrichment(path):
    """ load the enrichment results
    
    Args:
        path: path to table of de novo enrichment testing results.
    
    Returns:
        pandas DataFrame of enrichment results
    """
    
    enrich = pandas.read_table(path, sep="\t")
    
    # count how many loss-of-function and non-loss-of-function, but still
    # protein altering de novos the external studies had have per gene
    enrich["mis"] = enrich[["missense_indel", "missense_snv"]].sum(axis=1)
    enrich["lof"] = enrich[["lof_indel", "lof_snv"]].sum(axis=1)
    enrich["dnms"] = enrich[["mis", "lof"]].sum(axis=1)
    
    return enrich

def open_clustering(path):
    """ load the clustering results
    
    Args:
        path: path to table of de novo clustering testing results.
    
    Returns:
        pandas DataFrame of clustering results
    """
    
    cluster = pandas.read_table(path, sep="\t")
    
    # prepare the table for merging with the enrichment results, i.e. rename
    # columns and selct only the columns to be merged
    cluster = cluster[cluster["mutation_category"] == "missense"]
    cluster["p_clust"] = cluster["probability"]
    cluster["hgnc"] = cluster["gene_id"]
    cluster = cluster[["hgnc", "p_clust"]]
    
    return cluster
