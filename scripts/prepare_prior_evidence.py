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
from statsmodels.stats.multitest import fdrcorrection
import inflect

from ddd_4k.constants import ALPHA, NUM_GENES, THRESHOLD
from ddd_4k.convert_doi import doi_to_pubmed
from ddd_4k.combine_p_values import fishers_method

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to determine the prior"
        "evidence for candidate novel genes.")
    parser.add_argument("--external-enrichment", default=ENRICH, \
        help="Path to table of results from testing for de novo enrichment in"
            "the external studies alone.")
    parser.add_argument("--external-clustering", default=CLUSTER, \
        help="Path to table of results from testing for de novo clustering in"
            "the external studies alone.")
    parser.add_argument("--external-de-novos", default=EXTERNAL_DE_NOVOS, \
        help="Path to table of results from testing for de novo clustering in"
            "the external studies alone.")
    parser.add_argument("--subsets", default="intellectual_disability,epilepsy,autism,normal_iq_autism", \
        help="comma-separated list of phenotypes to pull external de novos for.")
    parser.add_argument("--results", default=RESULTS, \
        help="Path to table of results from de novo testing.")
    parser.add_argument("--output", default="prior_evidence.txt", \
        help="Path to write gene sentences to.")
    
    args = parser.parse_args()
    
    return args

def open_enrichment(path):
    """ load the enrichment results, and count how many loss-of-function and
    non-loss-of-function, but atill protein altering de novos the external
    studies had have per gene
    """
    
    enrich = pandas.read_table(path, sep="\t")
    enrich["mis"] = enrich[["missense_indel", "missense_snv"]].sum(axis=1)
    enrich["lof"] = enrich[["lof_indel", "lof_snv"]].sum(axis=1)
    enrich["dnms"] = enrich[["mis", "lof"]].sum(axis=1)
    
    return enrich

def open_clustering(path):
    """ load the clustering results, and prepare the table for merging with the
    enrichment results, i.e. rename columns and selct only the columns to be
    merged
    """
    
    cluster = pandas.read_table(path, sep="\t")
    cluster = cluster[cluster["mutation_category"] == "missense"]
    cluster["p_clust"] = cluster["probability"]
    cluster["hgnc"] = cluster["gene_id"]
    cluster = cluster[["hgnc", "p_clust"]]
    
    return cluster

def get_external_variants(path, subsets):
    """ load the de novos from the external cohorts
    """
    
    subsets = subsets.split(",")
    
    functional = ["missense_variant", "frameshift_variant", "stop_gained",
        "splice_donor_variant", "splice_acceptor_variant", "inframe_deletion",
        "inframe_insertion", "start_lost", "stop_lost",
        "protein_altering_variant", "stop_retained_variant",
        "coding_sequence_variant"]
    
    variants = pandas.read_table(path, sep="\t", compression="gzip")
    variants = variants[variants["consequence"].isin(functional)]
    variants = variants[variants["study_phenotype"].isin(subsets)]
    
    return variants

def get_novel_genes(results_path):
    """ get the HGNC symbols for the genes with genomewide significance
    
    Args:
        results_path: path to table of results, one row for each gene
    
    Returns:
        list of HGNC symbols for the candidate novel genes.
    """
    
    results = pandas.read_table(results_path, sep="\t")
    new_genes = results["hgnc"][results["p_min"] < THRESHOLD]
    new_genes = new_genes[~new_genes.isnull()]
    
    return new_genes

def format_pubmed(pubmed_ids):
    """ format a list of pubmed IDs into a text string suitable for a sentence
    
    Args:
        pubmed_ids: list of pubmed IDs
    
    Returns:
        correctly formatted stext string to insert in a sentence listing the
        pubmed IDs.
    """
    
    pubmed_ids = sorted(pubmed_ids)
    
    if len(pubmed_ids) == 0:
        pubmed = ""
    elif len(pubmed_ids) == 1:
        pubmed = " (PMID {})".format(pubmed_ids[0])
    else:
        final = pubmed_ids.pop()
        pubmed_ids = ", ".join(pubmed_ids)
        pubmed = " (PMID {} and {})".format(pubmed_ids, final)
        
    return pubmed

def initial_sentence(external, gene, pubmed, numbers):
    
    dnms = 0
    if len(external) > 0:
        dnms = external["dnms"]
    
    # figure out the correct type for the word "time", e.g. 1 "time", 2 "times"
    times = "times"
    if len(external) > 0 and external["dnms"] == 1:
        times = "time"
    
    initial_text = "In other large scale exome sequencing projects de novo " \
        "mutations in {} have been identified {} {}{}.".format(gene,\
        numbers.number_to_words(dnms), times, pubmed)
    
    return initial_text

def counts_sentence(external, numbers):
    
    lof_count = 0
    mis_count = 0
    dnms = 0
    if len(external) > 0:
        dnms = external["dnms"]
        lof_count = external["lof"]
        mis_count = external["mis"]
    
    if dnms == 1:
        if lof_count == 1:
            return "This mutation was loss of function."
        elif mis_count == 1:
            return "This mutation was missense/inframe."
    
    lof = "were"
    if len(external) > 0 and external["lof"] == 1:
        lof = "was"
    
    mis = "were"
    if len(external) > 0 and external["mis"] == 1:
        mis = "was"
    
    numbers_text = "{} of these {} missense/inframe and {} {} loss of " \
        "function.".format(\
            numbers.number_to_words(external["mis"]).capitalize(), mis, \
            numbers.number_to_words(external["lof"]), lof)
    
    return numbers_text

def strength_sentence(external):
    """ get the strength of evidence for a gene as a sentence
    
    Args:
        external: pandas Series of information for a gene. This includes entries
            for if the gene is "genomewide" significant, or "suggestive". If
            neither of these entries are true, then we consider the evidence to
            be negligible.
    
    Returns:
        strength of evidence in a sentence e.g. "strong", "suggestive", or
        "neglible"
    """
    
    if len(external) == 0:
        strength = "negligible"
    elif external["genomewide"]:
        strength = "strong"
    elif external["suggestive"]:
        strength = "suggestive"
    else:
        strength = "negligible"
    
    strength_text = "Based on this we would assess the prior evidence for "\
        "pathogenesis to be {}.".format(strength)
        
    return strength_text

def clustering_sentence(external):
    
    clustering = "is not"
    if len(external) > 0 and external["has_clustering"]:
        clustering = "is"
    
    clust_text = "There {} evidence for clustering of mutations in the gene " \
        "or encoded protein.".format(clustering)
    
    return clust_text

def main():
    args = get_options()
    
    numbers = inflect.engine()
    
    variants = get_external_variants(args.external_de_novos, args.subsets)
    novel_genes = get_novel_genes(args.results)
    enrich = open_enrichment(args.external_enrichment)
    cluster = open_clustering(args.external_clustering)
    output = open(args.output, "w")
    
    # figure out the pubmed IDs for all the variants from the doi codes
    dois = variants["publication_doi"].unique()
    pubmeds = [ doi_to_pubmed(x) for x in dois ]
    
    recode = dict(zip(dois, pubmeds))
    variants["pubmed"] = variants["publication_doi"].map(recode)
    
    # combine the enrichment and clustering results
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
    
    for gene in novel_genes:
        external = enrich[enrich["hgnc"] == gene]
        sites = variants[variants["hgnc"] == gene]
        
        # make sure the variant count matches the analyzed count
        if gene not in list(enrich["hgnc"]):
            assert len(sites) == 0
        else:
            assert list(external["dnms"])[0] == len(sites)
        
        pubmed = format_pubmed(sites["pubmed"].unique())
        
        if len(external) > 1:
            text = "too many possible rows for {}".format(gene)
        else:
            external = external.squeeze()
            initial_text = initial_sentence(external, gene, pubmed, numbers)
            counts_text = counts_sentence(external, numbers)
            clust_text = clustering_sentence(external)
            evidence_text = strength_sentence(external)
            
            if len(external) > 0 and external["dnms"] > 0:
                text = initial_text + " " + counts_text + " " + evidence_text+ " " + clust_text
            else:
                text = initial_text + " " + evidence_text
        
        output.write(gene + "\n")
        output.write(text + "\n")
        output.write("\n")
    
    output.close()

if __name__ == '__main__':
    main()
