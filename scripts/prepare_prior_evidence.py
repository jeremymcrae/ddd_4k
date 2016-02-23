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
import inflect

from ddd_4k.constants import THRESHOLD
from ddd_4k.convert_doi import doi_to_pubmed
from ddd_4k.combine_enrichment_and_clustering import get_gene_results

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to determine the prior"
        "evidence for candidate novel genes.")
    parser.add_argument("--external-enrichment", \
        help="Path to table of results from testing for de novo enrichment in"
            "the external studies alone.")
    parser.add_argument("--external-clustering", \
        help="Path to table of results from testing for de novo clustering in"
            "the external studies alone.")
    parser.add_argument("--external-de-novos", \
        help="Path to table of results from testing for de novo clustering in"
            "the external studies alone.")
    parser.add_argument("--subsets", default="intellectual_disability,epilepsy,autism,normal_iq_autism", \
        help="comma-separated list of phenotypes to pull external de novos for.")
    parser.add_argument("--results", \
        help="Path to table of results from de novo testing.")
    parser.add_argument("--output", default="prior_evidence.txt", \
        help="Path to write gene sentences to.")
    
    args = parser.parse_args()
    
    return args

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
    
    # figure out the pubmed IDs for all the variants from the doi codes
    dois = variants["publication_doi"].unique()
    pubmed_ids = [ doi_to_pubmed(x) for x in dois ]
    recode = dict(zip(dois, pubmed_ids))
    variants["pubmed"] = variants["publication_doi"].map(recode)
    
    return variants

def get_novel_genes(results_path):
    """ get the HGNC symbols for the novel genes with genomewide significance
    
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
        pubmed_ids: list (or numpy array) of pubmed IDs
    
    Returns:
        correctly formatted stext string to insert in a sentence listing the
        pubmed IDs.
    """
    
    pubmed_ids = sorted(set(pubmed_ids))
    
    if len(pubmed_ids) == 0:
        pubmed = ""
    elif len(pubmed_ids) == 1:
        pubmed = " (PMID {})".format(pubmed_ids[0])
    else:
        final = pubmed_ids.pop()
        pubmed_ids = ", ".join(pubmed_ids)
        pubmed = " (PMID {} and {})".format(pubmed_ids, final)
        
    return pubmed

def initial_sentence(external, gene, pubmed_ids, numbers):
    """ state how many de novos are in external cohorts (including pubmed IDs).
    
    Args:
        external: pandas Series of values for a gene. The Series can be empty if
            there aren't any de novos in that gene in the external cohorts.
        gene: hgnc symbol for the gene of interest
        pubmed_ids: list of pubmed IDs for the studies that the variants are
            found in.
        numbers: inflect engine, to convert python numbers into their word
            equivalent, e.g. 1 -> "one".
    
    Returns:
        text sentence stating the number of variants found in external cohorts.
    """
    
    pubmed_text = format_pubmed(pubmed_ids)
    
    dnms = 0
    if len(external) > 0:
        dnms = external["dnms"]
    
    # figure out the correct type for the word "time", e.g. 1 "time", 2 "times"
    times = "times"
    if len(external) > 0 and external["dnms"] == 1:
        times = "time"
    
    initial_text = "In other large scale exome sequencing projects de novo " \
        "mutations in {} have been identified {} {}{}.".format(gene,\
        numbers.number_to_words(dnms), times, pubmed_text)
    
    return initial_text

def counts_sentence(external, numbers):
    """ state the count of de novos from external cohorts by consequence type.
    
    Args:
        external: pandas Series of values for a gene. The Series can be empty if
            there aren't any de novos in that gene in the external cohorts.
        numbers: inflect engine, to convert python numbers into their word
            equivalent, e.g. 1 -> "one".
    
    Returns:
        text sentence stating the number of loss-of-function and
        missense/inframe de novos found in external cohorts.
    """
    
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
            be negligible. The Series can be empty if there aren't any de novos
            in that gene in the external cohorts.
    
    Returns:
        sentence stating strength of evidence e.g. "strong", "suggestive", or
        "neglible".
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
    """ state the evidence for clustering of de novos from external cohorts.
    
    Args:
        external: pandas Series of information for a gene. This includes a
            boolean entries for if the gene "has_clusering". The Series can be
            empty if there aren't any de novos in that gene in the external
            cohorts.
    
    Returns:
        text sentence stating the evidence for de novo clustering from de novos
        found in external cohorts.
    """
    
    clustering = "is not"
    if len(external) > 0 and external["has_clustering"]:
        clustering = "is"
    
    clust_text = "There {} evidence for clustering of mutations in the gene " \
        "or encoded protein.".format(clustering)
    
    return clust_text

def prepare_text(gene, enrich, variants, numbers):
    """ state the prior evidence for a gene in a few sentences
    
    Args:
        gene: hgnc symbol for a gene
        enrich: pandas dataframe of enrichment results from the external cohorts
        variants: pandas dataframe of variants included in the testing of the
            external cohorts.
        numbers: inflect engine, to convert python numbers into their word
            equivalent, e.g. 1 -> "one".
    
    Returns:
        text for a few sentences stating the prior evidence for a gene.
    """
    
    external = enrich[enrich["hgnc"] == gene]
    sites = variants[variants["hgnc"] == gene]
    
    # make sure the variant count matches the analyzed count
    if gene not in list(enrich["hgnc"]):
        assert len(sites) == 0
    else:
        assert list(external["dnms"])[0] == len(sites)
    
    if len(external) > 1:
        text = "too many possible rows for {}".format(gene)
    else:
        external = external.squeeze()
        initial_text = initial_sentence(external, gene, sites["pubmed"], numbers)
        counts_text = counts_sentence(external, numbers)
        clust_text = clustering_sentence(external)
        evidence_text = strength_sentence(external)
        
        if len(external) > 0 and external["dnms"] > 0:
            text = initial_text + " " + counts_text + " " + evidence_text+ " " + clust_text
        else:
            text = initial_text + " " + evidence_text
    
    return text

def main():
    args = get_options()
    
    # load the de novo variants from the external cohorts, and the combined
    # results from testing for enrichment and clustering in those datasets
    variants = get_external_variants(args.external_de_novos, args.subsets)
    enrich = get_gene_results(args.external_enrichment, args.external_clustering)
    
    numbers = inflect.engine()
    output = open(args.output, "w")
    for gene in get_novel_genes(args.results):
        text = prepare_text(gene, enrich, variants, numbers)
        
        output.write(gene + "\n")
        output.write(text + "\n")
        output.write("\n")
    
    output.close()

if __name__ == '__main__':
    main()
