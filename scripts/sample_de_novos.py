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

import pandas

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_gene import get_transcript_lengths, construct_gene_object, \
    get_transcript_ids_sorted_by_length
from denovonear.site_specific_rates import SiteRates
from denovonear.weighted_choice import WeightedChoice

from ddd_4k.get_hgnc_symbols import get_all_hgnc_symbols

mut_path = "/nfs/users/nfs_j/jm33/apps/denovonear/data/forSanger_1KG_mutation_rate_table.txt"

def get_transcript_for_gene(symbol, cache_dir="cache", genome_build="grch37"):
    """ identify and create a transcript object for a gene symbol
    
    Args:
        symbol: HGNC symbol for a gene
        cache_dir: path to folder for caching ensembl request information
        genome_build: genome build to request information for (eg "grch37")
    
    Returns:
        denovonear Transcript object for transcript, which contains coordinates,
        sequence, and methods to transform around the transcript.
    """
    
    ensembl = EnsemblRequest(cache_dir, genome_build)
    
    transcript_ids = get_transcript_ids_sorted_by_length(ensembl, symbol)
    
    # work through the transcript IDs in descending lengths
    transcript = None
    for (transcript_id, length) in transcript_ids:
        try:
            transcript = construct_gene_object(ensembl, transcript_id)
            break
        except ValueError:
            continue
    
    if transcript is None:
        raise IndexError("no suitable transcript for {}".format(symbol))
    
    return transcript

def get_rates_for_gene(gene_id, all_genes, mut_dict):
    """ find the missense and lof mutation rates for a gene
    
    For a given HGNC symbol, we want to find the longest complete protein coding
    transcript for that gene (some genes lack such transcripts, for example RNA
    genes, we ignore those genes). We want to be able to sample sites within the
    gene according to the null mutation model.
    
    Args:
        gene_id: HGNC symbol for a gene
        all_genes: dictionary of transcripts, mutation rates and objects to
            randomly sample functional types and sites within the transcript,
            indexed by HGNC symbols.
        mut_dict: dictionary of tri-nucleotide mutation rates for the null
            mutation model.
    
    Returns:
         the all_genes dictionary, likely with an extra entry for the current
         gene, if the gene has a complete protein coding transcript.
    """
    
    try:
        transcript = get_transcript_for_gene(gene_id)
    except IndexError:
        return all_genes
    
    # get the site specific mutation rates for missense and lof possibilities in
    # the gene
    site_weights = SiteRates(transcript, mut_dict)
    missense = site_weights.get_missense_rates_for_gene()
    lof = site_weights.get_lof_rates_for_gene()
    
    choices = [("missense", missense.cum_probs[-1]), ("lof", lof.cum_probs[-1])]
    gene_choices = WeightedChoice(choices)
    
    # add an entry to the all genes dictionary for the current gene
    all_genes[gene_id] = {}
    all_genes[gene_id]["transcript"] = transcript
    all_genes[gene_id]["missense"] = missense
    all_genes[gene_id]["lof"] = lof
    all_genes[gene_id]["rate"] = missense.cum_probs[-1] + lof.cum_probs[-1]
    all_genes[gene_id]["choices"] = gene_choices
    
    return all_genes

def exclude_readthrough_genes(symbols):
    """ some gene symbols are present because they are part of readthrough genes
    
    Rather than possibly sampling the same region twice (as part of the sole
    gene, and then as part of the readthrough gene), we will exclude the symbols
    where both parts exist separately in the list of symbols.
    
    Args:
        symbols: list of HGNC symbols
    """
    
    modified = []
    for x in symbols:
        if "-" in x:
            if sum([ i in symbols for i in x.split("-") ]) > 1:
                continue
            modified.append(x)
        else:
            modified.append(x)
    
    return modified

def sample_de_novos(gene_sampler, all_genes):
    """
    """
    
    proportion_covered = 0.9
    
    columns = ["person_id", "chrom", "pos", "ref", "alt", "hgnc", "consequence", "type"]
    de_novos = pandas.DataFrame(columns=columns)
    
    iteration = 1
    while True:
        person_id = "person{:0>8}".format(iteration)
        
        # sample a gene, where the chance of sampling each gene is proportional
        # to the mutation rate for the genes
        sampled_hgnc = gene_sampler.choice()
        
        # sample a consequence within the gene, where the chance of sampling
        # each functional type is proportional to the mutation rate for each
        # functional type
        consequence = all_genes[sampled_hgnc]["choices"].choice()
        
        # sample a position within the gene for the functional type
        cds_pos = all_genes[sampled_hgnc][consequence].choice()
        bp = all_genes[sampled_hgnc]["transcript"].get_position_on_chrom(cds_pos)
        chrom = all_genes[sampled_hgnc]["transcript"].get_chrom()
        
        ref = all_genes[sampled_hgnc]["transcript"].get_trinucleotide(bp)[1]
        alt = "X"
        var_type = "snv"
        
        de_novos = de_novos.append({"person_id": person_id, "chrom": chrom,
            "pos": bp, "ref": ref, "alt": alt, "hgnc": sampled_hgnc,
            "consequence": consequence, "type": var_type}, ignore_index=True)
        
        iteration += 1
        if iteration % 1000 == 0:
            print("{} de novos samples".format(iteration))
            if len(de_novos["hgnc"].unique())/len(all_genes) > proportion_covered:
                break
    
    return de_novos

def main():
    mut_dict = load_mutation_rates(mut_path)
    
    # print("extracting HGNC symbols from VCFs")
    # symbols = get_all_hgnc_symbols()
    # symbols = exclude_readthrough_genes(symbols)
    
    symbols = ["AES"]
    print("getting transcript information for genes")
    all_genes = {}
    for gene_id in symbols:
        print(gene_id)
        all_genes = get_rates_for_gene(gene_id, all_genes, mut_dict)
    
    rates = [ all_genes[x]["rate"] for x in all_genes ]
    gene_sampler = WeightedChoice(zip(all_genes.keys(), rates))
    
    print("randomly sampling de novos within genes")
    de_novos = sample_de_novos(gene_sampler, all_genes)
    de_novos.to_csv("sampled_de_novos.txt", sep="\t", index=False)

if __name__ == '__main__':
    main()
