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

import os
import argparse

from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import get_transcript_ids_sorted_by_length, \
    construct_gene_object

import pandas
from statsmodels.stats.multitest import fdrcorrection

from ddd_4k.constants import THRESHOLD

CNV_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/" \
    "ddd_4k.de_novo_cnvs.2015-10-12.txt"
ASSOCIATIONS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/" \
    "de_novos.ddd_4k.without_diagnosed.all.2015-11-24.txt"
DDG2P_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/" \
    "2015-04-13/DDG2P/dd_genes_for_clinical_filter"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check whether"
        "candidate de novo CNVs overlap candidate novel genes")
    parser.add_argument("--cnvs", default=CNV_PATH, \
        help="Path to table of candidate CNVs.")
    parser.add_argument("--associations", default=ASSOCIATIONS_PATH, \
        help="Path to table of association results.")
    parser.add_argument("--ddg2p", default=DDG2P_PATH, \
        help="Path to table of DDG2P genes.")
    parser.add_argument("--output", default="cnvs_in_interesting_genes.txt", \
        help="Path to plot graph to.")
    
    args = parser.parse_args()
    
    return args

def open_associations(path):
    """ load a file of association data
    
    Args:
        path: path to tab-sepaarred file containing combined association data
    
    Returns:
        pandas DataFrame of association data, with an extra column for
        fdr-corrected p-values.
    """
    
    results = pandas.read_table(path, sep="\t")
    results = results.sort(columns=["p_min"])
    results = results[~results["p_min"].isnull()]
    
    results["fdr"] = fdrcorrection(list(results["p_min"]))[1]
    results["genomewide"] = results["p_min"] < THRESHOLD
    
    return results

def open_cnvs(path):
    """ Open table of canddidate de novo CNVs
    
    Args:
        path: path to tab-separated file containing CNVs per proband
    
    Returns:
        pandas DataFrame of CNV table
    """
    
    cnvs = pandas.read_table(path, sep="\t")
    cnvs["start"] = cnvs["pos"]
    cnvs["end"] = cnvs["start"] + cnvs["SVLEN"]
    
    return cnvs

def open_ddg2p(path):
    """ open table of DDG2P genes
    
    Args:
        path: path to DDG2P data
    
    Returns:
        pandas DataFrame of DDG2P genes
    """
    
    ddg2p = pandas.read_table(path, sep="\t", skiprows=True, \
        names=["chrom", "start", "stop", "gene", "type", "mode", "mech", \
            "disorder"])
    
    return ddg2p

def get_transcript_for_gene(symbol, ensembl):
    """ obtain a Transcript object for a gene from ensembl coordinates
    
    Args:
        symbol: HGNC symbol for a gene
        ensembl: EnsemblRequest object, used to obtain information about the
            gene from Ensembl
    
    Returns:
        Transcript object for the gene, for the longest protein coding complete
        transcript for the gene. Or returns None if no such transcript can be
        found.
    """
    
    transcript_ids = get_transcript_ids_sorted_by_length(ensembl, symbol)
    
    # work through the transcript IDs in descending lengths
    transcript = None
    for transcript_id in sorted(transcript_ids, key=lambda k: -transcript_ids[k]):
        try:
            transcript = construct_gene_object(ensembl, transcript_id)
            break
        except ValueError:
            continue
    
    return transcript

def get_overlapping(cnvs, transcript):
    """ find the candidate de novos CNVs that overlap a given transcript
    
    Args:
        cnvs: pandas DataFrame of all candidate de novo CNVs
        transcript: Transcript object for a gene, containing gene, exon and CDS
            coordinates, along with gene sequences.
    
    Returns:
        DataFrame of CNVs that overlap the given transcript, or empty DataFrame
        if there aren't any overlapping CNVs.
    """
    
    # get the genome range for the CDS
    cds_range = (transcript.get_cds_start(), transcript.get_cds_end())
    lo = min(cds_range)
    hi = max(cds_range)
    
    chrom = transcript.get_chrom()
    
    overlap = cnvs[(cnvs["chrom"] == chrom) & (cnvs["start"] <= hi) & (cnvs["end"] >= lo)]
    
    return overlap

def filter_cnvs(cnvs, ddg2p):
    """ filter out canddiate de novo CNVs that are lower quality, in
    
    Args:
        cnvs: DataFrame of CNVs found to be overlapping a set of genes
        ddg2p: DataFrame of DDG2P genes
    
    Returns:
        dataframe of CNVs, with some CNVs filtered out, due to being lower
        quality, recurrent across different genes, in DDG2P genes, likely to be
        inherited (according to CIFER), or in probands with numerous candidate
        CNVs.
    """
    
    cnvs["madr"] = abs(cnvs["MEANLR2"]/cnvs["MADL2R"])
    
    # only include high confidence CNVs
    cnvs = cnvs[cnvs["madr"] > 15]
    
    # some CNVs occur multiple times for different genes. Strip out the duplicates
    cnvs = cnvs[~cnvs[["person_id", "chrom", "start"]].duplicated()]
    
    # only include the CNVs that are not for DDG2P genes
    cnvs = cnvs[~cnvs["associated_hgnc"].isin(ddg2p["gene"])]
    
    # drop out the CNVs with cifer inheritance classification of being inherited
    cnvs = cnvs[~cnvs["cifer"].isin(["maternal_inh", "paternal_inh"])]
    
    # drop out the CNVs that are in people with multiple candidate de novo CNVs
    dup_persons = cnvs["person_id"][cnvs["person_id"].duplicated()].value_counts()
    cnvs = cnvs[~cnvs["person_id"].isin(dup_persons.index[dup_persons > 5])]
    
    return cnvs

def main():
    """
    """
    
    args = get_options()
    
    ensembl = EnsemblRequest(cache_folder="cache", genome_build="grch37")
    cnvs = open_cnvs(args.cnvs)
    results = open_associations(args.associations)
    ddg2p = open_ddg2p(args.ddg2p)
    
    results = results[results["fdr"] < 0.05]
    
    overlaps = pandas.DataFrame(columns=list(cnvs.columns))
    
    for (pos, row) in results.iterrows():
        transcript = get_transcript_for_gene(row["hgnc"], ensembl)
        
        if transcript is None:
            continue
        
        de_novos = get_overlapping(cnvs, transcript)
        de_novos.loc[:, "associated_hgnc"] = row["hgnc"]
        de_novos.loc[:, "p_value"] = row["p_min"]
        de_novos.loc[:, "fdr_value"] = row["fdr"]
        de_novos.loc[:, "genomewide"] = row["genomewide"]
        
        overlaps = overlaps.append(de_novos, ignore_index=True)
    
    overlaps = filter_cnvs(overlaps, ddg2p)
    overlaps.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
