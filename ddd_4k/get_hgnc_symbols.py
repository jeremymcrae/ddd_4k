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

import sys
import gzip
import re

import pandas

from ddd_4k.constants import FAMILIES, TRIOS
from ddd_4k.parse_vcf import parse_vcf_format, parse_vcf_info, get_vcf_header, \
    exclude_vcf_header

IS_PYTHON3 = sys.version_info[0] == 3

def get_symbols_from_line(line):
    """ pulls out the HGNC symbols from the info field of a VCF line
    
    Args:
        line: string for line for variant in VCF.
    
    Returns:
        list of HGNC symbols from the INFO field
    """
    
    line = line.strip().split("\t")
    
    # only include autosomal genes, so as to not worry about different mutation
    # rates on chrX between males and females.
    chrom = line[1]
    if chrom in ["X", "chrX", "23", "Y", "chrY", "24", "MT"]:
        return []
    
    # ignore variants in CNVs, since those variants aren't accessible to our
    # de novo SNV mutation framework.
    alt = line[4]
    if alt in ["<DEL>", "<DUP>"]:
        return []
    
    info = parse_vcf_info(line[7])
    symbols = ""
    if "HGNC" in info:
        symbols = info["HGNC"]
    elif "HGNC_ALL" in info:
        symbols = info["HGNC_ALL"]
    
    symbols = re.split(',|&|\\|', symbols)
    
    return(symbols)

def get_genes(path):
    """ find the genes (as HGNC symbols) in a VCF
    
    Args:
        path: path to the VCF for a person e.g /path/to/vcf.vcf.gz
    
    Returns:
        set of HGNC symbols
    """
    
    vcf = gzip.open(path, "r")
    if IS_PYTHON3:
        vcf = gzip.open(path, "rt")
    
    header = get_vcf_header(vcf)
    exclude_vcf_header(vcf)
    
    genes = set([])
    for line in vcf:
        genes |= set(get_symbols_from_line(line))
    
    # remove the NA symbols
    genes.discard(".")
    genes.discard("")
    
    return genes

def get_all_hgnc_symbols(n_probands=50, trios=TRIOS, families=FAMILIES):
    """ obtains the HGNC symbols for all SNVs and indels from proband VCFs
    
    Args:
        n_probands: number of proband VCFs to find unique symbols from. we don't
            go through all of the child VCFs, since that could take days
        trios: path to file defining trios.
        families: path to ped file defining individuals, their relationships
            and their VCF paths.
    """
    
    trios = pandas.read_table(trios)
    families = pandas.read_table(families)
    
    probands = families[families["individual_id"].isin(trios["proband_stable_id"])]
    
    all_genes = set([])
    pos = 0
    while pos < n_probands:
        row = probands.iloc[pos]
        person_id = row["individual_id"]
        vcf_path = row["path_to_vcf"]
        
        genes = get_genes(vcf_path)
        all_genes |= genes
        
        pos += 1
    
    return sorted(all_genes)
