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
import sys
import gzip

import pandas

from ddd_4k.constants import FAMILIES, TRIOS
from ddd_4k.parse_vcf import parse_vcf_format, parse_vcf_info, get_vcf_header, \
    exclude_vcf_header

IS_PYTHON3 = sys.version_info[0] == 3

cnv_alts = set(["<DEL>", "<DUP>"])
required = ["person_id", "chrom", "pos", "ref", "alt", "ACGH_RC_FREQ50",
    "CALLP", "CALLSOURCE", "cifer", "CNS", "COMMONBACKWARDS", "CONVEXSCORE",
    "CQ", "HGNC", "INTERNALFREQ", "MADL2R", "MEANLR2",
    "NUMBERPROBESCNSOLIDATE", "NUMBERPROBESCONVEX", "RAREBACKWARDS",
    "RAREFORWARDS", "RC50INTERNALFREQ", "SVLEN", "vicar", "WSCORE"]

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to find candidate"
        "de novo CNVs in proband VCFs")
    parser.add_argument("--families", default=FAMILIES, \
        help="Path to family relationships, including sample VCF paths.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Path to table listing all trios.")
    parser.add_argument("--output", default="de_novos_cnvs.txt", \
        help="Path to send output to.")
    
    args = parser.parse_args()
    
    return args

def parse_line(line):
    """ parse a VCF line and return variant details if is a CNV, or None
    
    Args:
        line: line of VCF for a single variant.
    
    Returns:
        CNV details as dictionary, or None if VCF line isn't for a CNV.
    """
    
    # exclude variants that are not CNVs
    alt = line.split("\t", 5)[4]
    if alt not in cnv_alts:
        return None
    
    line = line.strip().split("\t")
    
    # only include CNVs which might be de novo mutations
    format_dict = parse_vcf_format(line[8], line[9])
    cifer = format_dict["CIFER_INHERITANCE"]
    vicar = format_dict["INHERITANCE"]
    
    if cifer != "not_inherited" and vicar != "deNovo":
        return None
    
    info = parse_vcf_info(line[7])
    
    cnv = {}
    for x in required:
        cnv[x] = "NA"
        if x in info:
            cnv[x] = info[x]
    
    cnv["chrom"] = line[0]
    cnv["pos"] = line[1]
    cnv["ref"] = line[3]
    cnv["alt"] = alt
    cnv["cifer"] = cifer
    cnv["vicar"] = vicar
    
    return cnv

def get_de_novo_cnvs(path, person_id):
    """ find the de novo CNVs in a VCF
    
    Args:
        path: path to the VCF for a person e.g /path/to/vcf.vcf.gz
        person_id: ID for the person e.g. DDDP100001
        required: list of columns required in the output
    
    Returns:
        pandas dataframe for the CNVs
    """
    
    vcf = gzip.open(path, "r")
    if IS_PYTHON3:
        vcf = gzip.open(path, "rt")
    
    header = get_vcf_header(vcf)
    exclude_vcf_header(vcf)
    
    cnv_alts = set(["<DEL>", "<DUP>"])
    
    cnvs = pandas.DataFrame(columns=required)
    
    for line in vcf:
        cnv = parse_line(line)
        if cnv is not None:
            cnvs = cnvs.append(cnv, ignore_index=True)
    
    cnvs["person_id"] = person_id
    
    return cnvs

def main():
    args = get_options()
    
    trios = pandas.read_table(args.trios)
    families = pandas.read_table(args.families)
    
    probands = families[families["individual_id"].isin(trios["proband_stable_id"])]
    
    cnv_alts = set(["<DEL>", "<DUP>"])
    
    cnvs = pandas.DataFrame(columns=required)
    for pos in range(len(probands)):
        row = probands.iloc[pos]
        person_id = row["individual_id"]
        vcf_path = row["path_to_vcf"]
        print(person_id, vcf_path)
        
        person_cnvs = get_de_novo_cnvs(vcf_path, person_id)
        cnvs = cnvs.append(person_cnvs, ignore_index=True)
    
    cnvs.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
