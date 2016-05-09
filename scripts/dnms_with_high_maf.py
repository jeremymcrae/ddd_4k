'''
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
'''

from __future__ import print_function

import sys
import os
import argparse

import pandas
import pysam

from denovonear.load_gene import get_transcript_ids, construct_gene_object
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import get_codon_info, get_gene_range, \
    get_mutated_aa

def get_options():
    ''' get the command line options
    '''
    
    parser = argparse.ArgumentParser(description='determines relative mutability '
        'of synonymous variants in exome with MAF < 0.01 vs exome with maf >= 0.01')
    parser.add_argument('--gencode',
        default='/lustre/scratch113/projects/ddd/users/jm33/all_gencode_genes_v19.txt.gz',
        help='path to gencode genes table')
    parser.add_argument('--rates',
        help='path to mutation rates at trinucleotide sequence contexts.',
        default='/nfs/users/nfs_j/jm33/apps/denovonear/data/forSanger_1KG_mutation_rate_table.txt')
    
    args = parser.parse_args()
    
    return args

def prepare_tabix():
    ''' get a dictionary of tabix handles, indexed by dataset key.
    '''
    
    resources = '/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2016-03-22/'
    
    paths = {
        'ddd_internal': {'path': os.path.join(resources,
            'ANNOTATE_TABIX/MAF/ddd-unaffected-parental-maf-2014-12-12.tsv.gz')},
        'esp': {'path': os.path.join(resources,
            'ANNOTATE_TABIX/MAF/ESP6500SI_21012013_MAF.bed.gz')},
        'uk10k': {'path': os.path.join(resources,
            'ANNOTATE_TABIX/MAF/UK10K_COHORT_20130116_AF.bed.gz')},
        'exac': {'path': os.path.join(resources,
            'ExAC_MAF/ExAC.r0.3.sites.vep.vcf.gz')},
        '1000g': {'path': os.path.join(resources,
            '1000G_MAF/ALL.chrXXX.vcf.gz')}
        }
    
    for key in paths:
        path = paths[key]['path']
        if key == '1000g':
            paths[key]['tabix'] = {}
            for chrom in range(1, 23) + ['X', 'Y']:
                chrom = str(chrom)
                paths[key]['tabix'][chrom] = pysam.TabixFile(path.replace('XXX', chrom))
        else:
            paths[key]['tabix'] = pysam.TabixFile(path)
    
    return paths

def parse_info(info):
    """ parse the info field from a VCF into a python dictionary
    
    Args:
        info: text string for an info field.
    
    Returns:
        dictionary of info values, indexed by the info keys.
    """
    
    info = info.strip().split(";")
    
    info_dict = {}
    for item in info:
        if "=" in item:
            item = item.split("=")
            key = item[0]
            value = item[1]
        else:
            key = item
            value = True
        
        info_dict[key] = value
    
    return info_dict

def fetch_rows(tbx, chrom, pos, key):
    ''' pull the correct rows out of the VCF, and handle errors appropriately
    
    Args:
        tbx: pysam TabixFile object for requesting rows.
        chrom: chromosome to be checked.
        pos: nucleotide position.
        key: text key for the population being checked, for reporting errors if
            anything goes wrong
    
    Returns:
        list of rows covering the required genome position.
    '''
    
    try:
        rows = list(tbx.fetch(str(chrom), pos-1, pos + 1))
    except ValueError as error:
        # the UK10K dataset lacks variants on the sex chromosomes. Just
        # ignore files where a chromosome is missing.
        if chrom not in tbx.contigs:
            rows = []
        else:
            raise error('unable to tabix access {0}:{1} in {2}'.format(chrom, pos, key))
    
    return rows

def parse_maf(row, key, alt):
    ''' get the minor allele frequency from a population row
    
    Args:
        row: list of text strings from a tab-separated line in a tabixed file
        key: dataset ID, which determines how we parse the data in the row
        alt: base (e.g. 'A') of the alternate allele.
    
    Returns:
        minor allele frequency extracted from the row. Where multiple MAFs are
        available, use the highest MAF value.
    '''
    
    if key in ['ddd_internal', 'esp', 'uk10k']:
        maf = float(row[5])
    else:
        idx = row[4].split(',').index(alt)
        info = parse_info(row[7])
        mafs = []
        if key == '1000g':
            for pop in ['SAS_AF', 'EUR_AF', 'AFR_AF', 'AMR_AF', 'EAS_AF']:
                mafs.append(float(info[pop].split(',')[idx]))
        elif key == 'exac':
            for pop in ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS']:
                ac = info['AC_{}'.format(pop)].split(',')[idx]
                an = info['AN_{}'.format(pop)]
                try:
                    mafs.append(float(ac)/float(an))
                except ZeroDivisionError:
                    mafs.append(0.0)
        
        maf = max(mafs)
    
    return maf

def lookup_maf(chrom, pos, alt):
    ''' get the maimum MAF from a variety of populations
    
    Args:
        chrom: chromosome for the variant.
        pos: nucleotide position of the variant.
        alt: alternate allele of the variant.
    
    Returns:
        maximum minor allele frequency as float.
    '''
    
    max_af = 0.0
    
    for key in MAFS:
        tbx = MAFS[key]['tabix']
        if key == '1000g':
            tbx = MAFS[key]['tabix'][chrom]
        
        rows = fetch_rows(tbx, chrom, pos, key)
        
        # make sure the rows are for the correct position and alt allele
        rows = [ x for x in rows if x.split('\t', 2)[1] == str(pos) ]
        rows = [ x for x in rows if alt in x.split('\t')[4].split(',') ]
        
        if len(rows) == 0:
            continue
        
        if len(rows) != 1:
            print('too many rows at {0}:{1} in {2} - {3}'.format(chrom, pos, key, rows))
        
        row = rows[0].split('\t')
        maf = parse_maf(row, key, alt)
        
        if maf > max_af:
            max_af = maf
    
    return max_af

def load_transcript(ensembl, hgnc):
    ''' load a Transcript object for a gene
    
    Args:
        ensembl: EsemblRequest object
        hgnc: HGNC symbol for a gene.
    
    Returns:
        Transcript object for gene (or None if not suitable gene found)
    '''
    
    transcript_ids = get_transcript_ids(ensembl, hgnc)
    
    # construct a Transcript object for a transcript ID, run though the
    # transcript IDs by descending transcript lengths until we get a suitable
    # transcript.
    for key in sorted(transcript_ids, key=transcript_ids.get, reverse=True):
        try:
            return construct_gene_object(ensembl, key)
        except ValueError:
            pass
    
    return None

class NonsynonymousError(Exception):
    def __str__(self):
        return 'not a synonymous site'

def get_maf(transcript, bp, codon, seq, alt):
    ''' get the minor allele frequency for a synonymous variant
    
    Args:
        transcript: Trnascript object for a gene
        bp: position of variant on chromosome
        codon: dictionary of information about the codon.
        seq: trinucleotide genome sequence centered on the variant position
        alt: alternate single base (e.g. 'G') for checking variant consequence,
            and minor allele frequency.
    
    Returns:
        minor allele frequency as float. Raises error if site isn't synonymous.
    '''
    
    mutated_aa = get_mutated_aa(transcript, alt, codon['codon_seq'],
        codon['intra_codon'])
    
    # only use synonymous sites
    if codon["initial_aa"] != mutated_aa:
        raise NonsynonymousError
    
    # get the ref and alt alleles, oriented to the '+' strand
    ref = seq[1]
    if transcript.get_strand() == '-':
        ref = transcript.reverse_complement(ref)
        alt = transcript.reverse_complement(alt)
    
    return lookup_maf(transcript.get_chrom(), bp, alt)

def exome_at_high_maf(transcript, mut_dict):
    ''' find the mutability adjusted proportion of the exome at high MAF sites.
    
    Args:
        transcript: Transcript object for a gene transcript
        mut_dict: dictionary of trinucleotide sequence based mutation rates.
    
    Returns:
        tuple of cumulative mutation rates at synonymous sites. The first entry
        is for sites with low MAF, the second entry is for sites with high MAF.
    '''
    
    bases = set(['A', 'C', 'G', 'T'])
    
    low_maf_rate = 0
    high_maf_rate = 0
    
    for bp in range(*get_gene_range(transcript)):
        
        # ignore noncoding positions, since we only want synonymous sites
        if not transcript.in_coding_region(bp):
            continue
        
        codon = get_codon_info(transcript, bp, boundary_dist=0)
        seq = transcript.get_trinucleotide(bp)
        
        if transcript.get_strand() == '-':
            seq = transcript.reverse_complement(seq)
        
        # drop the initial base, since we want to mutate to other bases
        for alt in bases - set(seq[1]):
            
            try:
                maf = get_maf(transcript, bp, codon, seq, alt)
            except NonsynonymousError:
                continue
            
            alt_seq = seq[0] + alt + seq[2]
            if maf >= 0.01:
                high_maf_rate += mut_dict[seq][alt_seq]
            else:
                low_maf_rate += mut_dict[seq][alt_seq]
    
    return low_maf_rate, high_maf_rate

def main():
    '''
    '''
    
    args = get_options()
    
    global MAFS
    MAFS = prepare_tabix()
    
    genes = pandas.read_table(args.gencode, compression='gzip')
    genes = genes[genes['gene_type'].isin(['protein_coding', 'polymorphic_pseudogene'])]
    ensembl = EnsemblRequest('cache', 'grch37')
    
    mut_dict = load_mutation_rates(args.rates)
    
    low_rate = 0
    high_rate = 0
    
    for hgnc in sorted(set(genes.gene)):
        
        transcript = load_transcript(ensembl, hgnc)
        
        if transcript is None:
            continue
        
        lo, hi = exome_at_high_maf(transcript, mut_dict)
        low_rate += lo
        high_rate += hi
        
        print('{0}\t{1}\t{2}'.format(hgnc, lo, hi))
    
    print('totals - low maf sites: {0}, high maf sites: {1}'.format(low_rate, high_rate))
    print('proportion in low maf sites: {}'.format(high_rate/low_rate))

if __name__ == '__main__':
    main()
