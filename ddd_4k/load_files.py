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
import numpy

from mupit.constants import LOF_CQ, MISSENSE_CQ

def open_de_novos(path, validations=None, exclude_invalid=True, exclude_synonymous=True):
    """ load the de novo dataset
    
    Args:
        path: path to known developmental disorder genes data file.
        validations: path to de novo validation results data file. If unused,
            then we don't exclude any variants that failed validation.
        exclude_invalid: whether to remove the sites which failed validation
        exclude_synonymous: whether to remove synonymous sites
    
    Returns:
        DataFrame for the de novo candidates.
    """
    
    de_novos = pandas.read_table(path)
    
    # figure out whether the sites have loss-of-function consequences
    recode = dict(zip(list(LOF_CQ) + list(MISSENSE_CQ), \
        ["loss-of-function"] * len(LOF_CQ) + ["functional"] * len(MISSENSE_CQ)))
    de_novos["category"] = de_novos["consequence"].map(recode)
    
    if exclude_synonymous:
        de_novos = de_novos[de_novos["consequence"].isin(LOF_CQ | MISSENSE_CQ)]
    
    # remove candidates which have been excluded by validation tests
    if validations:
        validations = pandas.read_table(validations)
        de_novos = de_novos.merge(validations, how="left",
            left_on=["person_stable_id", "chrom", "pos", "ref", "alt", \
                "symbol", "consequence"],
            right_on=["person_id", "chrom", "start_pos", "ref_allele", \
                "alt_allele", "hgnc", "consequence"])
    else:
        de_novos["status"] = numpy.nan
    
    if exclude_invalid:
        de_novos = de_novos[~de_novos["status"].isin(["false_positive", "inherited"])]
    
    de_novos["hgnc"] = de_novos["symbol"]
    
    return de_novos

def open_known_genes(path):
    """ open the dataset of known developmental disorder genes
    
    Args:
        path: path to known developmental disorder genes data file.
    
    Returns:
        DataFrame for the known genes.
    """
    
    with open(path) as handle:
        header = handle.readline()
    
    sep = "\t"
    if header.count("|") > 0:
        sep = "|"
    
    genes = pandas.read_table(path, sep=sep, index_col=False)
    if "ddg2p_status" in genes.columns:
        genes = genes[genes["ddg2p_status"] != "Possible DD Gene"]
    else:
        genes = genes[genes["type"] != "Possible DD Gene"]
    
    if "gencode_gene_name" not in genes.columns:
        genes["gencode_gene_name"] = genes["gene"]
    
    return genes

def open_phenotypes(pheno_path, alt_ids_path=None):
    """ load the phenotypes dataset, and merge the alternate IDs
    
    Args:
        pheno_path: path to phenotypes data file.
        alt_ids_path: path to data file containing alternate IDs.
    
    Returns:
        DataFrame for the phenotypes.
    """
    
    pheno = pandas.read_table(pheno_path, na_values=["NA", "Unknown"])
    
    if alt_ids_path is not None:
        # load the alternate IDs for the individuals. Remove the parental IDs, so we
        # can convert the decipher ID to an integer, so we can merge the dataframe
        # with the phenotypes dataframe, which stores the decipher ID column as an
        # integer
        alt_ids = pandas.read_table(alt_ids_path)
        alt_ids = alt_ids[~alt_ids["decipher_id"].str.contains(":")]
        alt_ids["decipher_id"] = alt_ids["decipher_id"].astype(int)
        
        # remove the duplicate rows (due to having multiple sanger IDs),
        # otherwise we get duplicate phenotype rows
        alt_ids = alt_ids.drop_duplicates(["decipher_id", "person_stable_id"])
        
        pheno = pheno.merge(alt_ids, how="left", left_on="patient_id", right_on="decipher_id")
    
    return pheno

def open_families(families_path, trios_path):
    """ load the families dataset, so we know the relationships between people
    
    Args:
        families_path: path to family relationships data file.
        datatypes_path: path to person datatypes data file.
    
    Returns:
        DataFrame for the families.
    """
    
    trios = pandas.read_table(trios_path)
    families = pandas.read_table(families_path)
    
    families = families.merge(trios, "left", left_on="individual_id", right_on="proband_stable_id")
    
    return families

def add_decipher_ids(de_novos, sanger_ids_path):
    """ map the DDD stables IDs to decipher IDs.
    
    open the sanger IDs table, strip out the parental IDs, and map the DDD
    stable IDs to decipher IDs, since we only want to report decipher IDs to
    clinicians.
    
    Args:
        de_novos: pandas dataframe of protein-altering (missense and
            loss-of-function) de novo mutations.
        sanger_ids_path: path to file mapping between DDD and decipher IDs.
    
    Returns:
        pandas dataframe of de novos with an extra column of decpiher IDs.
    """
    
    sanger_ids = pandas.read_table(sanger_ids_path, sep="\t")
    sanger_ids = sanger_ids[~sanger_ids["decipher_id"].str.contains(":")]
    recode = dict(zip(sanger_ids["person_stable_id"], sanger_ids["decipher_id"]))
    
    de_novos["decipher_id"] = de_novos["person_stable_id"].map(recode)
    
    return de_novos

def get_significant_results(results_path, threshold, column="p_min"):
    """ find the resuls for genes at genome-wide significance.
    
    Args:
        results_path: path to results
        threshold: p-value threshold for assigning genome-wide significance.
        column: name of column containing p-values.
    
    Returns:
        pandas Dataframe of results for genes at genome-wide significance.
    """
    
    results = pandas.read_table(results_path, sep="\t")
    
    # find the genes that exceed a multiple testing corrected genonmewide threshold
    table = results[results[column] < threshold]
    
    # ensure that we only select genes where the DDD has at least one de novo,
    # otherwise we will pick up genes where our data does not contribute to the
    # classififcation. Currently this only excludes one gene.
    table = table[~table["ddd.p_func"].isnull()]
    
    return table
