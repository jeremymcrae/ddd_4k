"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

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

def open_de_novos(path):
    """ load the de novo dataset
    
    Args:
        path: path to known developmental disorder genes data file.
    
    Returns:
        DataFrame for the de novo candidates.
    """
    
    de_novos = pandas.read_table(path, na_filter=False)
    
    # define the loss-of-function and "missense" consequences
    lof_cq = ["transcript_ablation", "splice_donor_variant",
        "splice_acceptor_variant", "stop_gained", "frameshift_variant",
        "coding_sequence_variant", "start_lost", "initiator_codon_variant"]
    missense_cq = ["stop_lost", "inframe_insertion", "inframe_deletion",
        "missense_variant", "transcript_amplification", "protein_altering_variant"]
    
    # figure out whether the sites have loss-of-function consequences
    de_novos["category"] = None
    de_novos["category"][de_novos["consequence"].isin(lof_cq)] = "loss-of-function"
    de_novos["category"][de_novos["consequence"].isin(missense_cq)] = "functional"
    
    de_novos = de_novos[de_novos["consequence"].isin(lof_cq + missense_cq)]
    
    # recode the sex codes to full names
    de_novos.sex[de_novos.sex == "M"] = "male"
    de_novos.sex[de_novos.sex == "F"] = "female"
    
    return de_novos

def open_known_genes(path):
    """ open the dataset of known developmental disorder genes
    
    Args:
        path: path to known developmental disorder genes data file.
    
    Returns:
        DataFrame for the known genes.
    """
    
    genes = pandas.read_table(path, sep="|", na_filter=False)
    genes = genes[genes["ddg2p_status"] != "Possible DD Gene"]
    
    return genes

def open_phenotypes(pheno_path, alt_ids_path=None):
    """ load the phenotypes dataset, and merge the alternate IDs
    
    Args:
        pheno_path: path to phenotypes data file.
        alt_ids_path: path to data file containing alternate IDs.
    
    Returns:
        DataFrame for the phenotypes.
    """
    
    pheno = pandas.read_table(pheno_path, na_filter=False)
    
    if alt_ids_path is not None:
        # load the alternate IDs for the individuals. Remove the parental IDs, so we
        # can convert the decipher ID to an integer, so we can merge the dataframe
        # with the phenotypes dataframe, which stores the decipher ID column as an
        # integer
        alt_ids = pandas.read_table(alt_ids_path, na_filter=False)
        alt_ids = alt_ids[~alt_ids["decipher_id"].str.contains(":")]
        alt_ids["decipher_id"] = alt_ids["decipher_id"].astype(int)
        
        # remove the duplicate rows (due to having multiple sanger IDs),
        # otherwise we get duplicate phenotype rows
        alt_ids = alt_ids.drop_duplicates(["decipher_id", "person_stable_id"])
        
        pheno = pheno.merge(alt_ids, how="left", left_on="patient_id", right_on="decipher_id")
    
    return pheno

def open_families(families_path, datatypes_path):
    """ load the families dataset, so we know the relationships between people
    
    Args:
        families_path: path to family relationships data file.
        datatypes_path: path to person datatypes data file.
    
    Returns:
        DataFrame for the families.
    """
    
    datatypes = pandas.read_table(datatypes_path, na_filter=None)
    families = pandas.read_table(families_path, na_filter=None)
    
    families = families.merge(datatypes, "left", left_on="individual_id", right_on="person_stable_id")
    
    return families
    
    
