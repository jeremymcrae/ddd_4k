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
