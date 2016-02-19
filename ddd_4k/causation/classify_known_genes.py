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

def classify_monoallelic_genes(known):
    """ classify monoallelic genes into haploinsufficient and
        nonhaploinsufficient sets
    
    Args:
        known: pandas DataFrame of known developmental dfisorder genes
    
    Returns:
        dictionary of gene sets, one entry for haploinsufficient genes, one
        entry for nonhaploinsufficient genes
    """
    
    # identify known haploinsufficient dominant genes
    monoallelic = known[known["mode"].isin(["Monoallelic", "X-linked dominant"])]
    haploinsufficient = set(monoallelic["gencode_gene_name"][
        monoallelic["mech"].isin(["Loss of function"])])
    
    # identify known nonhaploinsufficient dominant genes
    nonhaploinsufficient = set(monoallelic["gencode_gene_name"][
        monoallelic["mech"].isin(["Activating", "Dominant negative"])])
    
    # remove genes which fall into both categories, since de novos in those will
    # be a mixture of the two models, and we need to cleanly separate the
    # underlying models.
    overlap = haploinsufficient & nonhaploinsufficient
    haploinsufficient -= overlap
    nonhaploinsufficient -= overlap
    
    return {"haploinsufficient": haploinsufficient,
        "nonhaploinsufficient": nonhaploinsufficient}
