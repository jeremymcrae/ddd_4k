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

from __future__ import division

import pandas

def include_constraints(table, constraints):
    """ add a column for pLI constraints scores
    
    Args:
        table: pandas DataFrame, with one gene per row
        constraints: table of constraints scores, including columns for HGNC
            symbols ("gene") and constraints score ("pLI").
    
    Returns:
        same pandas DataFrame, with an extra pLI column (and minus genes without
        pLI scores.)
    """
    
    # load the pLI data, and append column to the table
    recode = dict(zip(constraints["gene"], constraints["pLI"]))
    table["pLI"] = table["hgnc"].map(recode)
    
    # drop genes which lacka constraint score
    table = table[~table["pLI"].isnull()]
    
    return table

def get_constraint_bins(table, bins=20, rate_correct=False):
    """ identify the quantile bins each gene falls into.
    
    Args:
        table: pandas DataFrame, with one gene per row.
        bins: number of quantiles to divide into, or list of quantiles.
        rate_correct: whether to adjust the bins to account for different
            mutation rates of genes depending on the pLI score. If used, each
            bin will expect equal numbers of synonymous mutations.
    
    Returns:
        pandas Series, indicating the bin each gene falls into.
    """
    
    if type(bins) == list:
        quantiles = bins
    else:
        # identify which pLI quantile each gene falls into
        quantiles = [ x/float(bins) for x in range(bins + 1) ]
    
    if not rate_correct:
        pLI_bin, bins = pandas.qcut(table["pLI"], q=quantiles,
            labels=quantiles[:-1], retbins=True)
    else:
        summed_synonymous = sum(table["synonymous_expected"])
        
        pLI_bin = []
        pos = 0
        current_sum = 0
        pli_sorted = table.sort("pLI")
        for key, row in pli_sorted.iterrows():
            if current_sum/summed_synonymous > quantiles[pos + 1]:
                pos += 1
            
            current_sum += row["synonymous_expected"]
            pLI_bin.append(quantiles[pos])
        
        pli_sorted["pLI_bin"] = pLI_bin
        pli_sorted = pli_sorted.sort_index()
        pLI_bin = pli_sorted["pLI_bin"]
    
    return pLI_bin
