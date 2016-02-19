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

import itertools

import pandas
import numpy

def aggregate(merged, consequences=None, normalise=True):
    """ get the difference between observed and expected counts across pLI bins.
    
    This returns the difference in the number of observed de novo mutations to
    expected numbers of de novo mutations. The difference is given as a delta
    and as a ratio.
    
    Args:
        merged: pandas DataFrame containing one row per gene, which includes
            columns for observed and expected counts, by lof and missense
            categories. This also includes a column indicating which pLI
            vigicile each genes falls into.
        consequences: list of consequence types to include, such as ["missense"],
            or ["missense", "lof"] or ["lof"]
        normalise: boolean for whether to normalise the differences so that the
            total sums to one.
    
    Returns:
        pandas DataFrame with columns for pLI bin, and observed to expected
        differences as delta and ratio.
    """
    
    if consequences is None:
        consequences = ["lof", "missense"]
        
    observed_columns = [ x + "_observed" for x in consequences ]
    expected_columns = [ x + "_expected" for x in consequences ]
    
    aggregated = merged.pivot_table(rows="pLI_bin",
        values=observed_columns + expected_columns, aggfunc=sum)
    
    aggregated["pLI_bin"] = aggregated.index
    
    aggregated["observed"] = aggregated[observed_columns].sum(axis=1)
    aggregated["expected"] = aggregated[expected_columns].sum(axis=1)
    
    aggregated["delta"] = aggregated["observed"] - aggregated["expected"]
    aggregated["ratio"] = aggregated["observed"]/aggregated["expected"]
    
    if normalise:
        aggregated["delta"] = aggregated["delta"]/sum(aggregated["delta"])
        aggregated["ratio"] = aggregated["ratio"]/sum(aggregated["ratio"])
    
    # figure out the number of bins
    combos = itertools.permutations(aggregated["pLI_bin"], 2)
    n_bins = 1/min([ abs(x[0] - x[1]) for x in combos ])
    
    # make sure we have all the quantile bins, even for tables that might lack
    # genes in a given pLI bin
    for x in [ x/float(n_bins) for x in range(int(n_bins)) ]:
        if not any([ almost_equal(x, y) for y in aggregated["pLI_bin"] ]):
            row = pandas.DataFrame({"pLI_bin": [x], "delta": [0], "ratio": [0]})
            aggregated = aggregated.append(row, ignore_index=True)
    
    aggregated = aggregated[["pLI_bin", "observed", "expected", "delta",
        "ratio"]].copy()
    aggregated = aggregated.sort("pLI_bin")
    aggregated.index = range(len(aggregated))
    
    return aggregated

def almost_equal(a, b, tol=0.0001):
    """ check whether two numbers are roughly equivalent.
    
    We have constructed a set of pLI bins, but some tables lack certain bins if
    there weren't any genes within a given pLI bin. In order to insert the
    missing bins, we need to check whether each possible bin already occurs in
    the table. First we have figured out the bin size (by getting the smallest
    increment between the known bins), then figured out how many bins that would
    allow. Unfortunately, this can give slightly different bin values than when
    the bins were constructed. We get around this by checking if there are bin
    values roughly close.
    
    Args:
        a: first number, for example 0.05
        b: second number, for example 0.0500001
        tol: percentage difference tolerated.
    
    Returns:
        true/false for equivalence.
    """
    
    try:
        return a==b or (max(a, b)/min(a, b)) - 1 < tol
    except ZeroDivisionError:
        return False
