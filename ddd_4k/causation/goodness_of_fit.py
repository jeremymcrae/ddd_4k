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

from __future__ import division, print_function

import pandas

def get_goodness_of_fit(missense_excess, lof_excess, gof_excess):
    """ identify optimal mixing proportion of HI and non-HI genes to reproduce
    observed frequencies across the pLI bins.
    
    Args:
        missense_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for all genes with observed candidate de novos.
        lof_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for LoF DNMs in dominant haploinsufficient genes.
        gof_excess: pandas DataFrame of observed to expected differences across
            the pLI bins, for missense DNMs in dominant nonhaploinsufficient genes
    
    Returns:
        proportion of loss-of-function variants required to best capture the
        observed frequencies at pLI bins.
    """
    
    increments = 200.0
    lof_freqs = [ x/increments for x in range(int(increments) + 1) ]
    difference = []
    for lof_frequency in lof_freqs:
        mis_frequency = 1 - lof_frequency
        
        mixed = lof_excess["delta"] * lof_frequency + \
            gof_excess["delta"] * mis_frequency
        
        difference.append(sum((mixed - missense_excess["delta"])**2))
    
    return pandas.DataFrame({"proportion": lof_freqs, "goodness_of_fit": difference})
