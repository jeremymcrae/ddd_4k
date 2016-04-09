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

import random

import pandas

from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.aggregate_pli_bins import aggregate

def sample_excess(variants, expected, consequences, target, mono, bins, get_variants=False):
    ''' sample from variant lists to get a target excess
    
    Args:
        variants: pandas DataFrame of variants to sample from
        expected: numbers of expected mutations per gene
        consequences: list of consequence types to include
        target: desired number of excess mutations
        mono: dictionary of haploinsufficient and nonhaploinsufficent dominant
            DD-associated genes
        bins: list of pLI quantiles.
    
    Returns:
        pandas DataFrame of excess by pLI bin.
    '''
    
    lo = 0
    hi = int(target * 2)
    delta = None
    
    if target == 0:
        if get_variants:
            return pandas.DataFrame(columns=variants.columns)
        else:
            bins = [ x for x in bins if x != 1.0 ]
            return pandas.DataFrame({'pLI_bin': bins, 'observed': [0] * len(bins),
                'expected': [0] * len(bins), 'delta': [0] * len(bins)})
    
    while delta is None or abs(delta) > 5:
        
        mid = int((lo + hi)/2.0)
        
        if mid == 0:
            delta = target
        else:
            idx = [ random.choice(variants.index) for x in range(mid) ]
            merged = merge_observed_and_expected(variants.ix[idx, ], expected)
            
            if consequences == ['missense']:
                merged = merged[merged["hgnc"].isin(mono["nonhaploinsufficient"])]
            else:
                merged = merged[merged["hgnc"].isin(mono["haploinsufficient"])]
            
            excess = aggregate(merged, consequences, normalise=True, bins=bins)
            delta = (excess['observed'] - excess['expected']).sum() - target
        
        if delta > 0:
            hi = mid
            if hi - lo < 2:
                lo -= 4
        else:
            lo = mid
            if hi - lo < 2:
                hi += 4
    
    if get_variants:
        return variants.ix[idx, ]
    else:
        return excess
