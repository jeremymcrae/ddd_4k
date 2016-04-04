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

import random

import numpy
import pandas

from mupit.constants import LOF_CQ, MISSENSE_CQ

from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.model_mixtures import get_goodness_of_fit

def variance_around_optimum(de_novos, expected, mono, optimum, iterations=1000):
    
    merged = merge_observed_and_expected(de_novos, expected)
    missense_excess = aggregate(merged, ["missense"], normalise=True)
    
    # define the variants to sample from
    hi_lof_variants = de_novos[de_novos['hgnc'].isin(mono['haploinsufficient']) &
        de_novos['consequence'].isin(LOF_CQ | MISSENSE_CQ)]
    nonhi_missense_variants = de_novos[de_novos['hgnc'].isin(mono['nonhaploinsufficient']) &
        de_novos['consequence'].isin(MISSENSE_CQ)]
    
    # count the number of missense across all genes, then calculate how many
    # of each type to sample.
    excess_target = (missense_excess['observed'] - missense_excess['expected']).sum()
    hi_excess_target = int(excess_target * optimum)
    non_hi_excess_target = int(excess_target * (1 - optimum))
    
    optimums = []
    for x in range(iterations):
        print(x)
        
        lof_excess = get_excess(hi_lof_variants, expected, ['lof', 'missense'],  hi_excess_target, mono)
        gof_excess = get_excess(non_hi_lof_variants, expected, ['missense'], non_hi_excess_target, mono)
        
        fits = get_goodness_of_fit(missense_excess, lof_excess, gof_excess)
        
        optimal = list(fits["proportion"])[numpy.argmin(fits["goodness_of_fit"])]
        optimums.append(optimal)

def get_excess(variants, expected, consequences, target, mono):
    '''
    '''
    
    idx = [ random.choice(variants.index) for x in range(target) ]
    sampled = variants.ix[idx, ]
    
    merged = merge_observed_and_expected(sampled, expected)
    
    excess = aggregate(merged, consequences, normalise=True)
    
    current = (excess['observed'] - excess['expected']).sum()
    
    if abs(current - target) < 5:
        return excess
    elif current - target < 0:
        # up the number of sampled variants
        diff = (target - current)/2.0
        idx = [ random.choice(variants.index) for x in range(diff) ]
        temp = sampled = variants.ix[idx, ]
        
        return 'low'
    else:
        # lower thenumber of sampled variants
        return 'high'
    
