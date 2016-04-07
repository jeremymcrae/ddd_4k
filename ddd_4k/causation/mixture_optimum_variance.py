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

import numpy
import pandas
from scipy.stats import gaussian_kde

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from mupit.constants import LOF_CQ, MISSENSE_CQ

from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.goodness_of_fit import get_goodness_of_fit
from ddd_4k.causation.sample_excess import sample_excess

import seaborn
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def variance_around_optimum(de_novos, expected, mono, optimum, bins, iterations=1000):
    '''
    
    Args:
        mono: dictionary of haploinsufficient and nonhaploinsufficent dominant
            DD-associated genes
    '''
    
    merged = merge_observed_and_expected(de_novos, expected)
    missense_excess = aggregate(merged, ["missense"], normalise=True, bins=bins)
    
    # define the variants to sample from
    hi_lof_variants = de_novos[de_novos['hgnc'].isin(mono['haploinsufficient']) &
        de_novos['consequence'].isin(LOF_CQ | MISSENSE_CQ)]
    non_hi_missense_variants = de_novos[de_novos['hgnc'].isin(mono['nonhaploinsufficient']) &
        de_novos['consequence'].isin(MISSENSE_CQ)]
    
    # # count the number of missense across all genes, then calculate how many
    # # of each type to sample.
    excess_target = (missense_excess['observed'] - missense_excess['expected']).sum()
    lof_target = int(excess_target * optimum)
    gof_target = int(excess_target * (1 - optimum))
    
    optimums = []
    for x in range(iterations):
        print(x)
        
        lof_excess = sample_excess(hi_lof_variants, expected, ['lof', 'missense'],  lof_target, mono, bins)
        gof_excess = sample_excess(non_hi_missense_variants, expected, ['missense'], gof_target, mono, bins)
        
        fits = get_goodness_of_fit(missense_excess, lof_excess, gof_excess)
        
        optimal = list(fits["proportion"])[numpy.argmin(fits["goodness_of_fit"])]
        optimums.append(optimal)
    
    plot_uncertainty(optimums)
    
    return numpy.median(optimums)

def plot_uncertainty(optimums, optimal):
    '''
    '''
    
    fig = pyplot.figure(figsize=(6,6))
    ax = fig.gca()
    
    e = ax.hist(optimums, bins=10)
    
    e = ax.axvline(optimal, linestyle='dashed', color='red')
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel('Proportion of missense as loss-of-function')
    e = ax.set_ylabel('Frequency')
    
    fig.savefig('results/proportion_uncertainty.pdf', format='pdf', bbox_inches='tight', pad_inches=0, transparent=True)
    
