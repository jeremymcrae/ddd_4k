"""
Copyright (c) 2016 Genome Research Ltd.

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

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.goodness_of_fit import get_goodness_of_fit

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})


def permute_fits(merged, increments=200):
    ''' create permuted fits for the optimal mixing proportions.
    
    Args:
        merged: pandas DataFrame of genes with observed and expected counts for
            Lof and missense de novo mutations. This is for all genes in the genome.
        increments: number of increments to sample across the proprotion range.
    '''
    
    hi_genes = merged[merged['hgnc'].isin(mono['haploinsufficient']) & ((merged['lof_observed'] > 0) | (merged['missense_observed'] > 0))]
    non_hi_genes = merged[merged['hgnc'].isin(mono['nonhaploinsufficient']) & ((merged['lof_observed'] > 0) | (merged['missense_observed'] > 0))]
    
    excess_target = (missense_excess['observed'] - missense_excess['expected']).sum()
    increments = 200
    proportions = [ x/float(increments) for x in range(increments + 1) ]
    fits = pandas.DataFrame({'proportion': [], 'optimal': []})
    for freq in proportions:
        lof_target = int(excess_target * freq)
        gof_target = int(excess_target * (1 - freq))
        
        print(freq)
        
        if freq == 0 or freq == 1.0:
            continue
        
        values = []
        for x in range(20):
            
            lof = [ random.choice(hi_genes.index) for x in range(lof_target) ]
            gof = [ random.choice(non_hi_genes.index) for x in range(gof_target) ]
            
            lof = hi_genes['pLI_bin'].ix[lof].value_counts()
            gof = non_hi_genes['pLI_bin'].ix[gof].value_counts()
            
            missing = set(bins[:-1]) - set(lof.index)
            missing = pandas.Series([0] * len(missing), index=missing)
            lof = lof.append(missing)
            
            missing = set(bins[:-1]) - set(gof.index)
            missing = pandas.Series([0] * len(missing), index=missing)
            gof = gof.append(missing)
            
            lof_excess = pandas.DataFrame({'delta': lof, 'pLI_bin': lof.index}).sort('pLI_bin').reset_index()
            gof_excess = pandas.DataFrame({'delta': gof, 'pLI_bin': gof.index}).sort('pLI_bin').reset_index()
            
            lof_excess['delta'] = lof_excess['delta']/sum(lof_excess['delta'])
            gof_excess['delta'] = gof_excess['delta']/sum(gof_excess['delta'])
            
            temp = get_goodness_of_fit(missense_excess, lof_excess, gof_excess)
            value = list(temp["proportion"])[numpy.argmin(temp["goodness_of_fit"])]
            
            values.append(value)
        
        estimate = sum(values)/float(len(values))
        fits = fits.append({'proportion': freq, 'optimal': estimate}, ignore_index=True)
    
    plot_permuted_fits(fits, output='set_proportion_vs_estimated_ptoportion.pdf')

def plot_permuted_fits(permuted_fits, output='permuted_fits.pdf'):
    ''' plot the permuted fits across the proportion range, to show where
    randomly sampled sits appear.
    
    Args:
        permuted_fits: pandas DataFrame
        output: path to save output pdf to.
    '''
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    # also mark the position of the fit we have for our true data.
    e = ax.plot(permuted_fits['proportion'], permuted_fits['optimal']), marker='.', linestyle='none')
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("estimated proportion HI")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format='pdf', bbox_inches='tight', pad_inches=0,
        transparent=True)
