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
from scipy.stats import gaussian_kde

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.causation.aggregate_pli_bins import aggregate
from ddd_4k.causation.goodness_of_fit import get_goodness_of_fit

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})


def permute_fits(merged, hi_merged, non_hi_merged, true_fit, draws=5000):
    ''' create permuted fits for the optimal mixing proportions.
    
    Args:
        merged: pandas DataFrame of genes with observed and expected counts for
            Lof and missense de novo mutations. This is for all genes in the genome.
        hi_merged: pandas DataFrame of genes with observed and expected counts
            for Lof and missense de novo mutations. This is for dominant
            haploinsufficient DD-associated genes.
        non_hi_merged: pandas DataFrame of genes with observed and expected
            counts for Lof and missense de novo mutations. This is for dominant
            nonhaploinsufficient DD-associated genes.
        true_fit: pandas Dataframe of goodness of fit for true data
        draws: number of permutations to run.
    '''
    
    permuted_fits = None
    for x in range(draws):
        print(x)
        
        hi_temp = merged.ix[random.sample(merged.index, len(hi_merged)), ]
        non_hi_temp = merged.ix[random.sample(merged.index, len(non_hi_merged)), ]
        lof_excess = aggregate(hi_temp, ["lof", "missense"], normalise=True)
        gof_excess = aggregate(non_hi_temp, ["missense"], normalise=True)
        
        temp_fits = get_goodness_of_fit(missense_excess, lof_excess, gof_excess)
        
        col_name = 'sample_{0:04d}'.format(x)
        temp_fits = temp_fits.rename(columns={'goodness_of_fit': col_name})
        
        if permuted_fits is None:
            permuted_fits = temp_fits
        else:
            permuted_fits[col_name] = temp_fits[col_name]
    
    plot_permuted_fits(permuted_fits, true_fit, output='permuted_fits.pdf')
    plot_permuted_optimum_densities(fits, output='permuted_densities.pdf')

def plot_permuted_fits(permuted_fits, true_fit, output='permuted_fits.pdf'):
    ''' plot the permuted fits across the proportion range, to show where
    randomly sampled sits appear.
    
    Args:
        permuted_fits: pandas DataFrame of fits, for samples 1 to n (5000ish)
        true_fit: pandas Dataframe of goodness of fit for true data
        output: path to save output pdf to.
    '''
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    # only plot the first 100 fits, otherwise the plot gets excessively overlaid.
    for x in range(100):
        col_name = 'sample_{0:04d}'.format(x)
        rescaled = numpy.log10(permuted_fits[col_name])
        e = ax.plot(permuted_fits['proportion'], rescaled, color='gray', alpha=0.5, marker="None")
    
    # also mark the position of the fit we have for our true data.
    e = ax.plot(true_fit['proportion'], numpy.log10(true_fit['goodness_of_fit']),
        color='red', marker="None")
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("log10(goodness of fit)")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format='pdf', bbox_inches='tight', pad_inches=0,
        transparent=True)

def plot_permuted_optimum_densities(fits, output='permuted_densities.pdf'):
    ''' plot the distribution of optimal mixing proportions from permuted data
    
    Args:
        fits: pandas DataFrame of fits, for samples 1 to n (5000ish)
        output: path to save output pdf to.
    '''
    
    optims = []
    for x in range(draws):
        col_name = 'sample_{0:04d}'.format(x)
        optimal = list(fits["proportion"])[numpy.argmin(fits[col_name])]
        optims.append(optimal)
    
    density = gaussian_kde(optims)
    x = numpy.arange(0.0, 1.0, 0.002)
    y = density(x)
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    e = ax.plot(x, y)
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("Density")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)
