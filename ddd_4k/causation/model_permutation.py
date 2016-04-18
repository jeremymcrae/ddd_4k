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
from scipy.stats import linregress, gaussian_kde, norm

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.causation.sample_excess import sample_excess
from ddd_4k.causation.goodness_of_fit import get_goodness_of_fit

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})


def permute_fits(de_novos, expected, mono, bins, missense_excess, lof_excess, gof_excess, prior, increments=100, permutations=20):
    ''' create permuted fits for the optimal mixing proportions.
    
    Args:
        merged: pandas DataFrame of genes with observed and expected counts for
            Lof and missense de novo mutations. This is for all genes in the genome.
        increments: number of increments to sample across the proprotion range.
    '''
    
    hi_variants = de_novos[de_novos['hgnc'].isin(mono['haploinsufficient'])]
    non_hi_variants = de_novos[de_novos['hgnc'].isin(mono['nonhaploinsufficient'])]
    
    data = pandas.DataFrame({'proportion': [], 'optimum': []})
    excess_target = (missense_excess['observed'] - missense_excess['expected']).sum()
    proportions = [ x/float(increments) for x in range(increments + 1) ]
    fits = pandas.DataFrame({'proportion': [], 'optimal': [], 'goodness_of_fit': []})
    for freq in proportions:
        lof_target = int(excess_target * freq)
        gof_target = int(excess_target * (1 - freq))
        
        print('checking frequency: {}'.format(freq))
        
        if freq == 0 or freq == 1.0:
            continue
        
        optimums = []
        for x in range(permutations):
            
            lof = sample_excess(hi_variants, expected, ['lof', 'missense'], lof_target, mono, bins)
            gof = sample_excess(non_hi_variants, expected, ['missense'], gof_target, mono, bins)
            
            mis = lof + gof
            mis['delta'] = mis['observed'] - mis['expected']
            mis['delta'] = mis['delta']/sum(mis['delta'])
            
            temp = get_goodness_of_fit(mis, lof_excess, gof_excess)
            
            optimal = list(temp["proportion"])[numpy.argmin(temp["goodness_of_fit"])]
            optimums.append(optimal)
        
        temp =  pandas.DataFrame({'proportion': [freq] * len(optimums), 'optimum': optimums})
        data = data.append(temp, ignore_index=True)
        
        optimum = sum(optimums)/float(len(optimums))
        fits = fits.append({'proportion': freq, 'optimal': optimum}, ignore_index=True)
    
    plot_permuted_fits(fits, output='set_proportion_vs_estimated_proportion.pdf')
    
    slope, intercept, r_value, p_value, std_err = linregress(fits['proportion'], fits['optimal'])
    print('slope: {0}, intercept: {1}'.format(slope, intercept))
    get_95_interval(data, slope, intercept, prior)
    
    return (slope, intercept)

def get_95_interval(data, slope, intercept, prior):
    '''
    '''
    
    data.to_csv('proportion_optimums.txt', sep='\t', index=False)
    
    data['adjusted'] = (data['optimum'] - intercept)/slope
    
    likelihoods = pandas.DataFrame({'proportion': [], 'density': []})
    for freq, group in data.groupby('proportion'):
        density = gaussian_kde(group['adjusted'])
        likelihoods = likelihoods.append({'proportion': freq, 'density': density(prior)}, ignore_index=True)
    
    likelihoods = likelihoods.sort('proportion')
    
    likelihoods.to_csv('proportion_likelihoods.txt', sep='\t', index=False)
    
    # scale the likelihoods so the points sum to 1
    likelihoods['scaled_density'] = likelihoods['density']/sum(likelihoods['density'])
    
    likelihoods['cumsum'] = likelihoods['scaled_density'].cumsum()
    
    below = likelihoods[likelihoods['cumsum'] < 0.025]
    above = likelihoods[likelihoods['cumsum'] > 0.975]
    
    lo = max(below['proportion'])
    hi = min(above['proportion'])
    
    print('95% CI around optimal proportion: {0}-{1}'.format(lo, hi))
    
    samples = []
    for i, row in likelihoods.iterrows():
        samples += [row['proportion']] * int(row['density'] * 1000)
    
    # determine a curve for the a normal distribution around the distribution
    # parameters
    mu, std = norm.fit(samples)
    x_values = [ x/float(200) for x in range(200) ]
    y_values = norm.pdf(x_values, mu, std)
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    # also mark the position of the fit we have for our true data.
    e = ax.plot(likelihoods['proportion'], likelihoods['density'], marker='.', linestyle='none')
    e = ax.plot(x_values, y_values, color='gray')
    
    e = ax.set_xlim(0.1, 0.65)
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("Density")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig('permutation.variance.pdf', format='pdf', bbox_inches='tight', pad_inches=0,
        transparent=True)

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
    e = ax.plot(permuted_fits['proportion'], permuted_fits['optimal'], marker='.', linestyle='none')
    
    e = ax.set_xlabel("proportion HI")
    e = ax.set_ylabel("estimated proportion HI")
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    fig.savefig(output, format='pdf', bbox_inches='tight', pad_inches=0,
        transparent=True)
