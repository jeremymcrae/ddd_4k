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
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes
from ddd_4k.constants import DENOVO_PATH, VALIDATIONS

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

PREVIOUS_COHORT = '/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv'

def plot_excess_by_pp_dnm_threshold(de_novos, expected, increments=100):
    ''' show the impact of changing the pp_dnm threshold on the excess DNMs
    
    Args:
        de_novos: pandas DataFrame of candidate de novo variants.
        expected: pandas DataFrame of expected mutation rates for genes.
        increments: number of increments to divide the pp_dnm range into
    '''
    
    pp_dnms = [ x/float(increments) for x in range(increments) ]

    burdens = []
    for threshold in pp_dnms:
        variants = de_novos[(~de_novos["pp_dnm"].isnull() & (de_novos["pp_dnm"] > threshold)) ]
        excess = get_consequence_excess(expected, variants)
        
        burden = excess['missense']['excess'] + excess['loss-of-function']['excess']
        burdens.append(burden)
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(pp_dnms, burdens, marker='.', markersize=10)
    
    e = ax.set_ylim(0, max(burdens) * 1.05)
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel('PP_DNM threshold')
    e = ax.set_ylabel('Excess de novo mutations')
    
    fig.savefig('excess_by_pp_dnm_threshold.pdf', format='pdf',
        bbox_inches='tight', pad_inches=0, transparent=True)

def get_pp_dnms(increments):
    """ make a list of pp_dnms to check ROC curves at.
    
    We need to span the pp_dnm range, but increase sampling at both ends of the
    pp_dnm range, since most variants occur at either end of the pp_dnm range.
    I could set this up as a non-uniform distribution of sampled points, but it
    was simpler just to add extra points within the very extremes. I had to do
    this twice, getting smaller and smaller increments, to fully capture the
    extremes of pp_dnm.
    
    Args:
        increments: number of increments to span the pp_dnm range across.
        
    Returns:
        list of pp_dnm values
    """
    
    low_2 = [ x/float(increments)** 3 for x in range(increments)  ]
    low_1 = [ x/float(increments)** 2 for x in range(increments)  ]
    
    high_1 = [ 1 - x for x in low_1 ][::-1]
    high_2 = [ 1 - x for x in low_2 ][::-1]
    
    pp_dnms = [ x/float(increments) for x in range(increments) ]
    
    pp_dnms = low_2 + low_1[1:] + pp_dnms[1:-1] + high_1[:-1] + high_2
    
    return pp_dnms

def plot_curve(false_positives, true_positives, pp_dnms, threshold=0.006345):
    ''' plot a ROC curve from validation data
    
    Args:
        false_positives: list of false positive rates
        true_positives: list of true positive rates
        pp_dnms: list of pp_dnm thresholds, matched to the lists of
            true_positives and false_positives.
        threshold: pp_dnm threshold to highlight on the plot.
    '''
    
    # find the x and y-value at a given pp_dnm threshold, so that we can
    # highlight this point on the plot
    delta = [ abs(x - threshold) for x in pp_dnms ]
    pos = delta.index(min(delta))
    x_val = false_positives[pos]
    y_val = true_positives[pos]
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(false_positives, true_positives, marker='.', markersize=10)
    
    e = ax.plot(ax.get_xlim(), ax.get_ylim(), linestyle="dashed", color="gray")
    e = ax.axhline(y_val, color='red', linestyle='dashed')
    e = ax.axvline(x_val, color='red', linestyle='dashed')
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel('false positive rate')
    e = ax.set_ylabel('true positive rate')
    
    fig.savefig('de_novo_roc_curve.pdf', format='pdf',
        bbox_inches='tight', pad_inches=0, transparent=True)

def get_roc_rates(de_novos, thresholds):
    ''' get true and false positive de novo validation rates by pp_dnm threshold
    
    Args:
        de_novos: candidate variants where the validation status is known
        thresholds: list of pp_dnm thresholds to get TPR and FPR at.
    
    Returns:
        tuple of ([true positive rates], [false positive rates]) matched to the
        pp_dnm thresholds list.
    '''
    
    counts = de_novos['status'].value_counts()
    condition_positive_sum = sum(de_novos['status'])
    condition_negative_sum = sum(~de_novos['status'])
    
    true_positive_rate = []
    false_positive_rate = []
    for threshold in thresholds:
        variants = de_novos[de_novos['pp_dnm'] > threshold]
        
        true_positive = 0
        false_positive = 0
        
        if any(~variants['status']):
            false_positive = sum(~variants['status'])
    
        if any(variants['status']):
            true_positive = sum(variants['status'])
        
        true_positive_rate.append(true_positive/condition_positive_sum)
        false_positive_rate.append(false_positive/condition_negative_sum)
    
    return true_positive_rate, false_positive_rate

def plot_roc_curve_for_validations(path=PREVIOUS_COHORT, increments=100):
    """ plot a ROC curve with varying pp_dnm from de novo validation data
    
    run for all de novos attempted for validation in the 1133 trios.
    
    Args:
        path: path to de novo validation data.
        increments: number of increments to spanb across the pp_dnm range (plus
            extra at extremes).
    """
    
    de_novos = pandas.read_table(path)
    
    # recode the 'status' column, which indicates whether a variant validated
    # as a de novo or not. Remove variants where this status is unknown.
    recode = {'DNM': True, 'FP': False, 'U': None, 'INH': False, 'P/U': None, 'R': None}
    de_novos['status'] = de_novos['validation_result'].map(recode)
    de_novos = de_novos[~de_novos['status'].isnull()]
    de_novos['status'] = de_novos['status'].astype(bool)
    
    pp_dnms = get_pp_dnms(increments)
    
    true_pr, false_pr = get_roc_rates(de_novos, thresholds)
    plot_curve(true_pr, false_pr, pp_dnms)
