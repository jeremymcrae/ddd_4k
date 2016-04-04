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

import argparse

import pandas
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

PREVIOUS_COHORT = '/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv'

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to probe the limits " \
        "of de novo causation in children with developmental disorders.")
    parser.add_argument("--validations", default=PREVIOUS_COHORT,
        help="Path to table of candidate de novo mutations.")
    parser.add_argument("--increments", default=100,
        help="Path to table of validation data.")
    parser.add_argument("--output", default='de_novo_roc_curve.pdf',
        help="Path to write pdf plot to.")
    
    args = parser.parse_args()
    
    return args

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

def plot_curve(x_values, y_values, pp_dnms, threshold=0.006345,
        output=None, xlabel='', ylabel=''):
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
    x_val = x_values[pos]
    y_val = y_values[pos]
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.plot(x_values, y_values, marker='.', markersize=10)
    
    # e = ax.plot(ax.get_xlim(), ax.get_ylim(), linestyle="dashed", color="gray")
    e = ax.axhline(y_val, color='red', linestyle='dashed')
    e = ax.axvline(x_val, color='red', linestyle='dashed')
    
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel(xlabel)
    e = ax.set_ylabel(ylabel)
    
    fig.savefig(output, format='pdf', bbox_inches='tight', pad_inches=0, transparent=True)

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
    ppv = []
    npv = []
    for threshold in thresholds:
        variants = de_novos[de_novos['pp_dnm'] > threshold]
        
        true_positive = 0
        false_positive = 0
        false_negative = 0
        true_negative = 0
        
        if any(~variants['status']):
            false_positive = sum(~variants['status'])
            true_negative = condition_negative_sum - false_positive
    
        if any(variants['status']):
            true_positive = sum(variants['status'])
            false_negative = condition_positive_sum - true_positive
        
        true_positive_rate.append(true_positive/condition_positive_sum)
        false_positive_rate.append(false_positive/condition_negative_sum)
        try:
            ppv.append(true_positive/(true_positive + false_positive))
            npv.append(true_negative/(false_negative + true_negative))
        except ZeroDivisionError:
            ppv.append(None)
            npv.append(None)
    
    return ppv, npv, true_positive_rate, false_positive_rate

def main():
    """ plot a ROC curve with varying pp_dnm from de novo validation data
    
    run for all de novos attempted for validation in the 1133 trios.
    """
    
    args = get_options()
    
    de_novos = pandas.read_table(args.validations)
    
    # recode the 'status' column, which indicates whether a variant validated
    # as a de novo or not. Remove variants where this status is unknown.
    recode = {'DNM': True, 'FP': False, 'U': None, 'INH': False, 'P/U': None, 'R': None}
    de_novos['status'] = de_novos['validation_result'].map(recode)
    de_novos = de_novos[~de_novos['status'].isnull()]
    de_novos['status'] = de_novos['status'].astype(bool)
    
    pp_dnms = get_pp_dnms(args.increments)
    
    ppv, npv, true_pr, false_pr = get_roc_rates(de_novos,  pp_dnms)
    plot_curve(ppv, true_pr, pp_dnms, threshold=0.0078125,
        output=args.output, xlabel='false positive rate', ylabel='true positive rate')

if __name__ == '__main__':
    main()
