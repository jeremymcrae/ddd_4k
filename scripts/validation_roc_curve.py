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

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

from ddd_4k.load_files import open_de_novos
from ddd_4k.validation_rates import open_previous_validations, get_rates
from ddd_4k.constants import DENOVO_PATH, PREVIOUS_VALIDATIONS

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to probe the limits " \
        "of de novo causation in children with developmental disorders.")
    parser.add_argument("--current-de-novos", default=DENOVO_PATH,
        help="Path to table of candidate de novo mutations.")
    parser.add_argument("--validations", default=PREVIOUS_VALIDATIONS,
        help="Path to table of candidate de novo mutations.")
    parser.add_argument("--increments", default=100,
        help="Path to table of validation data.")
    parser.add_argument("--threshold", default=0.0078125,
        help="PP_DNM threshold to label on the plot.")
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

def get_roc_rates(de_novos, thresholds):
    ''' get true and false positive de novo validation rates by pp_dnm threshold
    
    Args:
        de_novos: candidate variants where the validation status is known
        thresholds: list of pp_dnm thresholds to get TPR and FPR at.
    
    Returns:
        tuple of ([true positive rates], [false positive rates]) matched to the
        pp_dnm thresholds list.
    '''
    
    true_positive_rate = []
    false_positive_rate = []
    ppv = []
    npv = []
    for threshold in thresholds:
        rates = get_rates(de_novos, threshold)
        
        true_positive_rate.append(rates.tpr)
        false_positive_rate.append(rates.fpr)
        ppv.append(rates.ppv)
        npv.append(rates.npv)
    
    return ppv, npv, true_positive_rate, false_positive_rate

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

def main():
    """ plot a ROC curve with varying pp_dnm from de novo validation data
    
    run for all de novos attempted for validation in the 1133 trios.
    """
    
    args = get_options()
    
    current = open_de_novos(args.current_de_novos, exclude_synonymous=False)
    de_novos = open_previous_validations(args.validations, current)
    
    pp_dnms = get_pp_dnms(args.increments)
    
    ppv, npv, true_pr, false_pr = get_roc_rates(de_novos,  pp_dnms)
    plot_curve(ppv, true_pr, pp_dnms, threshold=args.threshold,
        output=args.output, xlabel='positive predictive value',
        ylabel='true positive rate')
    
    plot_curve(false_pr, true_pr, pp_dnms, threshold=args.threshold,
        output='de_novo_roc_curve.fpr_vs_tpr.pdf', xlabel='false positive rate',
        ylabel='true positive rate')

if __name__ == '__main__':
    main()
