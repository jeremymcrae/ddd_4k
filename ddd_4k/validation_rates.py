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

from ddd_4k.causation.excess_by_consequence import get_consequence_excess

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def open_previous_validations(path, current=None):
    ''' open the validation data from the previous datafreeze.
    
    Args:
        path: path to de novo validation data.
        
    Returns:
        a pandas DataFrame, containing de novo validation data, where the
        'status' column indicates whether the de novo validation was successful
        or not.
    '''
    
    de_novos = pandas.read_table(path)
    
    # recode the 'status' column, which indicates whether a variant validated
    # as a de novo or not. Remove variants where this status is unknown.
    recode = {'DNM': True, 'FP': False, 'U': None, 'INH': False, 'P/U': None, 'R': None}
    de_novos['status'] = de_novos['validation_result'].map(recode)
    de_novos = de_novos[~de_novos['status'].isnull()]
    de_novos['status'] = de_novos['status'].astype(bool)
    
    if current is not None:
        keys = current['person_stable_id'] + current['chrom'] + current['pos'].astype(str)
        de_novos['key'] = de_novos['person_stable_id'] + de_novos['chr'] + de_novos['pos'].astype(str)
        
        de_novos = de_novos[de_novos['key'].isin(set(keys))]
    
    return de_novos

def get_rates(de_novos, threshold):
    ''' estimate the true positive rate, false positive rate, positive
    predictive value and negative predictive value from de novo validation data
    at a given pp_dnm threshold.
    
    Args:
        de_novos: pandas DataFrame of de novo validarion data, one row per
            candidate
        threshold: pp_dnm threshold to filter on.
        
    Returns:
        Rates object, containing tpr, fpr, ppv and npv for the given threshold.
    '''
    
    variants = de_novos[de_novos['pp_dnm'] > threshold]
    
    condition_positive_sum = sum(de_novos['status'])
    condition_negative_sum = len(de_novos) - condition_positive_sum
    
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
    
    try:
        tpr = true_positive/(true_positive + false_negative)
    except ZeroDivisionError:
        tpr = 0
    
    try:
        fpr = false_positive/(false_positive + true_negative)
    except ZeroDivisionError:
        fpr = 0
    
    try:
        ppv = true_positive/(true_positive + false_positive)
    except ZeroDivisionError:
        ppv = None
    try:
        npv = true_negative/(false_negative + true_negative)
    except ZeroDivisionError:
        npv = None
        
    class Rates:
        def __init__(self, tpr, fpr, ppv, npv):
            self.tpr = tpr
            self.fpr = fpr
            self.ppv = ppv
            self.npv = npv
    
    return Rates(tpr, fpr, ppv, npv)

def plot_excess_by_pp_dnm_threshold(de_novos, expected, validations, increments=100):
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
        validation_rate = get_rates(validations, threshold)
        excess = get_consequence_excess(expected, variants, validation_rate.ppv, validation_rate.tpr)
        
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
