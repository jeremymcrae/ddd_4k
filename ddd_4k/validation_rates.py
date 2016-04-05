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
