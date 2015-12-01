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

import scipy.stats
import math
import numpy

def fishers_method(values):
    """ Combine the P values from different tests into a single value.
    
    If we have only one P value for the gene for the mutation type, we just use
    that value, if we don't have any data, we use
    "NA", otherwise we combine the P values from different transcripts using
    Fisher's combined Method.
    
    Args:
        values: lists of P values from different tests
    
    Returns:
        combined p-value
    """
    
    # drop out the NA values
    values = [ x for x, is_nan in zip(values, list(numpy.isnan(values))) if not is_nan ]
    
    if len(values) == 0:
        p_value = None
    elif len(values) == 1:
        p_value = values[0]
    else:
        # use Fisher's combined method to estimate the P value from multiple
        # P values. The chi square statistic is -2*sum(ln(P values))
        values = [ math.log(x) for x in values ]
        chi_square = -2 * sum(values)
        df = 2 * len(values)
        
        # estimate the P value using the chi square statistic and degrees of
        # freedom
        p_value = scipy.stats.chi2.sf(chi_square, df)
    
    return p_value
