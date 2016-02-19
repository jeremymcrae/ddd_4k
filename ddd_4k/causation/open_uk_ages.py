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

from __future__ import absolute_import

import pandas

def open_uk_parent_ages(path):
    """ gets tables of numbers of UK-wide father and mothers by age
    
    Args:
        path: path to excel spreadsheet containing UK-wide counts for mother's
            and father's by age.
    
    Returns:
        tuple of pandas DataFrames, one for fathers ages, one for mother's ages
    """
    
    ages = pandas.read_excel(path, sheetname="Mothers_Fathers", skiprows=1)
    ages.columns = ['fathers_age_raw', 'fathers_count_raw',
        'fathers_cumulative_frequency', 'mothers_age_raw', 'mothers_count_raw',
        'fathers_cumulative_frequency', 'none', 'none', 'age',
        'fathers_count', 'none', 'mothers_age', 'mothers_count']
    
    return ages[['age', 'fathers_count', 'mothers_count']].copy()
