"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

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

import pandas

def get_count_by_person(de_novos):
    """ count the number of functional and loss-of-function de novos per person
    
    Args:
        de_novos: DataFrame of de novo mutations
    
    Returns:
        DataFrame of counts per person, split by consequence type and whether
        the mutation is in a known developmental disorder gene.
    """
    
    try:
        by_category = de_novos.pivot_table(values="chrom",
            rows=["person_stable_id", "known"],
            cols=["category"], aggfunc=len)
    except TypeError:
        by_category = de_novos.pivot_table(values="chrom",
            index=["person_stable_id", "known"],
            columns=["category"], aggfunc=len)
    
    index = [ x[0] for x in zip(by_category.index) ]
    sample_ids = [ x[0] for x in index ]
    known = [ x[2] for x in index ]
    
    by_category = pandas.DataFrame({"person_stable_id": sample_ids, "known": known, "functional": by_category["functional"].values, "loss-of-function": by_category["loss-of-function"].values})
    
    counts = pandas.melt(by_category, id_vars=["person_stable_id", "known"], value_vars=["functional", "loss-of-function"])
    counts = counts[counts["value"].notnull()]
    counts = counts.rename(columns={"variable": "consequence"})
    
    return counts
