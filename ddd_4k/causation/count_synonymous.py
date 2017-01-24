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

import pandas

def count_synonymous_per_gene(de_novos):
    """ count the number of synonymous mutations per gene
    
    Args:
        de_novos: pandas DataFrame containing one row per candidate de novo
            mutation.
    
    Returns:
        pandas DataFrame of counts of synonymous variants per gene
    """
    
    synonymous = de_novos[de_novos["consequence"] == "synonymous_variant"]
    counts = synonymous.pivot_table(index="symbol", values="alt", aggfunc=len)
    counts = pandas.DataFrame({"hgnc": counts.index, "observed": counts})
    
    counts.index = range(len(counts))
    
    return counts
