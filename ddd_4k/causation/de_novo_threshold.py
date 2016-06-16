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

from ddd_4k.causation.count_synonymous import count_synonymous_per_gene

def get_pp_dnm_threshold(de_novos, expected):
    """ get the pp_dnm threshold where we expect as many synonymous de novos as we observed
    
    Args:
        de_novos: pandas DataFrame of all candidate coding de novos, including
            synonymous variants
        rates: pandas DataFrame of numbers of expected mutations within the
            cohort, across different categories. This includes a "synonymous_snv"
            column. The table includes all genes in the genome (or as close as
            possible), since the number of expected
    
    Returns:
        pp_dnm value, which if we exclude candidates below this value, then the
        ratio of observed synonymous variants to expected equals one. Applying
        the threshold to all candidate de novos threshold ensures we do not have
        more candidates due to errors than expected by chance in missense and
        protein-truncating categories.
    """
    
    rates = dict(zip(expected["hgnc"], expected["synonymous_snv"]))
    
    epsilon = 0.0009
    low = 0
    high = 1
    while True:
        # define a mid point to bisect at
        mid = (low + high)/2
        # select the candidates above the midpoint threshold
        candidates = de_novos[de_novos["pp_dnm"].isnull() | (de_novos["pp_dnm"] > mid) ]
        synonymous = count_synonymous_per_gene(candidates)
        
        # make sure we have expected numbers available. Some genes lack
        # expected numbers, so we exclude those genes, rather than counting
        # extra genes, for which we can never know if we are near expectations.
        synonymous["expected"] = synonymous["hgnc"].map(rates)
        synonymous = synonymous[~synonymous["expected"].isnull()]
        
        # we divide by the total expected number of synonymous variants across
        # all genes in the genome, rather than just the geens we have observed a
        # syonymous variant in.
        ratio = sum(synonymous["observed"])/sum(expected["synonymous_snv"])
        
        if ratio > 1:
            low = mid
        elif ratio < 1:
            high = mid
        
        if abs(ratio - 1) < epsilon:
            break
    
    return mid
