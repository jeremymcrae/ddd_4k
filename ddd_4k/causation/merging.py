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

from mupit.count_de_novos import get_de_novo_counts

def merge_observed_and_expected(de_novos, expected):
    """ get per gene counts of observed and expected de novo mutations
    
    Args:
        expected: pandas DataFrame of counts of expected de novo mutations for
            all genes (or nearly all) in the genome.
        de_novos: pandas DataFrame of candidate de novo mutations observed in
            the cohort.
    
    Returns:
        pandas DataFrame of rows per gene of observed and expected counts for
        loss-of-function and missense candidates.
    """
    
    # sum the observed truncating and missense candidates per gene
    observed = get_de_novo_counts(de_novos)
    observed["lof_observed"] = observed[["lof_snv", "lof_indel"]].sum(axis=1)
    observed["missense_observed"] = observed[["missense_snv", "missense_indel"]].sum(axis=1)
    
    # sum the expected truncating, missense and synonymous mutations per gene
    expected = expected.copy()
    expected["lof_expected"] = expected[["lof_snv", "lof_indel"]].sum(axis=1)
    expected["missense_expected"] = expected[["missense_snv", "missense_indel"]].sum(axis=1)
    expected['synonymous_expected'] = expected['synonymous_snv']
    
    # merge the observed counts with the expected counts, for genes without
    # counts to merge, replace NA values with 0.
    merged = expected.merge(observed, how="left", on="hgnc")
    merged["lof_observed"][merged["lof_observed"].isnull()] = 0
    merged["missense_observed"][merged["missense_observed"].isnull()] = 0
    
    return merged
