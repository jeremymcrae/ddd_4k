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
    
    # sum the observed loss-of-function and missense candidates per gene
    observed = get_de_novo_counts(de_novos)
    observed["lof_observed"] = observed[["lof_snv", "lof_indel"]].sum(axis=1)
    observed["missense_observed"] = observed[["missense_snv", "missense_indel"]].sum(axis=1)
    
    # restrict the observed counts to the columns for merging.
    observed = observed[["hgnc", "lof_observed", "missense_observed"]]
    
    # merge the observed counts with the expected counts, for genes without
    # counts to merge, replace NA values with 0.
    expected = expected[["hgnc", "lof_indel", "lof_snv", "missense_indel",
        "missense_snv", "synonymous_snv"]].copy()
    merged = expected.merge(observed, how="left", on="hgnc")
    merged["lof_observed"][merged["lof_observed"].isnull()] = 0
    merged["missense_observed"][merged["missense_observed"].isnull()] = 0
    
    # get the numbers of expected loss-of-function, missense and
    # synonymous mutations.
    merged["lof_expected"] = merged[["lof_snv", "lof_indel"]].sum(axis=1)
    merged["missense_expected"] = merged[["missense_snv", "missense_indel"]].sum(axis=1)
    merged = merged.rename(columns={'synonymous_snv': 'synonymous_expected'})
    
    # restrict the table to the minimum required columns
    merged = merged[["hgnc", "lof_observed", "missense_observed",
        "lof_expected", "missense_expected", "synonymous_expected"]].copy()
    
    return merged