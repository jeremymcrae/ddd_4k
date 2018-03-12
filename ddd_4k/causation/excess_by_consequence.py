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

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

from ddd_4k.causation.merging import merge_observed_and_expected

def get_consequence_excess(expected, de_novos, ppv, sensitivity):
    """ determine the ratio of observed to expected for synonymous, missense
    and protein-truncating candidate de novos
    
    Args:
        expected: pandas DataFrame of counts of expected de novo mutations for
            all genes (or nearly all) in the genome.
        de_novos: pandas DataFrame of candidate de novo mutations observed in
            the cohort.
        ppv: positive predictive value at a given de novo quality threshold.
        sensitivity: sensitivity to true positives at a given de novo quality
            threshold.
    
    Returns:
        dictionary of counts and ratios for synonymous, missense and
        truncating consequences.
    """
    
    # adjust the numbers of excess mutations by our sensitivity and specificity,
    # which are derived from validation data in an earlier dataset. The ppv
    # indicates the proportion of candidates which are true positive at a given
    # quality threshold, and the sensitivity indicates the sensitivity to true
    # positives at the same threshold. We also need to account for exome sequencing
    # not being 100% sensitive, as estimates are that it is around 95% sensitive.
    scale_factor = ppv * (1.0/0.95/sensitivity)
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    synonymous_ratio = 1.0 * scale_factor
    synonymous_ratio = sum(merged["synonymous_observed"])/sum(merged["synonymous_expected"]) * scale_factor
    missense_ratio = sum(merged["missense_observed"])/sum(merged["missense_expected"]) * scale_factor
    lof_ratio = sum(merged["lof_observed"])/sum(merged["lof_expected"]) * scale_factor
    
    synonymous_count = len(de_novos[de_novos["consequence"] == "synonymous_variant"]) * scale_factor
    lof_count = sum(merged["lof_observed"]) * scale_factor
    missense_count = sum(merged["missense_observed"]) * scale_factor
    
    synonymous_excess = synonymous_count * ((synonymous_ratio - 1)/synonymous_ratio)
    missense_excess = missense_count - sum(merged["missense_expected"]) * scale_factor
    lof_excess = lof_count - sum(merged["lof_expected"]) * scale_factor
    
    values = {
        "synonymous":
            {"ratio": synonymous_ratio,  "count": synonymous_count,  "excess": synonymous_excess},
        "truncating":
            {"ratio": lof_ratio,  "count": lof_count, "excess": lof_excess},
        "missense":
            {"ratio": missense_ratio,  "count": missense_count,  "excess": missense_excess}}
    
    return values

def plot_consequence_excess(ratios, output):
    """ plot the ratio of observed to expected de novo counts for overall consequences
    
    Args:
        ratios: dictionary of de novo counts and observed/expected ratios for
            "synonymous", "truncating" and "missense" candidate de novos.
        output: path to save plot as pdf to.
    """
    
    # format the ratios and counts
    temp = pandas.DataFrame(ratios).transpose()
    temp["consequence"] = temp.index
    temp = temp.reindex(["synonymous", "missense", "truncating"])
    temp.index = range(len(temp))
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    e = ax.bar(range(len(temp)), temp["ratio"], align="center")
    
    # annotate the plot, to show the baseline, and the numbers of candidate de
    # novos in each category
    e = ax.axhline(1.0, color="black", linestyle="dashed")
    for key, row in temp.iterrows():
        e = ax.text(key, row["ratio"]+0.01,
            "n={0:.0f}\nexcess={1:.0f}".format(row["count"], row["excess"]),
            horizontalalignment='center')
    
    # fix the axis limits and ticks
    e = ax.set_xlim((-0.5, len(temp) - 0.5))
    e = ax.set_xticks(range(len(temp)))
    e = ax.set_xticklabels(temp["consequence"])
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    e = ax.set_xlabel("consequence class")
    e = ax.set_ylabel("observed/expected")
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0)
