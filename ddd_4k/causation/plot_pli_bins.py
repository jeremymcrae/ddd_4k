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

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import seaborn

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def plot_by_hi_bin(aggregated, diff_type, output, expected=None, count_halves=False):
    """ plot difference in observed to expected across pLI vigiciles
    
    Args:
        aggregated: pandas DataFrame of observed to expected differences across
            the pLI bins.
        diff_type: list of consequence types to include, for example
            ["missense"] or ["missense", "lof"] or ["lof"]
        output: path to save the plot to.
        expected: position to plot horizontal line of expectation at, for
            example, at a ratio of 1.0.
        count_halves: whether to count the number of observed de novos in each
            half of the plot.
    """
    
    fig = pyplot.figure(figsize=(12, 6))
    ax = fig.gca()
    
    barwidth = 0.04
    e = ax.bar(aggregated["pLI_bin"], aggregated[diff_type], width=barwidth,
        align="center")
    
    # fix the axis limits and ticks
    quantiles = [ x/20.0 for x in range(20) ]
    e = ax.set_xlim((-0.05, max(quantiles) + 0.05))
    e = ax.set_xticks(quantiles)
    e = ax.set_xticklabels(quantiles, fontsize="medium")
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    
    if expected is not None:
        e = ax.axhline(expected, linestyle="dashed", color="black")
    
    if count_halves:
        e = ax.axvline(0.475, linestyle="dashed", color="black")
        lower = sum(aggregated["observed"][aggregated["pLI_bin"] <= 0.45])
        upper = sum(aggregated["observed"][aggregated["pLI_bin"] > 0.45])
        
        y_pos = ax.get_ylim()[1] * 0.9
        
        e = ax.text(0.35, y_pos, "n={0:.0f}".format(lower), horizontalalignment='center')
        e = ax.text(0.60, y_pos, "n={0:.0f}".format(upper), horizontalalignment='center')
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("normalised observed vs expected")
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0)
    pyplot.close()
