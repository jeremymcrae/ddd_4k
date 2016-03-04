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

import argparse

import numpy
import pandas
from scipy.stats import gaussian_kde

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot

import seaborn

from mupit.combine_analyses import fishersMethod

from ddd_4k.constants import THRESHOLD

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

RESULTS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.without_diagnosed.all.2015-11-24.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to .")
    parser.add_argument("--results", default=RESULTS_PATH,
        help="Path to table of consequence specific mutation rates per gene.")
    parser.add_argument("--output",
        default="results/phenotypic_similarity_comparison.pdf",
        help="Path to plot graph to.")
    
    args = parser.parse_args()
    
    return args

def combine_hpo(results):
    """ add a column to the table, for what would happen if we included the
    p-values from testing for similarity in HPO terms within DDD probands with
    de novos in each gene.
    
    Args:
        results: pandas DataFrame of results per gene.
    
    Returns:
        results DataFrame, but with extra columns for what would happen if we
        included the HPO similarity p-value.
    """
    
    results = results[~results["ddd.p_func"].isnull()]
    
    # combine the HPO similarity p-value with the functional enrichment and the
    # missense clustering.
    p_values = results[["p_min", "ddd.hpo_similarity_p_value"]]
    results["hpo_min"] = p_values.apply(fishersMethod, axis=1)
    
    return results

def plot_p_comparison(genomewide, output):
    """ plot a comparison of final p-values for genomewide genes with and
    without including the HPO similarity p-value.
    """
    
    fig = pyplot.figure(figsize=(6, 6))
    ax = fig.gca()
    
    delta = -numpy.log10(genomewide["hpo_min"]) - -numpy.log10(genomewide["p_min"])
    
    density = gaussian_kde(delta)
    x = numpy.arange(min(delta) - 2, max(delta) + 1.5, 0.01)
    e = ax.plot(x, density(x))
    
    e = ax.set_xlabel("Density")
    e = ax.set_ylabel("delta P (combined minus genotypic)")
    
    # add lines to indicate the threshold for genomewide significance
    e = ax.axvline(0, color="gray", linestyle='dashed')
    
    # fix the axis limits and ticks
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    # # plot a diagonal line of equality
    # x0, x1 = ax.get_xlim()
    # y0, y1 = ax.get_ylim()
    # lims = [max(x0, y0), min(x1, y1)]
    # e = ax.plot(lims, lims, linestyle='dashed', color='red')
    
    fig.savefig(output, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)

def main():
    args = get_options()
    
    # open the results dataset, and combine the HPO phenotype similarity p-value
    results = pandas.read_table(args.results)
    results = combine_hpo(results)
    
    # select the genes that exceed genomewide significance via either the standard
    # approach, or when combining the HPO similarity p-value
    genomewide = (results["hpo_min"] < THRESHOLD) | (results["p_min"] < THRESHOLD)
    genomewide = results[genomewide]
    
    gained = genomewide[genomewide["p_min"] >= THRESHOLD]
    lost = genomewide[genomewide["hpo_min"] >= THRESHOLD]
    print(gained[["hgnc", "p_min", "hpo_min", "ddd.hpo_similarity_p_value"]])
    print(lost[["hgnc", "p_min", "hpo_min", "ddd.hpo_similarity_p_value"]])
    
    plot_p_comparison(genomewide, args.output)

if __name__ == '__main__':
    main()
