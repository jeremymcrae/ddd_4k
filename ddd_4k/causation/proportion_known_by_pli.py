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
from matplotlib import pyplot, gridspec
import seaborn

from ddd_4k.causation.merging import merge_observed_and_expected
from ddd_4k.causation.classify_known_genes import classify_monoallelic_genes
from ddd_4k.causation.pli_functions import include_constraints, get_constraint_bins

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def plot_proportion_known_by_pLI(de_novos, expected,  constraints, known):
    """ plot the proportion of genes known to be dominant developmental genes
    at different pLI bins.
    """
    
    merged = merge_observed_and_expected(de_novos, expected)
    
    # identify known haploinsufficient dominant genes
    mono = classify_monoallelic_genes(known)
    
    # identify which pLI quantile each gene falls into
    merged = include_constraints(merged, constraints)
    merged["pLI_bin"] = get_constraint_bins(merged, bins=20)
    
    merged["haploinsufficient"] = merged["hgnc"].isin(mono["haploinsufficient"])
    merged["nonhaploinsufficient"] = merged["hgnc"].isin(mono["nonhaploinsufficient"])
    
    plot1(merged)
    plot2(merged)

def plot1(merged):
    """
    """
    
    hi = merged.pivot_table(index="pLI_bin", values="haploinsufficient", aggfunc=[sum, len])
    non_hi = merged.pivot_table(index="pLI_bin", values="nonhaploinsufficient", aggfunc=[sum, len])
    
    hi["ratio"] = hi["sum"]/hi["len"]
    non_hi["ratio"] = non_hi["sum"]/non_hi["len"]
    
    hi["pLI_bin"] = hi.index
    non_hi["pLI_bin"] = non_hi.index
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    e = ax.plot(hi["pLI_bin"], hi["ratio"])
    e = ax.plot(non_hi["pLI_bin"], non_hi["ratio"])
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("proportion known pathogenic")
    
    fig.savefig("proportion_XXX1.pdf", format="pdf", bbox_inches='tight', pad_inches=0)

def plot2(merged):
    """
    """
    
    hi_lof = merged[merged["haploinsufficient"]].pivot_table(index="pLI_bin",
        values="lof_observed", aggfunc=sum)
    non_hi_mis = merged[merged["nonhaploinsufficient"] | \
        merged["haploinsufficient"]].pivot_table(index="pLI_bin",
            values="missense_observed", aggfunc=sum)
    
    lof = merged.pivot_table(index="pLI_bin", values="lof_observed", aggfunc=sum)
    mis = merged.pivot_table(index="pLI_bin", values="missense_observed", aggfunc=sum)
    
    known_proportions = pandas.DataFrame({"pLI_bin": lof.index,
        "lof_ratio": hi_lof/lof, "mis_ratio": non_hi_mis/mis})
    
    fig = pyplot.figure()
    ax = fig.gca()
    
    e = ax.plot(known_proportions["pLI_bin"], known_proportions["lof_ratio"])
    e = ax.plot(known_proportions["pLI_bin"], known_proportions["mis_ratio"])
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('bottom')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("pLI bin (low to high)")
    e = ax.set_ylabel("proportion known pathogenic")
    
    fig.savefig("proportion_XXX2.pdf", format="pdf", bbox_inches='tight', pad_inches=0)
