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

# Extend the power analyses Matt did for haploinsufficient genes in the previous
# paper to compare the power of exome and genome sequencing for fixed budget,
# making some assumptions about the lower sensitivity of exome sequencing and
# the higher cost of genome sequencing - and, if we are feeling cheeky, a power
# curve for exomes run at lower cost on the X10. Rather than being based on the
# theoretical assumptions of my previous power analysis they could be based on
# the empirical observations of our 4000 trio study.

# see PLoS Genet 5(5): e1000477. doi:10.1371/journal.pgen.1000477

# the DDD 1K figure showing power to detect LoF enrichment by sample size:
# http://www.nature.com/nature/journal/v519/n7542/fig_tab/nature14135_SF6.html
# Figure method is in (search for "Loss of Function Saturation analysis"):
# http://www.nature.com/nature/journal/v519/n7542/extref/nature14135-s1.pdf

from __future__ import division, print_function

import os
import argparse

import pandas
from scipy.stats import poisson
from numpy import median

import matplotlib
matplotlib.use("Agg")

import seaborn

from ddd_4k.constants import ALPHA, NUM_GENES
from mupit.mutation_rates import get_default_rates

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(description=".")
    parser.add_argument("--output-haploinsufficiency", \
        default="haploinsufficiency_power.pdf", \
        help="Path to plot graph to.")
    parser.add_argument("--output-exome", \
        default="exome_vs_genome.simulated.pdf", \
        help="Path to plot graph to.")
    
    args = parser.parse_args()
    
    return args

def get_gene_probabilities(lof_rates, expected, threshold, cohort_n, population_n, disorder_freq):
    """ estimate probabilities of a significant result, at a cohort size
    
    Args:
        lof_rates: mutation rates per gene
        expected: list of number of mutations expected per gene in the UK +
            Ireland population.
        threshold: a significance threshold, genomewide corrected.
        cohort_n: size of the cohort being checked
        population_n: size of the population of the UK and Ireland.
        disorder_freq: expected frequency of developmental disorders
    
    Returns:
        pandas Series of probabilities for each gene
    """
    
    # estimate the expected number of mutations under the null mutation model,
    # given the cohort size.
    cohort_null = [ x * cohort_n for x in lof_rates ]
    
    # estimate individuals in UK + Ireland with developmental disorders
    disorder_n = population_n * disorder_freq
    
    # estimate the number of lof mutation per gene required to cross the
    # significance threshold, given the expected mutations in the cohort size
    mutations = [ x - 1 for x in poisson.isf(threshold, cohort_null) ]
    
    # determine the number of mutations expected for the cohort size, given the
    # proportion that the cohort size is to the number of individuals with
    # developmental disorders
    cohort_expected = [ x * cohort_n/disorder_n for x in expected ]
    
    return poisson.sf(mutations, cohort_expected)

def check_haploinsufficiency_power(rates, threshold, population_n, disorder_freq, plot_path):
    """ check power to detect lof enrichment of develop in UK
    
    Args:
        rates: dataframe of null mutation rates per gene, for different
            functional types.
        threshold: a significance threshold, genomewide corrected.
        population_n: size of the population of the UK and Ireland.
        disorder_freq: expected frequency of developmental disorders
        plot_path: path to send ouput plot to
    """
    
    # estimate the number of mutations expected per gene in the UK + Ireland population
    expected = [ x * population_n for x in rates["lof"] ]
    
    cohorts = range(1000, 13000, 1000)
    
    combined = pandas.DataFrame(columns=["hgnc", "cohort_n", "probability"])
    for cohort_n in cohorts:
        probabilities = get_gene_probabilities(rates["lof"], expected,
            threshold, cohort_n, population_n, disorder_freq)
        temp = pandas.DataFrame({"hgnc": list(rates["hgnc"]),
            "cohort_n": [cohort_n] * len(probabilities),
            "probability": list(probabilities)})
        
        combined = combined.append(temp, ignore_index=True)
    
    probabilities = get_gene_probabilities(rates["lof"], expected,
        threshold, 4295, population_n, disorder_freq)
    print(median(probabilities))
    
    combined["probability"] = combined["probability"].astype(float)
    fig = seaborn.factorplot(x="cohort_n", y="probability", data=combined, \
        kind="box", size=6, aspect=1.8, fliersize=0, color="gray")
    fig.savefig(plot_path, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)

def exome_vs_genome(rates, threshold, population_n, disorder_freq, plot_path):
    """ compare power of exome and genome sequencing
    
    Args:
        rates: dataframe of null mutation rates per gene, for different
            functional types.
        threshold: a significance threshold, genomewide corrected.
        population_n: size of the population of the UK and Ireland.
        disorder_freq: expected frequency of developmental disorders
        plot_path: path to send ouput plot to
    """
    
    # estimate the number of mutations expected per gene in the UK + Ireland
    # population
    expected = [ x * population_n for x in rates["lof"] ]
    
    # define a range of amounts of money for sequencing
    budgets = [1e6, 2e6, 3e6]
    genome_cost = 1000 * 3
    exome_relative_cost = [ (x+1)/10.0 for x in range(10) ]
    genome_sensitivity = [ (x/20.0)+1 for x in range(5) ]
    
    power = pandas.DataFrame(columns=["budget", "exome_cost", "sensitivity", "sequence", "power"])
    for budget in budgets:
        sensitivity = 1
        for relative_cost in exome_relative_cost:
            n_exomes = budget/(genome_cost * relative_cost)
            exome_probs = get_gene_probabilities(rates["lof"], expected, \
                threshold, n_exomes, population_n, disorder_freq)
            exome_median = median(exome_probs)
            power = power.append({"budget": budget,
                "exome_cost": relative_cost, "sensitivity": sensitivity,
                "sequence": "exome", "power": exome_median}, ignore_index=True)
        
        for sensitivity in genome_sensitivity:
            n_genomes = (budget/genome_cost) * sensitivity
            genome_probs = get_gene_probabilities(rates["lof"], expected, \
                threshold, n_genomes, population_n, disorder_freq)
            
            genome_median = median(genome_probs)
            power = power.append({"budget": budget,
                "exome_cost": relative_cost, "sensitivity": sensitivity,
                "sequence": "genome-{}".format(sensitivity),
                "power": genome_median}, ignore_index=True)
    
    fig = seaborn.factorplot(x="exome_cost", y="power", hue="sequence", \
        col="budget", data=power, size=8, aspect=0.4)
    fig.savefig(plot_path, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)

def main():
    
    args = get_options()
    
    # the DDD is sampling from the population of the UK and Ireland, which is about
    # 50 million people
    population_n = 50e6
    
    # the prevalence of developmental disorders is 0.5% of the general population
    disorder_freq = 0.005
    
    rates = get_default_rates()
    # determine the summed loss-of-function mutation rate for each gene
    rates["lof"] = rates[["non", "splice_site", "frameshift" ]].sum(axis=1)
    
    # estimate a genome-wide significance threshold
    threshold = ALPHA/NUM_GENES
    
    check_haploinsufficiency_power(rates, threshold, population_n, disorder_freq, \
        args.output_haploinsufficiency)
    exome_vs_genome(rates, threshold, population_n, disorder_freq, \
        args.output_exome)

if __name__ == '__main__':
    main()
