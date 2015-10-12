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

import urllib
import math
import tempfile

import pandas
from scipy.stats import poisson
from numpy import median

import matplotlib
matplotlib.use("Agg")

import seaborn

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})


RATES_URL = "http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S2.xls"

def get_mutation_rates(url):
    """ obtain the table of mutation rates from the Samocha et al paper
    
    Args:
        url: url to supplementary mutation rates table
    
    Returns:
        dataframe of mutation rates, with an extra column for summed lof rate
    """
    
    temp = tempfile.NamedTemporaryFile()
    # get a list of lof mutation rates per gene
    urllib.urlretrieve(url, temp.name)
    rates = pandas.read_excel(temp.name, sheetname="mutation_probabilities")
    
    # convert the log10 adjusted rates back to unscaled numbers
    rates["non"] = [ 10**x for x in rates["non"] ]
    rates["splice_site"] = [ 10**x for x in rates["splice_site"] ]
    rates["frameshift"] = [ 10**x for x in rates["frameshift"] ]
    
    # determine the summed loss-of-function mutation rate for each gene
    rates["lof"] = rates[["non", "splice_site", "frameshift" ]].sum(axis=1)
    
    return rates

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
        probabilities = get_gene_probabilities(rates["lof"], expected, threshold, cohort_n, population_n, disorder_freq)
        temp = pandas.DataFrame({"hgnc": list(rates["gene"]),
            "cohort_n": [cohort_n] * len(probabilities),
            "probability": list(probabilities)})
        
        combined = combined.append(temp, ignore_index=True)
    
    combined["probability"] = combined["probability"].astype(float)
    fig = seaborn.factorplot(x="cohort_n", y="probability", data=combined, kind="box", size=6, aspect=1.8, fliersize=0, color="gray")
    fig.savefig(plot_path, format="pdf")

def exome_vs_genome(rates, threshold, population_n, disorder_freq, plot_path):
    """
    """
    
    # estimate the number of mutations expected per gene in the UK + Ireland population
    expected = [ x * population_n for x in rates["lof"] ]
    
    # define a range of amounts of money for sequencing
    budgets = [1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8]
    genome_cost = 1000
    exome_relative_cost = [ x/10.0 for x in range(1, 11) ]
    genome_sensitivity = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    
    power = pandas.DataFrame(columns=["budget", "exome_cost", "sensitivity", "sequence_type", "power"])
    for budget in budgets:
        for relative_cost in exome_relative_cost:
            for sensitivity in genome_sensitivity:
                n_genomes = budget/genome_cost
                n_exomes = budget/(genome_cost * relative_cost)
                print(budget, relative_cost, sensitivity, n_genomes, n_exomes)
                
                exome_expected = [ x/(max(genome_sensitivity)/sensitivity) for x in expected ]
                
                genome_probs = get_gene_probabilities(rates["lof"], expected, \
                    threshold, n_genomes, population_n, disorder_freq)
                exome_probs = get_gene_probabilities(rates["lof"], exome_expected, \
                    threshold, n_genomes, population_n, disorder_freq)
                
                genome_median = median(genome_probs)
                exome_median = median(exome_probs)
                print(genome_median, exome_median)
                
                power = power.append({"budget": budget,
                    "exome_cost": relative_cost, "sensitivity": sensitivity,
                    "sequence_type": "genome", "power": genome_median}, ignore_index=True)
                power = power.append({"budget": budget,
                    "exome_cost": relative_cost, "sensitivity": sensitivity,
                    "sequence_type": "exome", "power": exome_median}, ignore_index=True)
    
    fig = seaborn.factorplot(x="exome_cost", y="power", hue="sequence_type", col="budget", row="sensitivity", data=power)
    fig.savefig("test_exome_vs_genome.pdf", format="pdf")

def main():
    # the DDD is sampling from the population of the UK and Ireland, which is about
    # 50 million people
    population_n = 50e6
    
    # the prevalence of developmental disorders is 0.5% of the general population
    disorder_freq = 0.005
    
    rates = get_mutation_rates(RATES_URL)
    
    # estimate a genome-wide significance threshold
    threshold = 0.05/18500
    
    check_haploinsufficiency_power(rates, threshold, population_n, disorder_freq, "test.pdf")
    # exome_vs_genome(rates, threshold, population_n, disorder_freq, "test_exome_vs_genome.pdf")

if __name__ == '__main__':
    main()
