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

from __future__ import division, print_function

import argparse
import getpass
import sys
import os
import math

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

from matplotlib import pyplot
import numpy
import pandas
import psycopg2
from scipy.stats import gaussian_kde
import seaborn

from ddd_4k.constants import PHENOTYPES, SANGER_IDS, DIAGNOSED
from ddd_4k.load_files import open_phenotypes
from ddd_4k.clinical_phenotypes import get_clinical_details, get_variant_details
from ddd_4k.ensembl_variant import EnsemblVariant

# define the plot style
seaborn.set_context("notebook", font_scale=6)
seaborn.set_style("white", {"ytick.major.size": 30, "xtick.major.size": 30})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--de-novos", \
        help="Path to table of variants in novel genes.")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to table of phenotypes.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS, \
        help="Path to table of alternate IDs for participants.")
    parser.add_argument("--diagnosed", default=DIAGNOSED, \
        help="Path to table of participants with diagnoses.")
    parser.add_argument("--output-dir", default="gene_reports", \
        help="Folder to send output tables to.")
    
    args = parser.parse_args()
    
    return args

def connect_to_database():
    """ connect to the DDD database
    
    This raises an error if the wrong password is entered
    """
    
    tries = 0
    while tries < 3:
        try:
            pwd = getpass.getpass("DDD database password: ")
            conn = psycopg2.connect(database="ddd_prod", user="ddd_login_ro", \
                    host="ddd-lims-db", port="5444", password=pwd)
            break
        except psycopg2.OperationalError:
            tries += 1
    
    if tries == 3:
        sys.exit("too many wrong password attempts")
    
    return conn

def x_samples(x, y, jump):
    while x < y:
        yield x
        x += jump

def plot_phenotypes(gene_rows, phenotypes, fig, spec):
    """
    """
    
    to_plot = {"birthweight_sd": {"color": "red", "name": "birthweight"},
               "height_sd": {"color": "green", "name": "height"},
               "weight_sd": {"color": "blue", "name": "weight"},
               "ofc_sd": {"color": "orange", "name": "OFC"} }
    
    subgrid = gridspec.GridSpecFromSubplotSpec(len(to_plot) + 1, 1, subplot_spec=spec)
    x_min = -10
    x_max = 10
    jump = (x_max - x_min)/50
    x_values = list(x_samples(x_min, x_max, jump))
    
    for pheno in sorted(to_plot):
        
        color = to_plot[pheno]["color"]
        
        ax = pyplot.Subplot(fig, subgrid[sorted(to_plot).index(pheno), :])
        fig.add_subplot(ax)
        
        # generate kernel densities for the full DDD population, plot the
        # population density,
        population_kde = gaussian_kde(phenotypes[pheno][~phenotypes[pheno].isnull()])
        population_y_values = list(population_kde(x_values))
        pop = ax.fill_between(x_values, population_y_values, facecolor="gray", edgecolor="gray", alpha=0.5)
        
        try:
            # if there are enough non-null values, plot a kernel density of the
            # values for the probands for the gene
            gene_kde = gaussian_kde(gene_rows[pheno][~gene_rows[pheno].isnull()])
            gene_y_values = list(gene_kde(x_values))
            gene = ax.fill_between(x_values, gene_y_values, facecolor=color, edgecolor=color, alpha=0.5)
        except ValueError:
            continue
        
        # set axis parameters
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')
        ticks = ax.yaxis.set_ticks(numpy.arange(min(x), max(x)+1))
        # ticks = ax.yaxis.set_ticks([0, 0.1, 0.2, 0.3, 0.4])
        # lab = ax.set_ylabel(to_plot[pheno]["name"])
        
        # only plot an x-axis for the last graph
        # if pheno == sorted(pheno)[-1]:
        #     ax.xaxis.set_ticks_position('bottom')

def plot_genes_by_phenotype(clinical_table, phenotypes, names):
    """ this plots an array of boxplots (and stripplots) for each phenotype.
    """
    
    genes = sorted(clinical_table["symbol"].unique())
    
    clinical = clinical_table[["symbol"] + names].append(
        phenotypes[["symbol"] + names], ignore_index=True)
    
    clinical = pandas.melt(clinical, id_vars=["symbol"], value_vars=names)
    
    # exclude extreme outliers, since these are less likely to be correct values
    clinical = clinical[abs(clinical["value"]) < 20]
    clinical = clinical.dropna()
    
    # the full DDD points span a much wider range than any other set. Since we
    # don't plot the DDD points, we don't want their full range shown in the
    # plots. Shrink the boxplot outlier values to just at the boxplot edge.
    for x in names:
        values = clinical["value"][(clinical["symbol"] == "full DDD") & (clinical["variable"] == x)]
        upper = values.quantile(0.75)
        lower = values.quantile(0.25)
        iqr = upper - lower
        values[(values > upper + 1.5 * iqr)] = upper
        values[(values < lower - 1.5 * iqr)] = lower
        clinical["value"][(clinical["symbol"] == "full DDD") & (clinical["variable"] == x)] = values
    
    # set the DDD values as None, so we don't get any points plotted on the
    # stripplot. We will still have the boxplot to show the distribution.
    clin_without = clinical.copy()
    clin_without["value"][clin_without["symbol"] == "full DDD"] = None
    
    kwargs = {"size": 15}
    strip = seaborn.factorplot(x="symbol", y="value", row="variable", data=clin_without,
        order=["full DDD"] + genes, kind="strip", alpha=0.5, jitter=True,
        size=6, aspect=6, sharey=False)
    
    strip = strip.set_xticklabels(rotation=90)
    strip.savefig("phenotypes_per_gene.strip.pdf", format="pdf")
    
    box = seaborn.factorplot(x="symbol", y="value", row="variable",
        data=clinical, order=["full DDD"] + genes, kind="box", fliersize=0,
        size=6, aspect=6, sharey=False)
    
    box = box.set_xticklabels(rotation=90)
    box.savefig("phenotypes_per_gene.box.pdf", format="pdf")

def main():
    
    args = get_options()
    args.de_novos = "/lustre/scratch113/projects/ddd/users/jm33/results/novel_gene_variants.ddd_4k.2015-11-24.txt"
    
    variants = pandas.read_table(args.de_novos, sep="\t")
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    
    phenotypes["symbol"] = "full DDD"
    
    # remove the variants in diagnosed individuals, since they didn't contribute
    # to identify the gene as being genomewide significant. These "diagnosed"
    # individuals still need to be inspected, to make sure they don't conflict
    # with the evidence from the other probands.
    variants = variants[~variants["person_stable_id"].isin(diagnosed["person_id"])]
    
    conn = connect_to_database()
    
    with conn.cursor() as cur:
        
        ensembl = EnsemblVariant(cache_folder="cache", genome_build="grch37")
        var_table = get_variant_details(variants, ensembl, cur)
        clinical_table = get_clinical_details(var_table, phenotypes)
        
        genes = sorted(clinical_table["symbol"].unique())
        size = math.ceil(math.sqrt(len(genes)))
        
        names = ["birthweight_sd", "height_sd", "weight_sd", "ofc_sd"]
        plot_genes_by_phenotype(clinical_table, phenotypes, names)
        
        fig = plt.figure()
        grid = gridspec.GridSpec(int(size), int(size))
        
        pos = 0
        for gene in genes:
            print(gene)
            gene_rows = clinical_table[clinical_table["symbol"] == gene]
            plot_phenotypes(gene_rows, phenotypes, fig, grid[pos])
            pos += 1
        
        fig.savefig("phenotypes_per_gene.alternate.pdf", format="pdf")

if __name__ == '__main__':
    main()
