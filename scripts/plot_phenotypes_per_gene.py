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

def plot_genes_by_phenotype(clinical_table, phenotypes, names, output_dir):
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
    
    strip = seaborn.factorplot(x="symbol", y="value", row="variable", data=clin_without,
        order=["full DDD"] + genes, kind="strip", alpha=0.5, jitter=True,
        size=6, aspect=6, sharey=False)
    
    strip = strip.set_xticklabels(rotation=90)
    strip.savefig(os.path.join(output_dir, "phenotypes_per_gene.strip.pdf"), format="pdf")
    
    box = seaborn.factorplot(x="symbol", y="value", row="variable",
        data=clinical, order=["full DDD"] + genes, kind="box", fliersize=0,
        size=6, aspect=6, sharey=False")
    
    box = box.set_xticklabels(rotation=90)
    box.savefig(os.path.join(output_dir, "phenotypes_per_gene.box.pdf"), format="pdf")

def main():
    
    args = get_options()
    
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
        plot_genes_by_phenotype(clinical_table, phenotypes, names, args.output_dir)

if __name__ == '__main__':
    main()
