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
from scipy.stats import gaussian_kde
import seaborn

from ddd_4k.constants import DENOVO_PATH, VALIDATIONS, PHENOTYPES, SANGER_IDS, \
    DIAGNOSED, THRESHOLD, TRIOS
from ddd_4k.load_files import open_phenotypes, open_de_novos, get_significant_results
from ddd_4k.clinical_phenotypes import get_clinical_details, get_variant_details
from ddd_4k.ensembl_variant import EnsemblVariant

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 30, "xtick.major.size": 30})

WITH_DIAGNOSED_RESULTS = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.all.2015-11-24.txt"
WITHOUT_DIAGNOSED_RESULTS = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.without_diagnosed.all.2015-11-24.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--results-with-diagnosed",
        default=WITH_DIAGNOSED_RESULTS,
        help="Path to table of association results.")
    parser.add_argument("--results-without-diagnosed",
        default=WITHOUT_DIAGNOSED_RESULTS,
        help="Path to table of association results.")
    parser.add_argument("--de-novos", default=DENOVO_PATH,
        help="Path to table of de novo mutations.")
    parser.add_argument("--validations", default=VALIDATIONS,
        help="Path to table of de novo validation results.")
    parser.add_argument("--trios", default=TRIOS,
        help="Path to table of trios.")
    parser.add_argument("--phenotypes", default=PHENOTYPES,
        help="Path to table of phenotypes.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS,
        help="Path to table of alternate IDs for participants.")
    parser.add_argument("--diagnosed", default=DIAGNOSED,
        help="Path to table of participants with diagnoses.")
    parser.add_argument("--output-dir", default="gene_reports",
        help="Folder to send output tables to.")
    
    args = parser.parse_args()
    
    return args

def plot_genes_by_phenotype(variants, phenotypes, names, output_dir):
    """ this plots an array of boxplots (and stripplots) for each phenotype.
    """
    
    genes = sorted(variants["symbol"].unique())
    
    proband_ids = variants["person_stable_id"]
    clinical = phenotypes[phenotypes["person_stable_id"].isin(proband_ids)]
    
    # make sure we have the hgnc symbol matched to each proband ID
    recode = dict(zip(variants["person_stable_id"], variants["symbol"]))
    values = list(clinical["person_stable_id"].map(recode))
    print("bbb")
    print(len(values), len(clinical))
    clinical["symbol"] = values
    print("ccc")
    
    clinical = pandas.melt(clinical, id_vars=["symbol"], value_vars=names)
    pheno = pandas.melt(phenotypes, id_vars=["symbol"], value_vars=names)
    
    clinical = clinical.append(pheno, ignore_index=True)
    
    # exclude extreme outliers, since these are less likely to be correct values
    clinical[abs(clinical["value"]) > 20] = numpy.nan
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
        size=6, aspect=6, sharey=False)
    
    box = box.set_xticklabels(rotation=90)
    box.savefig(os.path.join(output_dir, "phenotypes_per_gene.box.pdf"), format="pdf")

def plot_gene_point(gene, variants, phenotypes, names, fig, width=0.5):
    
    gene_vars = variants[ variants["symbol"] == gene ]
    probands = gene_vars["person_stable_id"]
    
    colors = {True: "gray", False: "red"}
    color = colors[gene_vars["known"].unique()[0]]
    
    data = phenotypes[phenotypes["person_stable_id"].isin(probands)]
    
    if len(data) > 2:
        x_mean = numpy.mean(data[names[0]])
        x_sd = numpy.std(data[names[0]])
        x_vals = [x_mean - x_sd, x_mean + x_sd]
        
        y_mean = numpy.mean(data[names[1]])
        y_sd = numpy.std(data[names[1]])
        y_vals = [y_mean - y_sd, y_mean + y_sd]
        
        e = fig.gca().plot(x_vals, [y_mean, y_mean], color=color,
            linestyle="dashed", linewidth=width)
        e = fig.gca().plot([x_mean, x_mean], y_vals, color=color,
            linestyle="dashed", linewidth=width)
        e = fig.gca().plot(x_mean, y_mean, marker="o", color=color,
            markersize=5)

def plot_bivariate_contours(variants, phenotypes, names, output_dir):
    """ this plots an array of boxplots (and stripplots) for each phenotype.
    """
    
    prev_genes = sorted(variants["symbol"][variants["known"]].unique())
    new_genes = sorted(variants["symbol"][~variants["known"]].unique())
    
    # exclude extreme outliers, since these are less likely to be correct values
    for col in names:
        values = phenotypes[col]
        values[abs(values) > 20] = numpy.nan
        phenotypes[col] = values
    
    phenotypes = phenotypes[["person_stable_id"] + names].dropna()
    
    ax = seaborn.kdeplot(phenotypes[names[0]], phenotypes[names[1]],
        cmap="Blues", shade=True, shade_lowest=False, n_levels=15)
    fig = ax.get_figure()
    
    for gene in prev_genes:
        plot_gene_point(gene, variants, phenotypes, names, fig)
    
    for gene in new_genes:
        plot_gene_point(gene, variants, phenotypes, names, fig)
    
    lim = fig.gca().set_xlim([-5, 4])
    lim = fig.gca().set_ylim([-7, 4])
    fig.savefig(os.path.join(output_dir, "phenotypes.{}.pdf".format("_".join(names))), format="pdf")
    
    pyplot.close()

def main():
    
    args = get_options()
    
    args.output_dir = ""
    
    de_novos = open_de_novos(args.de_novos, args.validations)
    phenotypes = open_phenotypes(args.phenotypes, args.sanger_ids)
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    trios = pandas.read_table(args.trios, sep="\t")
    phenotypes = phenotypes[phenotypes["decipher_id"].isin(trios["decipher_id"])]
    
    with_genomewide = get_significant_results(args.results_with_diagnosed, THRESHOLD)
    without_genomewide = get_significant_results(args.results_without_diagnosed, THRESHOLD)
    
    # Remove genes that genomewide significant in the with diagnoses results, but
    # not in the without diagnoses set, and are not known DD genes. When we remove
    # probands with diagnoses, probands in these genes are removed, which is why
    # they do not appear to be genomewide significant in the results without
    # diagnoses. We only want to include genes that are genomewide, and either
    # known, or novel.
    genes_to_remove = with_genomewide["hgnc"][~with_genomewide["hgnc"].isin(without_genomewide["hgnc"]) & ~with_genomewide["in_ddg2p"]]
    
    genes = set(with_genomewide["hgnc"]) | set(without_genomewide["hgnc"])
    genes -= set(genes_to_remove)
    hgnc = list(with_genomewide["hgnc"]) + list(without_genomewide["hgnc"])
    in_ddg2p = list(with_genomewide["in_ddg2p"]) + list(without_genomewide["in_ddg2p"])
    status = dict(zip(hgnc, in_ddg2p))
    
    # remove the genes we want to exclude, since they don't occur in the without
    # diagnosed results due to some of the probands having alternate, diagnostic
    # de novos.
    for x in genes_to_remove:
        del status[x]
    
    variants = de_novos[de_novos["symbol"].isin(genes)]
    variants["known"] = variants["symbol"].map(status)
    
    phenotypes["symbol"] = "full DDD"
    
    names = ["birthweight_sd", "height_sd", "weight_sd", "ofc_sd"]
    plot_genes_by_phenotype(variants, phenotypes, names, args.output_dir)
    pyplot.close()
    
    names = ["height_sd", "ofc_sd"]
    plot_bivariate_contours(variants, phenotypes, names, args.output_dir)

if __name__ == '__main__':
    main()
