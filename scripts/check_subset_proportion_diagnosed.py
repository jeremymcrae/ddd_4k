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
import os

import pandas

import matplotlib
matplotlib.use("Agg")
import seaborn

from ddd_4k.constants import KNOWN_GENES, DENOVO_PATH, VALIDATIONS, DIAGNOSED
from ddd_4k.load_files import open_de_novos

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check the impact"
        "of including different subsets on the ability to detect genes involved"
        "in developmental disorders.")
    parser.add_argument("--known-genes", default=KNOWN_GENES,
        help="Path to table of genes known to be involved in developmental disorders.")
    parser.add_argument("--de-novos", default=DENOVO_PATH,
        help="Path to table of variants from DDD.")
    parser.add_argument("--validations", default=VALIDATIONS,
        help="Path to table of results from DDD de novo validation experiments.")
    parser.add_argument("--diagnosed", default=DIAGNOSED,
        help="Path to table of results from DDD de novo validation experiments.")
    parser.add_argument("--external",
        help="Path to table of variants from external sequencing studies.")
    parser.add_argument("--output", default="diagnosed_per_subset.pdf",
        help="Path to plot output to.")
    
    args = parser.parse_args()
    
    return args

def get_ddd_diagnosed(de_novo_path, validations_path, diagnosed_path):
    """ get the DDD probands with a de novo, annotated with whether the have a
    likely diagnosis from a dominant de novo.
    
    Args:
        de_novo_path: path to table of candidate de novos identified in the DDD
            cohort.
        validations_path: path to table of validation results for candidate de
            novos identified in the DDD cohort.
        diagnosed_path: path to table of probands with diagnoses, or likely to
            have dominant de novo diagnoses in the DDD.
    
    Returns:
        pandas dataframe of unique probands with columns for person_id,
        study_phenotype, and diagnosed status.
    """
    
    # open the list of DDD probands with a diagnosis from a de novo mutation
    diagnosed = pandas.read_table(diagnosed_path)
    diagnosed = diagnosed[diagnosed["inheritance"] == "de_novo"]
    
    # determine which DDD probands have a diagnosis from a de novo
    internal = open_de_novos(de_novo_path, validations_path)
    internal["diagnosed"] = internal["person_stable_id"].isin(diagnosed["person_id"])
    
    # restrict to the unique set of probands
    internal = internal[~internal["person_stable_id"].duplicated()]
    
    # format the table to be able to join to the external probands
    internal["person_id"] = internal["person_stable_id"]
    internal["study_phenotype"] = "ddd_4k"
    internal = internal[["person_id", "study_phenotype", "diagnosed"]]
    
    return internal

def get_external_diagnosed(external_path, known_genes):
    """ get a table of external probands, along with whether they have a
    dominant diagnostic de novo.
    
    Args:
        external_path: path to table of de novos identified in external
            sequencing studies.
        known_genes: pandas dataframe of known developmental disorder genes,
            which includes gene symbols and the mode of inheritance expected for
            each gene.
    
    Returns:
        pandas dataframe of unique probands with columns for person_id,
        study_phenotype, and diagnosed status.
    """
    
    # identify the genes which might contribute to a dominant diagnosis
    dominant_modes = ["Monoallelic", "X-linked dominant"]
    symbols = known_genes["gene"]
    dominant = symbols[known_genes["mode"].isin(dominant_modes)]
    hemizygous = symbols[known_genes["mode"].isin(["Hemizygous"])]
    
    # load the de novos from the external cohorts, and identify which might have
    # a dominant diagnosis
    external = pandas.read_table(external_path, sep="\t", compression="gzip")
    diagnosed_ids = external["person_id"][external["hgnc"].isin(dominant) | \
        ((external["sex"] == "male") & external["hgnc"].isin(hemizygous))]
    external["diagnosed"] = external["person_id"].isin(diagnosed_ids)
    
    # restrict to the unique set of probands
    external = external[~external[["person_id", "study_code"]].duplicated()]
    
    external = external[["person_id", "study_phenotype", "diagnosed"]]
    
    return external

def plot_proportion_diagnosed(internal, external, plot_path):
    """ plot the proportion of probands with a dominant diagnostic de novo per subset
    
    Args:
        internal: pandas dataframe of probands in DDD, with columns for
            proband ID, whether the proband hasn a diagnostic dominant de
            novo, and a column naming the DDD study.
        external: pandas dataframe of probands in external studies, with columns
            for proband ID, whether the proband hasn a diagnostic dominant
            de novo, and a column for the phenotype studied in the cohort.
        plot_path: path to write a pdf of the plot to.
    """
    
    probands = internal.append(external, ignore_index=True)
    
    # count the number of diagnoised and undiagnosed probands for each
    # phenotypic subset
    counts = []
    for name, group in probands.groupby("study_phenotype"):
        diagnosed = group["person_id"][group["diagnosed"]].unique()
        not_diagnosed = group["person_id"][~group["diagnosed"]].unique()
        counts.append([name, len(diagnosed), len(not_diagnosed)])
    
    # format the counts as a dataframe, and determine the proportion of individuals with a
    counts = pandas.DataFrame(counts, columns=["phenotype", "diagnosed", "undiagnosed"])
    counts["total"] = counts["diagnosed"] + counts["undiagnosed"]
    counts["proportion"] = counts["diagnosed"]/counts["total"]
    
    counts = counts.sort("proportion", ascending=False)
    
    fig = seaborn.factorplot(x="phenotype", y="proportion", data=counts, \
        color="gray", kind="bar", size=6, aspect=1.2)
    fig = fig.set_xticklabels(rotation=90)
    ax = fig.set_ylabels("Proportion with dominant de novo")
    ax = fig.set_xlabels("Subset")
    ax = fig.ax.axhline(y=0.1, color="red", linestyle="dashed")
    fig.savefig(plot_path, format="pdf")

def main():
    
    args = get_options()
    
    known_genes = pandas.read_table(args.known_genes, sep="\t", index_col=False)
    
    internal = get_ddd_diagnosed(args.de_novos, args.validations, args.diagnosed)
    external = get_external_diagnosed(args.external, known_genes)
    
    plot_proportion_diagnosed(internal, external, args.output)

if __name__ == '__main__':
    main()
