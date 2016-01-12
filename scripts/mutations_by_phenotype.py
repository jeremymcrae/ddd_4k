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

import os
import argparse

import numpy

import matplotlib
matplotlib.use('Agg')
import seaborn
import pandas
from scipy.stats import fisher_exact, mannwhitneyu
from matplotlib import pyplot

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes, \
    open_families
from ddd_4k.count_hpo import count_hpo_terms
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, PHENOTYPES, SANGER_IDS, \
    FAMILIES, TRIOS, VALIDATIONS
from ddd_4k.count_mutations_per_person import get_count_by_person
from ddd_4k.convert_durations import get_duration
from ddd_4k.scale_durations import autoscale_durations

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check whether"
        "probands with longer autozygous regions are less likely to have"
        "diagnoses")
    parser.add_argument("--de-novos", default=DENOVO_PATH, \
        help="Path to file of canddidate de novo variants.")
    parser.add_argument("--ddg2p", default=KNOWN_GENES, \
        help="path to known genes (DDG2P).")
    parser.add_argument("--phenotypes", default=PHENOTYPES, \
        help="Path to table of phenotypes for all probands.")
    parser.add_argument("--sanger-ids", default=SANGER_IDS, \
        help="Path to file mapping between Decipher IDs and DDD IDs.")
    parser.add_argument("--families", default=FAMILIES, \
        help="Path to listing family relationships.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Path to file listing all exome-sequenced trios.")
    parser.add_argument("--validations", default=VALIDATIONS, \
        help="Path to file listing validation results.")
    parser.add_argument("--output-folder", default="results", \
        help="Path to plot graph to.")
    
    args = parser.parse_args()
    
    return args

def get_categorical_summary(merged, cohort_counts, value):
    """ count the number of individuals by known gene status and functional
    de novo status. Include the totoal number of individuals in each group, so
    we can estimate ratios, and perform Fisher's exact test to compare groups.
    
    Args:
        merged: pandas DataFrame listing the indivioduals with de novos by
            functional status and known  gene status, along with their
            phenotypic data, so we can group by the required "value".
        cohort_counts: pandas Series listing the total number of individuals in
            the cohort for the different groups split by the "value".
        value: the phenotypic value that we are grouping by.
        folder: the folder to save graphs into.
    
    Returns:
        pandas DataFrame containing tallies of individuals by functional de novo
        status and known gene status, for the different "value" groups.
    """
    
    summary = pandas.DataFrame({"known": [], value: [], "consequence": [], \
        "tally": [], "total": []})
    
    for known in [True, False]:
        for x in cohort_counts.keys():
            temp = merged[(merged.known == known) & (merged[value] == x)]
            tally = temp["consequence"].value_counts()
            total = [cohort_counts[x]] * len(tally)
            cq = list(tally.keys())
            
            temp = pandas.DataFrame({"known": [known] * len(tally), \
                value: [x] * len(tally), "consequence": cq, "tally": tally, \
                "total": total})
            summary = summary.append(temp)
    
    summary = summary.reset_index(drop=True)
    summary["ratio"] = summary["tally"]/summary["total"]
    
    return summary

def plot_categorical(counts, pheno, value, folder):
    """ plot probands per category by consequence by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (loss-of-function/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        value: the phenotypic value that we are grouping by.
        folder: the folder to save graphs into.
    """
    
    merged = counts.merge(pheno[["person_stable_id", value]], on="person_stable_id")
    
    cohort_counts = pheno[value].value_counts()
    summary = get_categorical_summary(merged, cohort_counts, value)
    
    # test whether differences exist between the groups for the value column
    results = []
    for known in [True, False]:
        for cq in summary.consequence.unique():
            table = summary[(summary.known == known) & (summary.consequence == cq)]
            table = table[["tally", "total"]]
            odds_ratio, p_value = fisher_exact(table)
            results.append([known, cq, p_value])
    
    fig = seaborn.factorplot(x="consequence", y="ratio", hue=value, col="known", data=summary, kind="bar")
    fig.set_ylabels("Frequency")
    pyplot.table(cellText=results, colLabels=["known", "cq", "P"], loc="top right")
    fig.savefig("{}/{}_by_consequence.pdf".format(folder, value), format="pdf")
    
    matplotlib.pyplot.close()

def plot_quantitative(counts, pheno, value, folder, y_label):
    """ plot quantitative metric by functional category by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (loss-of-function/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        value: the phenotypic value that we are grouping by.
        folder: the folder to save graphs into.
        y_label: the y-axis label for the phenotype.
    """
    
    merged = counts.merge(pheno[["person_stable_id", value]], on="person_stable_id")
    
    results = []
    for known in [True, False]:
        lof = merged[value][(merged.known == known) & (merged.consequence == "loss-of-function")]
        func = merged[value][(merged.known == known) & (merged.consequence == "functional")]
        u, p_value = mannwhitneyu(lof, func)
        results.append([known, p_value])
    
    fig = seaborn.factorplot(x="known", y=value, hue="consequence", size=6, data=merged, kind="violin")
    fig.set_ylabels(y_label)
    pyplot.table(cellText=results, colLabels=["known", "cq", "P"], loc="top right")
    fig.savefig("{}/{}_by_consequence.pdf".format(folder, value), format="pdf")
    
    matplotlib.pyplot.close()

def plot_hpo_by_consequence(counts, pheno, folder):
    """ Plot number of HPO terms by functional category by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (loss-of-function/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        value: the phenotypic value that we are grouping by.
    """
    
    pheno["child_hpo_n"] = count_hpo_terms(pheno, "child")
    hpo_counts = counts.merge(pheno[["person_stable_id", "child_hpo_n"]], on="person_stable_id")
    
    fig = seaborn.factorplot(x="known", y="child_hpo_n", hue="consequence", size=6, data=hpo_counts, kind="violin")
    fig.set_ylabels("HPO terms per proband (n)")
    fig.savefig("{}/hpo_by_consequence.pdf".format(folder), format="pdf")
    matplotlib.pyplot.close()

def plot_achievement(counts, pheno, achievement, folder):
    """ plot developmental milestone achievement ages by consequence by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (loss-of-function/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        achievement: column name for the developmental milestone e.g.
            "social_smile", "first_words".
    """
    
    # Convert the achievement age to number of seconds since birth (rather than
    # having values like 5 weeks, 6 months etc), then log10 transform the
    # duration, so that the values are more normally distributed.
    durations, unit = autoscale_durations(pheno[achievement].apply(get_duration))
    pheno[achievement] = numpy.log10(durations)
    plot_quantitative(counts, pheno, achievement, folder, "{} log10({}s)".format(achievement, unit))

def main():
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations)
    known = open_known_genes(args.ddg2p)
    monoallelic = known[known["mode"].isin(["Monoallelic", "X-linked dominant"])]
    de_novos["known"] = de_novos["symbol"].isin(monoallelic["gencode_gene_name"])
    
    # For each proband, count the number of functional de novos (split by lof
    # and missense), in known developmental disorder genes (and other genes).
    counts = get_count_by_person(de_novos)
    
    # load the phenotype data
    pheno = open_phenotypes(args.phenotypes, args.sanger_ids)
    
    # restrict the phenotype dataset to the probands for whom we could have de novo
    # candidates, that is the set of probands where all the members of their trio
    # have exome sequence data available.
    families = open_families(args.families, args.trios)
    probands = families["individual_id"][~families["mother_stable_id"].isnull()]
    pheno = pheno[pheno["person_stable_id"].isin(probands)]
    
    plot_hpo_by_consequence(counts, pheno, args.output_folder)
    plot_categorical(counts, pheno, "gender", args.output_folder)
    plot_categorical(counts, pheno, "scbu_nicu", args.output_folder)
    plot_categorical(counts, pheno, "feeding_problems", args.output_folder)
    plot_categorical(counts, pheno, "maternal_illness", args.output_folder)
    plot_categorical(counts, pheno, "bleeding", args.output_folder)
    plot_categorical(counts, pheno, "abnormal_scan", args.output_folder)
    plot_categorical(counts, pheno, "assisted_reproduction", args.output_folder)
    
    # birthweight (or birthweight corrected for duration of gestation, or
    # birthweight_percentile (if that corrects for duration of gestation))
    plot_quantitative(counts, pheno, "decimal_age_at_assessment", args.output_folder, "Age at assessment (years)")
    plot_quantitative(counts, pheno, "birthweight", args.output_folder, "Birthweight (grams)")
    plot_quantitative(counts, pheno, "gestation", args.output_folder, "Gestation duration (weeks)")
    plot_quantitative(counts, pheno, "height_percentile", args.output_folder, "Height (percentile)")
    plot_quantitative(counts, pheno, "weight_percentile", args.output_folder, "weight (percentile)")
    plot_quantitative(counts, pheno, "ofc_percentile", args.output_folder, "OFC (percentile)")
    
    # plot distributions of time to achieve developmental milestones by
    # the functional categories
    plot_achievement(counts, pheno, "social_smile", args.output_folder)
    plot_achievement(counts, pheno, "sat_independently", args.output_folder)
    plot_achievement(counts, pheno, "walked_independently", args.output_folder)
    plot_achievement(counts, pheno, "first_words", args.output_folder)
    

if __name__ == '__main__':
    main()
