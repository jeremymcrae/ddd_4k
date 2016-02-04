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
import math

import numpy

import matplotlib
matplotlib.use('Agg')
import seaborn
import pandas
from scipy.stats import fisher_exact, mannwhitneyu, norm, mstats
from statsmodels.genmod.generalized_linear_model import GLM
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Binomial
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
    
    ci_95 = norm.ppf(0.975)
    
    merged = counts.merge(pheno[["person_stable_id", value]], on="person_stable_id")
    
    cohort_counts = pheno[value].value_counts()
    summary = get_categorical_summary(merged, cohort_counts, value)
    
    # test whether differences exist between the groups for the value column
    ratios = {}
    results = []
    for known in [True, False]:
        for cq in summary.consequence.unique():
            table = summary[(summary.known == known) & (summary.consequence == cq)]
            table = table[["tally", "total"]]
            odds_ratio, p_value = fisher_exact(table)
            
            standard_error = math.sqrt((1/table).sum().sum())
            log_odds = math.log(odds_ratio)
            upper_ci = math.exp(log_odds + ci_95 * standard_error)
            lower_ci = math.exp(log_odds - ci_95 * standard_error)
            results.append([known, cq, p_value])
            
            if known == True and cq == "loss-of-function":
                ratios["name"] = value
                ratios["odds_ratio"] = odds_ratio
                ratios["upper"] = upper_ci
                ratios["lower"] = lower_ci
    
    fig = seaborn.factorplot(x="consequence", y="ratio", hue=value, col="known", data=summary, kind="bar")
    fig.set_ylabels("Frequency")
    pyplot.table(cellText=results, colLabels=["known", "cq", "P"], loc="top right")
    fig.savefig("{}/{}_by_consequence.pdf".format(folder, value), format="pdf")
    
    matplotlib.pyplot.close()
    
    return ratios

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
    
    merged = counts.merge(pheno[["person_stable_id", "gender", value]], on="person_stable_id")
    
    ratios = {}
    results = []
    for known in [True, False]:
        lof = merged[value][(merged.known == known) & (merged.consequence == "loss-of-function")]
        func = merged[value][(merged.known == known) & (merged.consequence == "functional")]
        u, p_value = mannwhitneyu(lof, func)
        results.append([known, p_value])
        
        if known:
            data = merged[[value, "gender", "consequence"]][(merged["known"] == known) &
                (merged["consequence"].isin(["functional", "loss-of-function"]))].copy()
            recode = {"functional": 0, "loss-of-function": 1}
            data["consequence"] = data["consequence"].map(recode)
            
            data = data.dropna()
            data[value] = mstats.zscore(data[value])
            
            data["gender"] = data["gender"].map({"Female": 1, "Male": 0})
            model = smf.glm(formula='consequence ~ {}:gender'.format(value),
                data=data, family=Binomial())
            result = model.fit()
            
            ratios["name"] = value
            ratios["beta"] = result.params[1]
            ratios["upper"] = result.conf_int()[1][1]
            ratios["lower"] = result.conf_int()[0][1]
    
    fig = seaborn.factorplot(x="known", y=value, hue="consequence", size=6, data=merged, kind="violin")
    fig.set_ylabels(y_label)
    pyplot.table(cellText=results, colLabels=["known", "cq", "P"], loc="top right")
    fig.savefig("{}/{}_by_consequence.pdf".format(folder, value), format="pdf")
    
    matplotlib.pyplot.close()
    
    return ratios

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
    ratios = plot_quantitative(counts, pheno, achievement, folder, "{} log10({}s)".format(achievement, unit))
    
    return ratios

def forest_plot(data):
    """make a forest-like plot of 95% CI of odds ratios or betas.
    """
    
    if "beta" in data.columns:
        value = "beta"
        h_pos = 0
    elif "odds_ratio" in data.columns:
        value = "odds_ratio"
        h_pos = 1
    
    fig = pyplot.figure(figsize=(6, len(data) * 0.5))
    ax = fig.gca()
    
    # sort by the value, to group variables by effects
    data = data.sort(value)
    data["name"] = data["name"].str.replace("_", " ")
    
    data.index = range(len(data))
    data["y"] = data.index
    
    ax.plot(data[value], data["y"], linestyle='None', marker="o", color="gray")
    for key, x in data.iterrows():
        ax.plot([x["lower"], x["upper"]], [x["y"]] * 2, color="gray")
    
    # add a vertical line which the CI should intersect
    ax.axvline(h_pos, color="black", linestyle="dashed")
    ax.set_ylim((min(data["y"]) - 0.5), max(data["y"]) + 0.5)
    
    # Hide the right, top and bottom spines
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['left'].set_visible(False)
    e = ax.tick_params(direction='out')
    
    # Only show ticks on the left and bottom spines
    e = ax.yaxis.set_ticks_position('left')
    e = ax.xaxis.set_ticks_position('bottom')
    
    e = ax.set_yticks(data["y"])
    e = ax.set_yticklabels(data["name"])
    e = ax.set_xlabel(value)
    
    fig.savefig("{}.pdf".format(value), format="pdf")

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
    pheno["child_hpo_n"] = count_hpo_terms(pheno, "child")
    
    ratios = []
    ratios.append(plot_categorical(counts, pheno, "gender", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "scbu_nicu", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "feeding_problems", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "maternal_illness", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "bleeding", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "abnormal_scan", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "assisted_reproduction", args.output_folder))
    
    odds_ratios = pandas.DataFrame(ratios)
    forest_plot(odds_ratios)
    
    # birthweight (or birthweight corrected for duration of gestation, or
    # birthweight_percentile (if that corrects for duration of gestation))
    betas = []
    betas.append(plot_quantitative(counts, pheno, "child_hpo_n", args.output_folder, "Number of proband HPO terms"))
    betas.append(plot_quantitative(counts, pheno, "decimal_age_at_assessment", args.output_folder, "Age at assessment (years)"))
    betas.append(plot_quantitative(counts, pheno, "birthweight_sd", args.output_folder, "Birthweight (SD)"))
    betas.append(plot_quantitative(counts, pheno, "gestation", args.output_folder, "Gestation duration (weeks)"))
    betas.append(plot_quantitative(counts, pheno, "height_sd", args.output_folder, "Height (SD)"))
    # betas.append(plot_quantitative(counts, pheno, "weight_sd", args.output_folder, "weight (SD)"))
    betas.append(plot_quantitative(counts, pheno, "ofc_sd", args.output_folder, "OFC (SD)"))
    betas.append(plot_quantitative(counts, pheno, "fathers_age", args.output_folder, "Father's age (years)"))
    
    # and add in the values from the logistic regression of autozygosity length
    # vs having a dominant diagnostic de novo. These values are determined in
    # the autozygosity_vs_diagnosed.py script
    betas.append({'upper': -0.13772607320296332, 'beta': -0.25055052100431563, 'lower': -0.36337496880566794, 'name': 'autozygosity_length'})
    
    # plot distributions of time to achieve developmental milestones by
    # the functional categories
    betas.append(plot_achievement(counts, pheno, "social_smile", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "sat_independently", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "walked_independently", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "first_words", args.output_folder))
    
    betas = pandas.DataFrame(betas)
    forest_plot(betas)
    

if __name__ == '__main__':
    main()
