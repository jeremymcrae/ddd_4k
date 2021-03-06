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
from scipy.stats import fisher_exact, mannwhitneyu, norm, mstats, pearsonr, spearmanr
from statsmodels.genmod.generalized_linear_model import GLM
import statsmodels.formula.api as smf
from statsmodels.genmod.families import Binomial
from matplotlib import pyplot

from mupit.open_ddd_data import open_known_genes
from ddd_4k.load_files import open_de_novos, open_phenotypes, open_families
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

def plot_categorical(counts, pheno, value, folder):
    """ plot probands per category by consequence by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (truncating/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        value: the phenotypic value that we are grouping by.
        folder: the folder to save graphs into.
    """
    
    ci_95 = norm.ppf(0.975)
    
    # identify the probands with a pathogenic de novo
    merged = pheno[["person_stable_id", value]].copy()
    probands = counts["person_stable_id"][counts["known"]]
    merged["has_causal"] = merged["person_stable_id"].isin(probands)
    
    table = merged.pivot_table(index="has_causal", columns=value, values="person_stable_id", aggfunc=len)
    odds_ratio, p_value = fisher_exact(table)
    
    standard_error = math.sqrt((1/table).sum().sum())
    log_odds = math.log(odds_ratio)
    upper_ci = math.exp(log_odds + ci_95 * standard_error)
    lower_ci = math.exp(log_odds - ci_95 * standard_error)
    
    ratios = {}
    ratios["name"] = value
    ratios["odds_ratio"] = odds_ratio
    ratios["upper"] = upper_ci
    ratios["lower"] = lower_ci
    ratios["p_value"] = p_value
    
    print(ratios)
    
    return ratios

def plot_quantitative(counts, pheno, value, folder, y_label, delta_to_median=False, covariate=None):
    """ plot quantitative metric by functional category by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (truncating/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        value: the phenotypic value that we are grouping by.
        folder: the folder to save graphs into.
        y_label: the y-axis label for the phenotype.
        delta_to_median: whether to adjust the data, so we test for the abolute
            difference from the median value, which could be useful for
            phenotypes where we expect genetic variants to spread affected
            probands to the extremes of the distributions.
        covariate: name of addiditional covariate to adjust for in the logistic
            regression.
    """
    
    columns = ["person_stable_id", "gender", value]
    if covariate is not None:
        columns.append(covariate)
    
    data = pheno[columns].copy()
    probands = counts["person_stable_id"][counts["known"]]
    data["has_causal"] = data["person_stable_id"].isin(probands)
    data["has_causal"] = data["has_causal"].map({True: 1, False: 0})
    
    if delta_to_median:
        median = numpy.median(data[value])
        data[value] = abs(data[value] - median)
    
    ratios = {}
    data = data.dropna()
    data["z_score"] = mstats.zscore(data[value])
    
    formula = 'has_causal ~ z_score*gender'
    if covariate is not None:
        formula += ' + {}:gender'.format(covariate)
        data[covariate] = mstats.zscore(data[covariate])
    
    data["gender"] = data["gender"].map({"Female": False, "Male": True})
    model = smf.glm(formula=formula, data=data, family=Binomial())
    result = model.fit()
    
    ratios["name"] = value
    ratios["beta"] = result.params["z_score"]
    ratios["upper"] = result.conf_int()[1]["z_score"]
    ratios["lower"] = result.conf_int()[0]["z_score"]
    ratios["p_value"] = result.pvalues["z_score"]
    
    print(ratios)
    
    fig = seaborn.factorplot(x="has_causal", y=value, size=6, data=data,
        aspect=0.6, kind="violin")
    fig.set_ylabels(y_label)
    
    fig.savefig("{}/{}_by_consequence.pdf".format(folder, value), format="pdf",
        bbox_inches='tight', pad_inches=0, transparent=True)
    
    matplotlib.pyplot.close()
    
    return ratios

def plot_achievement(counts, pheno, achievement, folder):
    """ plot developmental milestone achievement ages by consequence by known gene status
    
    Args:
        counts: dataframe of number of de novos per proband per consequence type
            (truncating/functional) by known developmental gene status.
        pheno: dataframe of phenotypic values for probands
        achievement: column name for the developmental milestone e.g.
            "social_smile", "first_words".
    """
    
    
    cap_ages = {'social_smile': 0.25, 'sat_independently': 0.75,
        'walked_independently': 1.5, 'first_words': 2.0}
    cap_age = cap_ages[achievement]
    
    # Convert the achievement age to number of seconds since birth (rather than
    # having values like 5 weeks, 6 months etc), then log10 transform the
    # duration, so that the values are more normally distributed.
    durations, unit = autoscale_durations(pheno[[achievement, 'decimal_age']].apply(get_duration, axis=1, cap_age=cap_age))
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
    data = data.reindex(index=data.index[::-1])
    data["name"] = data["name"].str.replace("_", " ")
    
    data.index = range(len(data))
    data["y"] = data.index
    
    # format the p-values to 3 significant figures
    data["p_value"] = [ "P = {0:.3g}".format(x) for x in data["p_value"] ]
    
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
    
    # and include the p-values on the opposing y-axis
    new_ax = ax.twinx()
    new_ax.set_ylim((min(data["y"]) - 0.5), max(data["y"]) + 0.5)
    e = new_ax.spines['right'].set_visible(False)
    e = new_ax.spines['top'].set_visible(False)
    e = new_ax.spines['left'].set_visible(False)
    
    e = new_ax.yaxis.set_ticks_position('right')
    e = new_ax.xaxis.set_ticks_position('bottom')
    
    e = new_ax.set_yticks(data["y"])
    e = new_ax.set_yticklabels(data["p_value"])
    
    fig.savefig("{}.pdf".format(value), format="pdf", bbox_inches='tight',
        pad_inches=0, transparent=True)

def main():
    args = get_options()
    
    de_novos = open_de_novos(args.de_novos, args.validations)
    known = open_known_genes(args.ddg2p)
    monoallelic = set(known["gene"][known["mode"].isin(["Monoallelic"])])
    de_novos["known"] = de_novos["hgnc"].isin(monoallelic)
    
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
    
    recode = {'NA': False, 'Neither': False, 'Unknown': False, '': False, 'None': False, None: False, 'No': False,
        'Both': True, 'Mother': True, 'Father': True, '1': True, '2': True,
        '3 or more': True, 'Yes': True}
    pheno['similar_phenotype_parents'] = pheno['similar_phenotype_parents'].map(recode)
    pheno['similar_phenotype_siblings'] = pheno['similar_phenotype_siblings'].map(recode)
    pheno['similar_phenotype_relatives'] = pheno['similar_phenotype_relatives'].map(recode)
    
    recode = {'Singleton': 1, 'Twin': 2, 'Triplet': 3, '4+': 4}
    pheno['multiple_births'] = pheno['multiple_births'].map(recode)
    
    recode = {'None': 0, '1': 1, '2': 2, '3 or more': 3, '3 or More': 3}
    pheno['pregnancy_loss_history'] = pheno['pregnancy_loss_history'].map(recode)
    
    ratios = []
    ratios.append(plot_categorical(counts, pheno, "gender", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "feeding_problems", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "anticonvulsant_drugs", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "maternal_illness", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "bleeding", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "abnormal_scan", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "scbu_nicu", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "assisted_reproduction", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "similar_phenotype_parents", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "similar_phenotype_siblings", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "similar_phenotype_relatives", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "consanguinity", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "cranial_mri_scan_outcome", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "maternal_diabetes", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "only_patient_affected", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "x_linked_inheritance", args.output_folder))
    ratios.append(plot_categorical(counts, pheno, "increased_nt", args.output_folder))
    
    odds_ratios = pandas.DataFrame(ratios)
    forest_plot(odds_ratios)
    
    betas = []
    # plot distributions of time to achieve developmental milestones by
    # the functional categories
    betas.append(plot_achievement(counts, pheno, "first_words", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "walked_independently", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "sat_independently", args.output_folder))
    betas.append(plot_achievement(counts, pheno, "social_smile", args.output_folder))
    
    # birthweight (or birthweight corrected for duration of gestation, or
    # birthweight_percentile (if that corrects for duration of gestation))
    betas.append(plot_quantitative(counts, pheno, "child_hpo_n", args.output_folder, "Number of proband HPO terms"))
    betas.append(plot_quantitative(counts, pheno, "height_sd", args.output_folder, "Height (SD)", delta_to_median=True))
    betas.append(plot_quantitative(counts, pheno, "birthweight_sd", args.output_folder, "Birthweight (SD)", delta_to_median=True))
    betas.append(plot_quantitative(counts, pheno, "ofc_sd", args.output_folder, "OFC (SD)", delta_to_median=True))
    betas.append(plot_quantitative(counts, pheno, "weight_sd", args.output_folder, "weight (SD)"))
    betas.append(plot_quantitative(counts, pheno, "decimal_age_at_assessment", args.output_folder, "Age at assessment (years)"))
    betas.append(plot_quantitative(counts, pheno, "gestation", args.output_folder, "Gestation duration (weeks)"))
    betas.append(plot_quantitative(counts, pheno, "fathers_age", args.output_folder, "Father's age (years)"))
    betas.append(plot_quantitative(counts, pheno, "mothers_age", args.output_folder, "Mother's age (years)"))
    betas.append(plot_quantitative(counts, pheno, "multiple_births", args.output_folder, "Multiple births"))
    betas.append(plot_quantitative(counts, pheno, "pregnancy_loss_history", args.output_folder, "History of preganacy loss"))
    
    # and add in the values from the logistic regression of autozygosity length
    # vs having a dominant diagnostic de novo. These values are determined in
    # the autozygosity_vs_diagnosed.py script
    betas.append({'upper': -0.1154436448891767, 'beta': -0.1845241470915473, 'lower': -0.25360464929391791, 'name': 'autozygosity_length', 'p_value': 1.646715410585663e-07})
    
    betas = pandas.DataFrame(betas)
    forest_plot(betas)
    

if __name__ == '__main__':
    main()
