# Copyright (c) 2015 Genome Research Ltd.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import random

from scipy.stats import t
from numpy import sqrt, mean, std
import pandas
import matplotlib
matplotlib.use("Agg")
import seaborn

from mupit.open_ddd_data import standardise_ddd_de_novos, get_ddd_rates, \
    open_known_genes
from mupit.gene_enrichment import analyse_enrichment

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """
    """
    
    parser = argparse.ArgumentParser(
        description="test the power to detect dominant genes with exome and genome sequencing")
    parser.add_argument("--rates", help="Path to table of mutation rates.",
        default="/nfs/users/nfs_j/jm33/apps/denovonear/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt")
    parser.add_argument("--de-novos", help="Path to DDD de novo dataset.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.txt")
    parser.add_argument("--validations", help="Path to validation results.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-10-12.txt")
    parser.add_argument("--families", help="Path to families PED file.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
    parser.add_argument("--trios", help="Path to file listing complete trios.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt")
    parser.add_argument("--known", help="Path to table of known developmental disorder genes.",
        default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
    parser.add_argument("--iterations", type=int,
        help="number of iterations to run per condition.", default=100)
    parser.add_argument("--output", help="Path to plot pdf to file.",
        default="exome_vs_genome.pdf")
    
    args = parser.parse_args()
    
    return args


#' get the probands from trios in the current DDD dataset
#'
#' @param families_path path to family relationships ped file, containing sample
#'        IDs and sex
#' @param trios_path path to table of information on the trios only
#'
#' @return dataframe of proband IDs and sex for the trio-based probands
def get_probands(families_path, trios_path):
    
    families = pandas.read_table(families_path, sep="\t")
    
    # determine the trios with exome data available
    trios = pandas.read_table(trios_path, sep="\t")
    probands = families[families["individual_id"].isin(trios["proband_stable_id"])]
    
    probands = probands[["individual_id", "sex"]].copy()
    
    return probands

#' get the number of male and female probands in a cohort, or cohort subset
#'
#' @param probands path to table of DDD de novos
#'
#' @return list, with entries for counts of trios with male and female probands
def get_trio_counts(probands):
    # get the number of trios studied in our data for each sex
    sex = probands["sex"].value_counts()
    male = sex[["M"]]
    female = sex[["F"]]
    
    return {"male": int(male), "female": int(female)}

#' open de novo mutations from the DDD cohort, removing variants that failed validation
#'
#' @param de_novos_path path to table of DDD de novos
#' @param validations_path path to table of results from validation experiments
#'
#' @return dataframe of de novos in the DDD dataset
def get_de_novos(de_novos_path, validations_path):
    
    variants = standardise_ddd_de_novos(de_novos_path)
    
    validations = pandas.read_table(validations_path, sep="\t")
    variants = variants.merge(validations, how="left", on=["person_id", "chrom",
        "start_pos", "end_pos", "ref_allele", "alt_allele", "hgnc",
        "consequence"])
    
    # drop out the variants that failed to validate (i.e. were false positives,
    # or inherited)
    variants = variants[~variants["status"].isin(["false_positive", "inherited"])]
    del variants["status"]
    
    return variants

#' get the dominant DDG2P genes reaching genomewide significance
#'
#' @param n_trios number of trios to sample from among the probands
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param genes vector of dominant DDG2P genes
#'
#' @return vector of HGNC symbols for DDG2P genes reaching genomewide significance
def get_enrichment_in_sample(n_trios, probands, de_novos, rates, threshold, genes):
    
    # Get a sampled subset of probands, then figure out the male and female
    # sex counts in the subset. Restrict the de novos to ones from probands
    # in the subset.
    try:
        sampled_ids = random.sample(probands["individual_id"], int(n_trios))
    except ValueError:
        # if we try to sample more trios than exists in the current DDD cohort,
        # then return NA, rather than using the full cohort. This will drop out
        # these data points
        return None
    
    probands = probands[probands["individual_id"].isin(sampled_ids)]
    trios = get_trio_counts(probands)
    de_novos = de_novos[de_novos["person_id"].isin(sampled_ids)]
    
    # figure out the enrichment in the
    enrich = analyse_enrichment(de_novos, trios, rates=rates)
    
    enrich["p_min"] = enrich[["p_lof", "p_func"]].min(axis=1, skipna=True)
    enrich["genomewide"] = enrich["p_min"] < threshold
    enrich["dominant"] = enrich["hgnc"].isin(genes)
    
    dominant = enrich["hgnc"][enrich["genomewide"] & enrich["dominant"]]
    dominant = dominant[~dominant.isnull()]
    
    return dominant

#' detect genes reaching genomewide significance, under specific conditions
#'
#' Given a budget, relative cost of exome sequencing to genome sequencing, and
#' relative sensitivity of genome sequencing to detect de novo mutations,
#' determine the number of genes that could be detected by sequencing as many
#' trios as the budget allows. This function separately checks many random
#' samples of probands used for each test, so as to capture the variation in the
#' population.
#'
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param dominant vector of dominant DDG2P genes
#' @param genomewide vector of genes reaching genomewide significance in
#'        in the complete DDD dataset
#' @param genome_cost reference cost of performing genome sequencing
#' @param budget budget for the current conditions
#' @param relative_cost relative cost of exome sequencing versus genome sequencing
#' @param sensitivity relative sensitivity of genome sequencing for detecting de
#'        novo mutations compared to exome sequencing.
#' @param iterations number of iterations to run for the current conditions
#'
#' @return dataframe of test conditions, numbers of genes reaching genomewide
#'         significance from genome and exome sequencing
def run_iterations (probands, de_novos, rates, threshold, dominant,
    genomewide, genome_cost, budget, relative_cost, sensitivity,
    iterations):
    
    # Determine the number of trios that could be sequenced, given the budget.
    # Adjust the number of genome trios upwards by the relative sensitivity of
    # genome sequencing to exome sequencing, to account for the additional de
    # novo mutations that would be detected by genome sequencing.
    n_genome_trios = budget/(genome_cost * 3) * sensitivity
    n_exome_trios = budget/((genome_cost * relative_cost) * 3)
    
    power = pandas.DataFrame(columns=["budget", "relative_cost",
        "sensitivity", "genome", "exome"])
    
    for n in range(iterations):
        # We only need to test the exome trios at a sensitivity of 1.0, since
        # variation in sensitity only affects the genome results.
        n_exome = None
        if sensitivity == 1.0:
            exome_genes = get_enrichment_in_sample(n_exome_trios, probands,
                de_novos, rates, threshold, dominant)
            if exome_genes is not None:
                n_exome = sum(exome_genes.isin(genomewide))
        
        n_genome = None
        if relative_cost == 1.0:
            genome_genes = get_enrichment_in_sample(n_genome_trios, probands,
                de_novos, rates, threshold, dominant)
            if genome_genes is not None:
                n_genome = sum(genome_genes.isin(genomewide))
        
        temp = pandas.DataFrame({"budget": [budget],
            "relative_cost": [relative_cost], "sensitivity": [sensitivity],
            "genome": [n_genome], "exome": [n_exome]})
        
        power = power.append(temp, ignore_index=True)
    
    return power

#' simulate power of exome and genome sequencing to detect dominant genes, under various conditions
#'
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param dominant vector of dominant DDG2P genes
#' @param genomewide vector of genes reaching genomewide significance in
#'        in the complete DDD dataset
#' @param iterations number of iterations to run for the current conditions
#'
#' @return dataframe of power simulations for each condition, containing numbers
#'         of genes reaching genomewide significance from genome and exome sequencing
def simulate_power(probands, de_novos, rates, threshold, dominant, genomewide, iterations):
    budgets = [1e6, 2e6, 5e6]
    genome_cost = 1000
    exome_relative_cost = [ (x+1)/5.0 for x in range(5) ]
    genome_sensitivity = [ (x/20.0)+1 for x in range(5) ]
    
    for budget in budgets:
        for relative_cost in exome_relative_cost:
            for sensitivity in genome_sensitivity:
                temp = run_iterations(probands, de_novos, rates, threshold,
                    dominant, genomewide, genome_cost, budget,
                    relative_cost, sensitivity, iterations)
                
                if "power" not in locals():
                    power = temp
                else:
                    power = power.append(temp, ignore_index=True)
    
    return power

#' reshape the power dataframe, so that we have the mean number of genes
#' reaching genomewide significance for each condition, along with the
#' confidence intervals.
#'
#' @param power dataframe of numbers for genome and exome sequencing
#' @param conf_interval confidence interval
#'
#' @return
def get_mean_and_ci(power, conf_interval=0.95):
    
    power = pandas.melt(power, id_vars=["budget", "relative_cost", "sensitivity"])
    power = power[~power["value"].isnull()]
    power["value"] = power["value"].astype(int)
    
    reshaped = pandas.pivot_table(power, rows="variable",
        cols=["budget", "relative_cost", "sensitivity"],
        values="value", aggfunc=[mean, std, len])
    
    means = reshaped["mean"].transpose()
    stdevs = reshaped["std"].transpose()
    lengths = reshaped["len"].transpose()
    
    # find the t-statistic given the number of points for each group
    genome_ci_mult = t.interval(conf_interval, lengths["genome"])[1]
    exome_ci_mult = t.interval(conf_interval, lengths["exome"])[1]
    
    # find the confidence intervals
    conf = means.copy()
    conf["genome"] = stdevs["genome"]/sqrt(lengths["genome"]) * genome_ci_mult
    conf["exome"] = stdevs["exome"]/sqrt(lengths["exome"]) * exome_ci_mult
    
    index = means.index
    for x, name in enumerate(index.names):
        levels = index.levels[x]
        labels = index.labels[x]
        means[name] = [ levels[x] for x in labels ]
        conf[name] = [ levels[x] for x in labels ]
    
    # reshape the dataset so we have the mean and confidence intervals in a
    # single dataframe
    means = pandas.melt(means, id_vars=["budget", "relative_cost", "sensitivity"])
    conf = pandas.melt(conf, id_vars=["budget", "relative_cost", "sensitivity"])
    means["ci"] = conf["value"]
    means = means[~means["value"].isnull()]
    
    return means

#' plot the results from simulation power of exome and genome sequencing
#'
#' @param power dataframe of power simulations for each condition, containing
#'        numbers of genes reaching genomewide significance from genome and
#'        exome sequencing
#' @param output_path path to save output plot to
def plot_power(power, output_path):
    power = get_mean_and_ci(power)
    
    power["sequence"] = power["variable"] + "-" + power["sensitivity"].astype(str)
    
    fig = seaborn.factorplot(x="relative_cost", y="value", hue="sequence",
        col="budget", data=power, size=8, aspect=0.4)
    
    for ax in fig.axes[0]:
        e = ax.set_ylim((0, max(power["value"])))
    
    fig.savefig(output_path, format="pdf")

def main():
    args = get_options()
    
    rates = get_ddd_rates(args.rates)
    known = open_known_genes(args.known)
    dominant = sorted(set(known["gene"][known["mode"].isin(["Monoallelic", "X-linked dominant"])]))
    
    # set the multiple testing corrected threshold for genomewide significance
    alpha = 0.01
    num_genes = 18500
    num_tests = num_genes * 2
    threshold = alpha/num_tests
    
    probands = get_probands(args.families, args.trios)
    de_novos = get_de_novos(args.de_novos, args.validations)
    genomewide = get_enrichment_in_sample(len(probands), probands,
        de_novos, rates, threshold, dominant)
    
    power = simulate_power(probands, de_novos, rates, threshold, dominant,
        genomewide, args.iterations)
    
    plot_power(power, args.output)

if __name__ == '__main__':
    main()
