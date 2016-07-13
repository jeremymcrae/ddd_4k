"""
Copyright (c) 2016 Genome Research Ltd.

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
    parser.add_argument("--seizures", help="Path to table of probands with seizures.",
        default="/nfs/users/nfs_j/jm33/seizure_probands.txt")
    parser.add_argument("--iterations", type=int,
        help="number of iterations to run per condition.", default=100)
    parser.add_argument("--output", help="Path to plot pdf to file.",
        default="seizure_vs_complete_power.pdf")
    
    args = parser.parse_args()
    
    return args

def get_probands(families_path, trios_path):
    """ get the probands from trios in the current DDD dataset
    
    Args:
        families_path: path to family relationships ped file, containing sample
            IDs and sex
        trios_path: path to table of information on the trios only
    
    Returns:
        dataframe of proband IDs and sex for the trio-based probands
    """
    
    families = pandas.read_table(families_path, sep="\t")
    
    # determine the trios with exome data available
    trios = pandas.read_table(trios_path, sep="\t")
    probands = families[families["individual_id"].isin(trios["proband_stable_id"])]
    
    probands = probands[["individual_id", "sex"]].copy()
    
    return probands

def get_trio_counts(probands):
    """ get the number of male and female probands in a cohort, or cohort subset
    
    Args:
        probands: path to table of DDD de novos
    
    Returns:
        dictionary of counts of trios with male and female probands
    """
    
    # get the number of trios studied in our data for each sex
    sex = probands["sex"].value_counts()
    male = sex[["M"]]
    female = sex[["F"]]
    
    return {"male": int(male), "female": int(female)}

def get_de_novos(de_novos_path, validations_path):
    """ open de novo mutations from the DDD cohort, removing variants that failed validation
    
    Args:
        de_novos_path: path to table of DDD de novos
        validations_path: path to table of results from validation experiments
    
    Returns:
        dataframe of de novos in the DDD dataset
    """
    
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

def get_enrichment_in_sample(n_trios, probands, de_novos, rates, threshold, genes):
    """ get the dominant DDG2P genes reaching genomewide significance
    
    Args:
        n_trios: number of trios to sample from among the probands
        probands: dataframe of proband IDs and sex info for DDD probands in trios
        de_novos: dataframe of de novos in the DDD dataset
        rates: dataframe of mutation rates per gene
        threshold: multiple testing corrected threshold for genomewide significance
        genes: vector of dominant DDG2P genes
    
    Returns:
        list of HGNC symbols for DDG2P genes reaching genomewide significance
    """
    
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

def run_iterations(probands, de_novos, rates, threshold, dominant,
        size, iterations, subset=None):
    """ detect genes reaching genomewide significance, under specific conditions
    
    Given a budget, relative cost of exome sequencing to genome sequencing, and
    relative sensitivity of genome sequencing to detect de novo mutations,
    determine the number of genes that could be detected by sequencing as many
    trios as the budget allows. This function separately checks many random
    samples of probands used for each test, so as to capture the variation in the
    population.
    
    Args:
        probands: dataframe of proband IDs and sex info for DDD probands in trios
        de_novos: dataframe of de novos in the DDD dataset
        rates: dataframe of mutation rates per gene
        threshold: multiple testing corrected threshold for genomewide significance
        dominant: vector of dominant DDG2P genes
        size: size of cohort to randomly sample
        iterations: number of iterations to run for the current conditions
        subset: list of person IDs to subset to, or None
    
    Returns:
        dataframe of test conditions, numbers of genes reaching genomewide
        significance from genome and exome sequencing
    """
    
    if subset is not None:
        probands = probands[probands['individual_id'].isin(subset)].copy()
        de_novos = de_novos[de_novos['person_id'].isin(subset)].copy()
    
    power = []
    
    for n in range(iterations):
        n = None
        genes = get_enrichment_in_sample(size, probands, de_novos, rates,
            threshold, dominant)
        if genes is not None:
            n = len(genes)
        
        power.append(n)
    
    return power

def plot_power(power, output_path):
    """ plot the results from simulation power of exome and genome sequencing
    
    Args:
        power: dataframe of power simulations for each condition, containing
            numbers of genes reaching genomewide significance from genome and
            exome sequencing
        output_path: path to save output plot to
    """
    
    power['count'] = power['count'].astype(int)
    
    fig = seaborn.factorplot(x="size", y="count", hue="subset", kind='box',
        data=power, size=8)
    
    fig.savefig(output_path, format="pdf", bbox_inches='tight', pad_inches=0,
        transparent=True)

def main():
    args = get_options()
    
    rates = get_ddd_rates(args.rates)
    known = open_known_genes(args.known)
    dominant = sorted(set(known["gene"][known["mode"].isin(["Monoallelic", "X-linked dominant"])]))
    seizure_subset = pandas.read_table(args.seizures)
    seizure_subset = set(seizure_subset['person_stable_id'])
    
    # set the multiple testing corrected threshold for genomewide significance
    alpha = 0.01
    num_genes = 18500
    num_tests = num_genes * 2
    threshold = alpha/num_tests
    
    probands = get_probands(args.families, args.trios)
    de_novos = get_de_novos(args.de_novos, args.validations)
    
    power = pandas.DataFrame(columns=['subset', 'size', 'count'])
    for size in [200, 400, 600, 800]:
        print(size)
        
        complete_power = run_iterations(probands, de_novos, rates, threshold,
            dominant, size, args.iterations)
        seizure_power = run_iterations(probands, de_novos, rates, threshold,
            dominant, size, args.iterations, seizure_subset)
        
        seizure_power = pandas.DataFrame({'subset': ['seizure'] * len(seizure_power),
            'size': [size] * len(seizure_power), 'count': seizure_power})
        complete_power = pandas.DataFrame({'subset': ['complete'] * len(complete_power),
            'size': [size] * len(complete_power), 'count': complete_power})
        
        power = power.append(seizure_power, ignore_index=True)
        power = power.append(complete_power, ignore_index=True)
    
    plot_power(power, args.output)

if __name__ == '__main__':
    main()
