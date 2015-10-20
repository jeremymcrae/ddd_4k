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

import os
import argparse

import matplotlib
matplotlib.use('Agg')
import seaborn
import pandas
from scipy.stats import fisher_exact, linregress
from matplotlib import pyplot
from numpy import median, log10, mean

from ddd_4k.constants import TRIOS, DIAGNOSED

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

AUTOZYGOSITY_DIR = "/nfs/users/nfs_j/jm33/apps/recessiveStats/data-raw/autozygosity"
CONSANG_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/kinship_and_pca_trios.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check whether"
        "probands with longer autozygous regions are less likely to have"
        "diagnoses")
    parser.add_argument("--autozygosity-dir", default=AUTOZYGOSITY_DIR, \
        help="Path to folder with autozygosity predictions for all probands.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Path to file of probands with diagnoses or likely diagnoses.")
    parser.add_argument("--diagnosed", default=DIAGNOSED, \
        help="Path to file of probands with diagnoses or likely diagnoses.")
    parser.add_argument("--consanguinous", default=CONSANG_PATH, \
        help="Path to all probands table that includes kinship statistics.")
    parser.add_argument("--output-groups", default="results/autozygosity_by_diagnosed.groups.pdf", \
        help="Path to plot regression graph to.")
    parser.add_argument("--output-regression", default="results/autozygosity_by_diagnosed.regression.pdf", \
        help="Path to plot regression graph to.")
    
    args = parser.parse_args()
    
    return args

def get_consanguinous(path):
    """ get a list of probands with consanguinous parents
    
    Args:
        path: path to table that includes kinship statistics.
    
    Returns:
        pandas Series of proband IDs for the proabnds from consanguinous unions.
    """
    
    consang = pandas.read_table(path, sep="\t")
    
    consang = consang[consang["king_kinship"] > 0]
    
    return consang["proband_stable_id"]

def get_autozygous_lengths(trios_path, folder):
    """ get the total length of the autozygous regions in each proband
    
    Args:
        trios_path: path to PED file for the cohort, detailing family
            relationships.
        folder: folder containing files per proband, each listing autozygous
            regions for the proband.
    
    Returns:
        pandas dataframe of person IDs and autozygosity lengths for all probands.
    """
    
    trios = pandas.read_table(trios_path, sep="\t")
    
    probands = trios["proband_stable_id"]
    
    sample_ids = []
    lengths = []
    for proband in probands:
        path = os.path.join(folder, proband)
        
        if os.path.exists(path):
            sample = pandas.read_table(path, sep="\t")
            sample["chrom"] = sample["chrom"].astype(str)
            
            # figure out the length of each autozygous region
            sample = sample[sample["chrom"] != "X"]
            sample["delta"] = sample["end_pos"] - sample["start_pos"]
            
            sample_ids.append(proband)
            lengths.append(sum(sample["delta"]))
    
    return pandas.DataFrame({"person_id": sample_ids, "length": lengths})

def define_quintiles(lengths, n_bins):
    """ find the break points for various quantiles of autozygosity length
    
    Args:
        lengths: DataFrame of autozygosity lengths and person IDs
        n_bins: number of bins to split the dataset into
    
    Returns:
        list of quantile values (including zero as the first value)
    """
    
    # get the non-zero autozygous lengths. We use non-zero, as the bulk of the
    # cohort has zero autozygous regions, which would distort a quantile-based
    # assessment if they were included.
    regions = lengths[lengths["length"] > 0]
    
    # figure out the quantiles for the number of bins, then get the length at
    # each quantile.
    quantiles = [ x/n_bins for x in range(n_bins + 1) ]
    quintiles = [ regions["length"].quantile(x) for x in quantiles ]
    
    return quintiles

def classify_by_quintile(lengths, quintiles):
    """ identify which quintile of autozygosity length each proband falls into
    
    Args:
        lengths: DataFrame of autozygosity lengths and person IDs
        quintiles: break points for the required quintiles
    
    Returns:
        dataframe, with an extra column added, to indicate which quintile each
        proband falls into.
    """
    
    # split the probands into a set who don't have any autozygous regions, vs
    # those who do. this is because the bulk of probands don't have any predicted
    # autozygosity, which would interfere with splitting the cohort into quintiles
    # based on their autozygosity.
    no_region = lengths[lengths["length"] == 0]
    has_region = lengths[lengths["length"] > 0]
    
    # define the quintile for the probands without autozygosity as one that will
    # sort ahead of the quintiles from the probands with autozygosity
    no_region["quintile"] = "no region"
    
    value = pandas.Series([0.0] * len(has_region), index=has_region.index)
    for pos, low in enumerate(quintiles[:-1]):
        high = quintiles[pos + 1]
        
        label = "{0:.1f}-{1:.1f}Mb".format(low/1e6, high/1e6)
        value[(has_region["length"] >= low) & (has_region["length"] <= high)] = label
    
    has_region["quintile"] = value
    
    lengths = no_region.append(has_region)
    
    return lengths

def autozygosity_vs_diagnosed(lengths, diagnosed_ids, plot_path):
    """ plot rates of autozygosity versus the likelihood of having a diagnosis
    
    Args:
        lengths: DataFrame of autozygosity lengths, person IDs, quintile
        diagnosed_ids: pandas Series of person IDs for the probands who have
            diagnoses.
        plot_path: path to plot graph to as pdf.
    """
    
    # figure out the correct order of the quintiles
    labels = lengths["quintile"].unique()
    values = [ float(x.split("-")[0]) if "-" in x else -1 for x in labels ]
    positions = [ values.index(x) for x in sorted(values) ]
    labels = [ labels[x] for x in positions ]
    
    # assess whether each proband has a diagnosis
    lengths["diagnosed"] = lengths["person_id"].isin(diagnosed_ids)
    lengths["diagnosed"] = lengths["diagnosed"].map({True: 1, False: 0})
    
    fig = seaborn.factorplot(x="quintile", y="diagnosed", \
        data=lengths, order=labels, join=False, ci=95, kind="point", size=6, \
        aspect=1.2, legend_out=False)
    fig.set_xticklabels(rotation=90, fontsize="x-large")
    fig.savefig(plot_path, format="pdf")

def plot_regression(lengths, diagnosed_ids, plot_path):
    """ plot regression of autozygosity length by diagnostic likelihood
    
    Args:
        lengths: DataFrame of autozygosity lengths, person IDs, quintile
        diagnosed_ids: pandas Series of person IDs for the probands who have
            diagnoses.
        plot_path: path to plot graph to as pdf.
    """
    
    # assess whether each proband has a diagnosis
    lengths["diagnosed"] = lengths["person_id"].isin(diagnosed_ids)
    lengths["diagnosed"] = lengths["diagnosed"].map({True: 1, False: 0})
    
    lengths = lengths[lengths["length"] > 0]
    new_lengths = pandas.DataFrame(columns=lengths.columns + ["median_length"])
    for (key, x) in lengths.groupby("quintile"):
        if key != "no region":
            median_length = log10(median(x["length"]))
            x["median_length"] = median_length
            new_lengths = new_lengths.append(x, ignore_index=True)
    
    # get the linear regression parameters
    groups = new_lengths.groupby("median_length")
    data = [ (key, sum(x["diagnosed"])/len(x)) for (key, x) in groups ]
    model = linregress(zip(*data))
    
    slope = model[0]
    r_squared = model[2]**2
    
    fig = seaborn.lmplot(x="median_length", y="diagnosed", \
        data=new_lengths, x_estimator=mean, size=6, aspect=1.2)
        
    text = "For every log10-unit increase in autozygous\nlength, the " \
        "diagnostic probability drops by {:.1f}%".format(abs(slope)*100)
    fig.fig.text(0.4, 0.9, text, fontsize="large")
    fig.fig.text(0.8, 0.85, "r^2={0:.3f}".format(r_squared), fontsize="large")
    fig.savefig(plot_path, format="pdf")

def main():
    args = get_options()
    lengths = get_autozygous_lengths(args.trios, args.autozygosity_dir)
    
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    diagnosed = diagnosed[diagnosed["inheritance"] == "de_novo"]
    diagnosed = diagnosed["person_id"].unique()
    
    # breakpoints = define_quintiles(lengths, 5)
    breakpoints = [0, 1e7, 1e8, 1e9]
    
    # classify each proband's autozygosity length into the various quintiles,
    # for all probands
    lengths = classify_by_quintile(lengths, breakpoints)
    
    autozygosity_vs_diagnosed(lengths, diagnosed, args.output_groups)
    plot_regression(lengths, diagnosed, args.output_regression)


if __name__ == '__main__':
    main()
