"""
Copyright (c) 2015 Wellcome Trust Sanger Institute

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

import numpy

import matplotlib
matplotlib.use('Agg')
import seaborn

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes, \
    open_families
from ddd_4k.count_hpo import count_hpo_terms
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, PHENOTYPES, SANGER_IDS, \
    FAMILIES, DATATYPES
from ddd_4k.count_mutations_per_person import get_count_by_person
from ddd_4k.convert_durations import get_duration
from ddd_4k.scale_durations import autoscale_durations

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def plot_categorical(counts, pheno, value):
    """ plot probands per category by consequence by known gene status
    """
    
    merged = counts.merge(pheno[["person_stable_id", value]], on="person_stable_id")
    
    fig = seaborn.factorplot(x="consequence", hue=value, col="known", data=merged, size=6, kind="count")
    fig.set_ylabels("Frequency")
    fig.savefig("results/{}_by_consequence.pdf".format(value), format="pdf")
    matplotlib.pyplot.close()

def plot_quantitative(counts, pheno, value, y_label):
    """ plot quantitative metric by functional category by known gene status
    """
    
    merged = counts.merge(pheno[["person_stable_id", value]], on="person_stable_id")
    
    fig = seaborn.factorplot(x="known", y=value, hue="consequence", size=6, data=merged, kind="violin")
    fig.set_ylabels(y_label)
    fig.savefig("results/{}_by_consequence.pdf".format(value), format="pdf")
    matplotlib.pyplot.close()

def plot_hpo_by_consequence(counts, pheno):
    """ Plot number of HPO terms by functional category by known gene status
    """
    
    pheno["child_hpo_n"] = count_hpo_terms(pheno, "child")
    hpo_counts = counts.merge(pheno[["person_stable_id", "child_hpo_n"]], on="person_stable_id")
    
    fig = seaborn.factorplot(x="known", y="child_hpo_n", hue="consequence", size=6, data=hpo_counts, kind="violin")
    fig.set_ylabels("HPO terms per proband (n)")
    fig.savefig("results/hpo_by_consequence.pdf", format="pdf")
    matplotlib.pyplot.close()

def plot_achievement(counts, pheno, achievement):
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
    plot_quantitative(counts, pheno, achievement, "{} log10({}s)".format(achievement, unit))

def main():
    de_novos = open_de_novos(DENOVO_PATH)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    
    # For each proband, count the number of functional de novos (split by lof
    # and missense), in known developmental disorder genes (and other genes).
    counts = get_count_by_person(de_novos)
    
    # load the phenotype data
    pheno = open_phenotypes(PHENOTYPES, SANGER_IDS)
    
    # restrict the phenotype dataset to the probands for whom we could have de novo
    # candidates, that is the set of probands where all the members of their trio
    # have exome sequence data available.
    families = open_families(FAMILIES, DATATYPES)
    probands = families["individual_id"][(families["dng"] == 1)]
    pheno = pheno[pheno["person_stable_id"].isin(probands)]
    
    plot_hpo_by_consequence(counts, pheno)
    plot_categorical(counts, pheno, "gender")
    plot_categorical(counts, pheno, "scbu_nicu")
    plot_categorical(counts, pheno, "feeding_problems")
    plot_categorical(counts, pheno, "maternal_illness")
    
    # birthweight (or birthweight corrected for duration of gestation, or
    # birthweight_percentile (if that corrects for duration of gestation))
    plot_quantitative(counts, pheno, "decimal_age_at_assessment", "Age at assessment (years)")
    plot_quantitative(counts, pheno, "birthweight", "Birthweight (grams)")
    plot_quantitative(counts, pheno, "gestation", "Gestation duration (weeks)")
    plot_quantitative(counts, pheno, "height_percentile", "Height (percentile)")
    plot_quantitative(counts, pheno, "weight_percentile", "weight (percentile)")
    plot_quantitative(counts, pheno, "ofc_percentile", "OFC (percentile)")
    
    # plot distributions of time to achieve developmental milestones by
    # the functional categories
    plot_achievement(counts, pheno, "social_smile")
    plot_achievement(counts, pheno, "sat_independently")
    plot_achievement(counts, pheno, "walked_independently")
    plot_achievement(counts, pheno, "first_words")
    

if __name__ == '__main__':
    main()
