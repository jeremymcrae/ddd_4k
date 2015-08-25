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

from __future__ import division

import numpy
import pandas

import matplotlib
matplotlib.use('Agg')
import seaborn

from ddd_4k.load_files import open_de_novos, open_known_genes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES
from ddd_4k.count_mutations_per_person import get_count_by_person

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def open_constraints():
    """ open the constraints dataset
    """
    
    # We have received a constraints dataset from Kaitlin Samocha and Mark Daly.
    # This data is largely identical to a file available from ExAC
    # ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/forweb_cleaned_exac_r03_2015_03_16_z_data.txt
    # but with one additional column, pLI, which is their metric to predict
    # loss-of-function intolerance.
    
    constraints_path = "data/cleaned_exac_with_pLI_march16.txt"
    constraints = pandas.read_table(constraints_path)
    
    return constraints

def open_haploinsufficiency():
    """ load the haploinsufficiency dataset
    """
    
    # check haploinsufficiency enrichment by constraint, to determine constraint
    # priors
    # get haploinsufficiency predictions genomewide
    # PLOS Genet 6:e1001154 - doi:10.1371/journal.pgen.1001154
    # TODO: look for updated list of haploinsufficiency, with newer gene symbols
    # and get better threshold to classify haploinsufficiency.
    url = "http://files.figshare.com/410746/Dataset_S1.txt"
    hi = pandas.read_table(url, skiprows=1, header=None)
    hi.columns = ["chrom", "start", "end", "name", "score", "strand",
        "thick_start", "thick_end", "color"]
    
    # reformat the gene and chromosome from the table
    hi["gene"] = [ x.split("|")[0] for x in hi["name"] ]
    hi["chrom"] = [ x.strip("chr") for x in hi["chrom"] ]
    
    return hi

def get_enrichment_ratios(constraints, haploinsufficiency):
    """ figure out the fraction of genes which are haploinsufficient, relative to the lowest bin.
    """
    
    # classify each gene as belonging to one of 20 evenly spaced bins
    constraints["bin"], bins = pandas.cut(constraints["pLI"], bins=20, retbins=True)
    
    merged = constraints.merge(haploinsufficiency, on="gene")
    groups = merged.groupby("bin")
    
    # determine the enrichment ratio in each pLI bin
    threshold = 0.9
    values = [ sum(x["score"] > threshold)/len(x.index) for i, x in groups ]
    enrichment = [ x/values[0] for x in values ]
    enrichment = pandas.DataFrame({"pLI": bins[:-1], "enrichment": enrichment})
    
    return enrichment

def main():
    de_novos = open_de_novos(DENOVO_PATH)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    
    # For each proband, count the number of functional de novos (split by lof
    # and missense), in known developmental disorder genes (and other genes).
    counts = get_count_by_person(de_novos)
    
    constraints = open_constraints()
    haploinsufficiency = open_haploinsufficiency()
    enrichment = get_enrichment_ratios(constraints, haploinsufficiency)
    
    # plot the enrichment ratio for each bin
    fig = seaborn.lmplot(x="pLI", y="enrichment", data=enrichment, size=6, lowess=True)
    fig.savefig("results/enrichment_ratio.pdf", type="pdf")

if __name__ == '__main__':
    main()
