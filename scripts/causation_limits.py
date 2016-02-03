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

from __future__ import division

import numpy
import pandas

import matplotlib
matplotlib.use('Agg')
import seaborn

from ddd_4k.load_files import open_de_novos, open_known_genes
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, VALIDATIONS, CONSTRAINTS_URL
from ddd_4k.count_mutations_per_person import get_count_by_person

from mupit.mutation_rates import get_default_rates, get_expected_mutations
from mupit.count_de_novos import get_de_novo_counts
from mupit.gene_enrichment import gene_enrichment

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

HAPLOINSUFFICIENCY_URL = "http://files.figshare.com/410746/Dataset_S1.txt"

def open_haploinsufficiency(url):
    """ load the haploinsufficiency dataset
    
    Args:
        url: URL for haploinsufficiency dataset.
    
    Returns:
        pandas DataFrame of haploinsufficiency scores by gene symbol.
    """
    
    # get haploinsufficiency predictions genomewide. The data is from a report
    # in PLOS Genetics, formatted as a bed file:
    # PLOS Genet 6:e1001154 - doi:10.1371/journal.pgen.1001154
    # TODO: look for updated list of haploinsufficiency, with newer gene symbols
    # and get better threshold to classify haploinsufficiency.
    hi = pandas.read_table(url, skiprows=1, header=None)
    hi.columns = ["chrom", "start", "end", "name", "score", "strand",
        "thick_start", "thick_end", "color"]
    
    # reformat the gene and chromosome from the table
    hi["gene"] = [ x.split("|")[0] for x in hi["name"] ]
    hi["chrom"] = [ x.strip("chr") for x in hi["chrom"] ]
    
    return hi

def get_enrichment_ratios(constraints, haploinsufficiency):
    """ figure out the haploinsufficiency enrichment by constraint bin
    
    Args:
        constraints: pandas DataFrame of constraint scores by gene
        haploinsufficiency: pandas DataFrame of haploinsufficiency scores by gene
    
    Returns:
        pandas DataFrame of enrichment ratios for each pLI bin.
    """
    
    # classify each gene as belonging to one of 20 evenly spaced bins
    constraints["bin"], bins = pandas.cut(constraints["pLI"], bins=20, retbins=True)
    
    merged = constraints.merge(haploinsufficiency, on="gene")
    groups = merged.groupby("bin")
    
    # determine the enrichment ratio in each pLI bin
    threshold = 0.9
    values = [ sum(x["score"] > threshold)/len(x.index) for i, x in groups ]
    ratios = [ x/values[0] for x in values ]
    enrichment = pandas.DataFrame({"pLI": bins[:-1], "enrichment": ratios})
    
    return enrichment

def main():
    de_novos = open_de_novos(DENOVO_PATH, VALIDATIONS)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])
    
    counts = get_de_novo_counts(de_novos)
    
    # # For each proband, count the number of functional de novos (split by lof
    # # and missense), in known developmental disorder genes (and other genes).
    # counts = get_count_by_person(de_novos)
    
    male = 2407
    female = 1887
    rates = get_default_rates()
    expected = get_expected_mutations(rates, male, female)
    
    enriched = gene_enrichment(expected, counts)
    
    ratios = expected[["hgnc", "lof_indel", "lof_snv", "missense_indel", "missense_snv"]].copy()
    ratios["lof_ratio"] = numpy.nan
    ratios["mis_ratio"] = numpy.nan
    
    for (key, gene_exp) in ratios.iterrows():
        print(key)
        hgnc = gene_exp["hgnc"]
        # gene_obs = counts[counts["hgnc"] == hgnc]
        gene_obs = counts[counts["hgnc"] == hgnc]
        
        if len(gene_obs) == 0:
            gene_obs = pandas.DataFrame({"lof_indel": [0], "lof_snv": [0],
                "missense_indel": [0], "missense_snv": [0]})
        
        lof_exp = gene_obs[["lof_indel", "lof_snv"]]
        lof_obs = gene_exp[["lof_indel", "lof_snv"]]
        lof_ratio = float(lof_exp.sum(axis=1))/float(lof_obs.sum(axis=1))
        
        mis_exp = gene_obs[["missense_indel", "missense_snv"]]
        mis_obs = gene_exp[["missense_indel", "missense_snv"]]
        mis_ratio = float(mis_exp.sum(axis=1))/float(mis_obs.sum(axis=1))
        
        ratios.ix[key, "lof_ratio"] = lof_ratio
        ratios.ix[key, "mis_ratio"] = mis_ratio
    
    constraints = pandas.read_table(CONSTRAINTS_URL)
    recode = dict(zip(constraints["gene"], constraints["pLI"]))
    
    ratios["pLI"] = ratios["hgnc"].map(recode)
    
    quantiles = [ x/20 for x in range(21) ]
    ratios["bin"], bins = pandas.qcut(ratios["pLI"], q=quantiles, labels=quantiles[:-1], retbins=True)
    
    lof_enrich = pandas.pivot_table(ratios, values="lof_ratio", rows=["bin"], aggfunc=numpy.mean, fill_value=0)
    mis_enrich = pandas.pivot_table(ratios, values="mis_ratio", rows=["bin"], aggfunc=numpy.mean, fill_value=0)
    
    # haploinsufficiency = open_haploinsufficiency(HAPLOINSUFFICIENCY_URL)
    # enrichment = get_enrichment_ratios(constraints, haploinsufficiency)
    
    # # plot the enrichment ratio for each bin
    # fig = seaborn.lmplot(x="pLI", y="enrichment", data=enrichment, size=6, lowess=True)
    # fig.savefig("results/enrichment_ratio.pdf", type="pdf")

if __name__ == '__main__':
    main()
