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
import sys
import tempfile

import pandas
import scipy
from numpy import log10
import matplotlib
matplotlib.use("Agg")

import seaborn

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import CalculateSimilarity

from ddd_4k.constants import VALIDATIONS, DENOVO_PATH
from ddd_4k.load_files import open_de_novos
from ddd_4k.convert_doi import open_url

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def load_neurodevelopmental(path):
    """ load the dominant LoF neurodevelopmental genes, with recognizability
    """
    
    neurodev = pandas.read_excel(path, sheetname="neurodevelopmental.dominant_lof")
    
    return neurodev

def load_rates(url):
    """ load the mutation rates for these genes
    """
    
    # the version of pandas I have cannot natively load excel files from urls.
    # So we have to downalod the file, save it, then read from the file on disk.
    response, status, header = open_url(url, headers={})
    
    if status != 200:
        sys.exit("Failed to load Samocha supplementary rates from url")
    
    with tempfile.NamedTemporaryFile(mode="wb") as handle:
        handle.write(response)
        handle.flush()
        rates = pandas.read_excel(handle.name, sheetname="mutation_probabilities")
    
    return rates

def count_de_novos(de_novo_path, validations_path):
    """ count the loss-of-function mutations for each gene.
    """
    
    de_novos = open_de_novos(de_novo_path, validations_path)
    # de_novos = de_novos[de_novos["category"] == "loss-of-function"]
    
    counts = {}
    for (gene, rows) in de_novos.groupby("hgnc"):
        counts[gene] = len(rows)
    
    return counts

def main():
    rates_url = "http://www.nature.com/ng/journal/v46/n9/extref/ng.3050-S2.xls"
    neurodev_path = "/nfs/users/nfs_j/jm33/neurodevelopmental.dominant_lof_DRF.xlsx"
    
    rates = load_rates(rates_url)
    rates["splice_site"][rates["splice_site"].isnull()] = -300
    rates["rate"] = log10(10**rates["mis"] + 10**rates["non"] +
        10**rates["splice_site"] + 10**rates["frameshift"])
    rates = rates[["gene", "rate"]]
    
    counts = count_de_novos(DENOVO_PATH, VALIDATIONS)
    neurodev = load_neurodevelopmental(neurodev_path)
    
    neurodev["count"] = 0
    for hgnc in counts:
        if hgnc in neurodev["hgnc"].values:
            neurodev["count"][neurodev["hgnc"] == hgnc] = counts[hgnc]
    
    # plot LoF counts by mutation rates, shade by recognizability
    neurodev = neurodev.merge(rates, how="left", left_on=["hgnc"], right_on=["gene"])
    
    fig = seaborn.lmplot(x="rate", y="count", hue="Recognisable", data=neurodev,
        size=6, hue_order=[1, 2, 3, 4, 5], y_jitter=0.1)
    fig.savefig("temp.pdf", format="pdf")

if __name__ == '__main__':
    main()
