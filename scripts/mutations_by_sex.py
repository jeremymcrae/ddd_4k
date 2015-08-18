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

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import seaborn
import pandas

from ddd_4k.load_de_novos import open_known_genes

user_dir = os.path.expanduser("~")

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

DENOVO_PATH = "{}/apps/denovoFilter/de_novos.ddd_4k.ddd_only.txt".format(user_dir)
KNOWN_GENES = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/ddg2p_freeze_Jul15_corrected2_with_gencode_v19_coordinates_fixed.txt"

def open_known_genes(path):
    """ open the dataset of known developmental disorder genes
    
    Args:
        path: path to known developmental disorder genes data file.
    
    Returns:
        DataFrame for the known genes.
    """
    
    genes = pandas.read_table(path, sep="|", na_filter=False)
    genes = genes[genes["ddg2p_status"] != "Possible DD Gene"]
    
    return genes

def get_count_by_person(de_novos):
    """ count the number of functional and loss-of-function de novos per person
    
    Args:
        de_novos: DataFrame of de novo mutations
    
    Returns:
        DataFrame of counts per person, split by consequence type and whether
        the mutation is in a known developmental disorder gene.
    """
    
    by_category = de_novos.pivot_table(values="chrom", rows=["person_stable_id", "sex", "known"], cols=["category"], aggfunc=len)
    
    index = [ x[0] for x in zip(by_category.index) ]
    sample_ids = [ x[0] for x in index ]
    sex = [ x[1] for x in index ]
    known = [ x[2] for x in index ]
    
    by_category = pandas.DataFrame({"person_stable_id": sample_ids, "sex": sex, "known": known, "functional": by_category["functional"].values, "loss-of-function": by_category["loss-of-function"].values})
    
    by_category = pandas.melt(by_category, id_vars=["person_stable_id", "sex", "known"], value_vars=["functional", "loss-of-function"])
    by_category = by_category[by_category["value"].notnull()]
    
    return by_category

def main():
    de_novos = open_de_novos(DENOVO_PATH)
    known = open_known_genes(KNOWN_GENES)
    de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])

    person_counts = get_count_by_person(de_novos)

    # seaborn.factorplot(x="variable", y="value", hue="sex", col="known", data=person_counts, size=6, kind="bar")
    seaborn.factorplot("variable", hue="sex", col="known", data=person_counts, size=6, kind="count")
    pyplot.savefig("results/test.pdf", format="pdf")

if __name__ == '__main__':
    main()
