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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import seaborn
import pandas

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes
from ddd_4k.count_hpo import count_hpo_terms
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, PHENOTYPES, SANGER_IDS

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

de_novos = open_de_novos(DENOVO_PATH)
known = open_known_genes(KNOWN_GENES)
de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])

pheno = open_phenotypes(PHENOTYPES, SANGER_IDS)
pheno["child_hpo_n"] = count_hpo_terms(pheno, "child")

# figure out if each proband has a loss-of-function de novo in a known
# developmental disorder gene

# plot age by functional category by known gene status

# count number of HPO terms per proband. Plot number of HPO terms by functional
# category by known gene status
