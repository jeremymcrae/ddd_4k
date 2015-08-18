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

from ddd_4k.load_de_novos import open_de_novos, open_known_genes, open_phenotypes

user_dir = os.path.expanduser("~")

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

DENOVO_PATH = "{}/apps/denovoFilter/de_novos.ddd_4k.ddd_only.txt".format(user_dir)
KNOWN_GENES = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/ddg2p_freeze_Jul15_corrected2_with_gencode_v19_coordinates_fixed.txt"
PHENOTYPES = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/phenotypes_and_patient_info.txt"
SANGER_IDS = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/person_sanger_decipher.txt"

de_novos = open_de_novos(DENOVO_PATH)
known = open_known_genes(KNOWN_GENES)
de_novos["known"] = de_novos["symbol"].isin(known["gencode_gene_name"])

pheno = open_phenotypes(PHENOTYPES, SANGER_IDS)



# count number of HPO terms per proband
