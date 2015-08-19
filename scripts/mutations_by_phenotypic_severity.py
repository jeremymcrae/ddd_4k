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
import seaborn
import pandas

from ddd_4k.load_files import open_de_novos, open_known_genes, open_phenotypes, \
    open_families
from ddd_4k.count_hpo import count_hpo_terms
from ddd_4k.constants import DENOVO_PATH, KNOWN_GENES, PHENOTYPES, SANGER_IDS, \
    FAMILIES, DATATYPES
from ddd_4k.count_mutations_per_person import get_count_by_person

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

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
    
    # plot age by functional category by known gene status
    age_counts = counts.merge(pheno[["person_stable_id", "decimal_age_at_assessment"]], on="person_stable_id")
    fig = seaborn.factorplot(x="known", y="decimal_age_at_assessment", hue="variable", size=6, data=age_counts, kind="box")
    fig.savefig("results/age_by_consequence.pdf", format="pdf")
    
    # Plot number of HPO terms by functional category by known gene status
    pheno["child_hpo_n"] = count_hpo_terms(pheno, "child")
    hpo_counts = counts.merge(pheno[["person_stable_id", "child_hpo_n"]], on="person_stable_id")
    fig = seaborn.factorplot(x="known", y="child_hpo_n", hue="variable", size=6, data=hpo_counts, kind="box")
    fig.savefig("results/hpo_by_consequence.pdf", format="pdf")

if __name__ == '__main__':
    main()
