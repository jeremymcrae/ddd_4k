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

from __future__ import print_function, division

import os
import argparse

import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot
import pandas
import seaborn

from ddd_4k.constants import TRIOS, DIAGNOSED
from ddd_4k.rank_hpo import rank_terms

from hpo_similarity.ontology import Ontology
from hpo_similarity.similarity import ICSimilarity
from hpo_similarity.load_files import load_participants_hpo_terms

# define the plot style
seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--de-novos", \
        help="Path to table of variants in novel genes.")
    parser.add_argument("--phenotypes", \
        default="/lustre/scratch113/projects/ddd/users/jm33/" \
            "de_novos.ddd_4k.phenotypes_by_proband.json", \
        help="JSON file of HPO terms per proband.")
    parser.add_argument("--trios", default=TRIOS, \
        help="Path to table of alternate IDs for participants.")
    parser.add_argument("--diagnosed", default=DIAGNOSED, \
        help="Path to table of probands with diagnoses.")
    parser.add_argument("--output-dir", \
        help="Path to folder to save plots to.")
    
    args = parser.parse_args()
    
    return args

def plot_gene(graph, probands, proband_hpo, hgnc, trios, table, output_dir, \
    max_count=10, include_key=False):
    """ plot the terms relevant for a gene as a heatmap
    
    Args:
        graph: ICSimilarity object, used to look up and down the HPO graph.
        probands: list of proband IDs for the gene.
        proband_hpo: dictionary of HPO terms by proband for the full ddd 4k
            cohort, indexed by proband IDs.
        hgnc: HGNC symbol for the gene, used to determine a pdf filename.
        trios: table of probands in exome-sequenced, including decipher and DDD IDs.
        table: pandas DataFrame of relevant terms for the gene, one row per term.
        output_dir: folder to save the plot to.
        max_count: the highest number of probands in a gene, for the highest
            value that a heatmap can have. We pass this variable in to
            standardise the heatmap colors between genes, by using the same
            value across genes.
        include_key: whether to include a key for the heatmap.
    """
    
    decipher = trios[trios["proband_stable_id"].isin(probands)]
    decipher_ids = dict(zip(decipher["proband_stable_id"], decipher["decipher_id"]))
    
    columns = sorted(decipher_ids.values())
    data = pandas.DataFrame(columns=columns)
    
    for num, row in table.iterrows():
        term = row["term"]
        name = row["name"]
        count = row["count"]
        
        desc = set([term]) | graph.get_descendants(term)
        has_term = [ len(set(proband_hpo[x]) & desc) > 0 for x in probands ]
        has_term = [ [count] if x else [0] for x in has_term ]
        
        decipher = [ decipher_ids[x] for x in probands ]
        
        temp = pandas.DataFrame(dict(zip(decipher, has_term)), index=[name])
        data = data.append(temp)
    
    data = data[columns]
    data = data[data.columns].astype(float)
    
    # plot the heatmap
    ax = seaborn.heatmap(data, cmap="Blues", cbar=include_key, square=True, \
        vmin=0, vmax=max_count)
    fig = ax.get_figure()
    fig.savefig(os.path.join(output_dir, "{}_terms.pdf").format(hgnc), format="pdf")
    pyplot.close()

def main():
    args = get_options()
    
    # args.de_novos = "/lustre/scratch113/projects/ddd/users/jm33/results/novel_gene_variants.ddd_4k.2015-11-24.txt"
    variants = pandas.read_table(args.de_novos, sep="\t")
    diagnosed = pandas.read_table(args.diagnosed, sep="\t")
    
    variants = variants[~variants["person_stable_id"].isin(diagnosed["person_id"])]
    
    # open the phenotype data, and restrict it to the probands with complete trios
    trios = pandas.read_table(args.trios, sep="\t")
    proband_ids = set(trios["proband_stable_id"])
    
    # open the HPO ontology, so we can get the set of terms which are relevant
    # to each disorder
    hpo_ontology = Ontology(None)
    graph = hpo_ontology.get_graph()
    
    alt_node_ids = hpo_ontology.get_alt_ids()
    obsolete_ids = hpo_ontology.get_obsolete_ids()
    
    # load the HPO terms for the probands, and subset the probands to the ones
    # in the DDD 4K trios
    proband_hpo = load_participants_hpo_terms(args.phenotypes, \
        alt_node_ids, obsolete_ids)
    proband_hpo = { x: proband_hpo[x] for x in proband_hpo if x in proband_ids }
    
    graph.tally_hpo_terms(proband_hpo)
    
    max_count = max(variants["symbol"].value_counts())
    
    for hgnc in sorted(variants["symbol"].unique()):
        gene = variants[variants["symbol"] == hgnc]
        probands = list(gene["person_stable_id"])
        table = rank_terms(graph, probands, proband_hpo)
        
        if len(table) > 0:
            plot_gene(graph, probands, proband_hpo, hgnc, trios, table, \
                args.output_dir, max_count)
    
    # replot the final gene, but this time include the key.
    plot_gene(graph, probands, proband_hpo, hgnc, trios, table, \
        args.output_dir, max_count, include_key=True)

if __name__ == '__main__':
    main()
