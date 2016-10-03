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

import pandas
import numpy
from scipy.stats import norm

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, gridspec
import seaborn

from ddd_4k.load_files import open_de_novos

seaborn.set_context('notebook', font_scale=2)
seaborn.set_style('white', {'ytick.major.size': 10, "xtick.major.size": 10})

de_novos_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt'
validation_path = '/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-11-24.txt'
seizure_probands_path = '/nfs/users/nfs_j/jm33/seizure_probands.txt'
trios_path = '/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt'

def load_seizure_genes():
    GENE_PANEL = "/nfs/users/nfs_j/jm33/apps/seizure/data-raw/epilepsy_gene_tests.txt"
    EPILEPSY_GENE_TESTS = "/nfs/users/nfs_j/jm33/apps/seizure/data-raw/GOSH Epilepsy Panel DDG2P coded.tab"
    OHTAHARA_GENES = "/nfs/users/nfs_j/jm33/apps/seizure/data-raw/known_ohtahara_genes.txt"
    
    panel = pandas.read_table(GENE_PANEL)
    epilepsy = pandas.read_table(EPILEPSY_GENE_TESTS, header=None)
    ohtahara = pandas.read_table(OHTAHARA_GENES)
    
    ohtahara = list(ohtahara['gene_name'])
    epilepsy = list(epilepsy[0][~epilepsy[0].isnull()])
    panel = list(panel['hgnc_symbol'])
    
    return set(panel + epilepsy + ohtahara)

def main():
    de_novos = open_de_novos(de_novos_path, validation_path)
    probands = pandas.read_table(seizure_probands_path)
    trios = pandas.read_table(trios_path)
    
    genes = load_seizure_genes()
    de_novos['seizure_gene'] = de_novos['hgnc'].isin(genes)
    seizure_de_novos = de_novos[de_novos['seizure_gene']]
    
    trios['has_seizures'] = trios['proband_stable_id'].isin(probands['person_stable_id'])
    trios['has_seizure_de_novo'] = trios['proband_stable_id'].isin(seizure_de_novos['person_stable_id'])
    
    tab = pandas.crosstab(trios['has_seizure_de_novo'], trios['has_seizures'])
    print(tab)
    
    fig = seaborn.factorplot(x='has_seizures', y='has_seizure_de_novo',
        data=trios, join=False, ci=95, kind='point', size=6,
        aspect=1.2, legend_out=False, capsize=0.2)
    fig.set_xlabels('Has seizures')
    fig.set_ylabels('Fraction with a DNM in a\nseizure-associated gene')
    ax = fig.ax
    
    ax.set_ylim(0, max(ax.get_ylim()))
    
    fig.savefig('seizure_probands_with_seizure_de_novos.pdf', format='pdf', bbox_inches='tight', pad_inches=0, transparent=True)

if __name__ == '__main__':
    main()
