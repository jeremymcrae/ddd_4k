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

import os
import argparse

import pandas
import numpy

from denovonear.gene_plot.diagram_plotter import DiagramPlotter
from denovonear.gene_plot.interpro_requestor import InterproRequest
from denovonear.transcript import Transcript
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_gene import load_gene

from ddd_4k.constants import THRESHOLD, DENOVO_PATH, VALIDATIONS, DIAGNOSED
from ddd_4k.load_files import open_de_novos

ASSOCIATIONS = "/lustre/scratch113/projects/ddd/users/jm33/results/" \
    "de_novos.ddd_4k.without_diagnosed.all.2015-10-12.txt"
EXTERNAL_DE_NOVOS = "/nfs/users/nfs_j/jm33/apps/publishedDeNovos/data-raw/" \
    "variants.txt.gz"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--results", default=ASSOCIATIONS, \
        help="Path to folder with association results.")
    parser.add_argument("--de-novos", default=DENOVO_PATH, \
        help="Path to file of candidate de novos.")
    parser.add_argument("--validations", default=VALIDATIONS, \
        help="Path to table of de novo validation results.")
    parser.add_argument("--diagnosed", default=DIAGNOSED, \
        help="Path to table of probands with other diagnoses.")
    parser.add_argument("--external-sites", default=EXTERNAL_DE_NOVOS, \
        help="Path to table of de novos from reported exome or" \
            "genome-sequencing studies.")
    parser.add_argument("--output-dir", default="gene_plots", \
        help="path to folder to place gene plots in.")
    
    args = parser.parse_args()
    
    return args

def get_domains(protein_sequence):
    """ get protein domains from the protein sequence
    """
    
    requestor = InterproRequest("jeremy.mcrae@sanger.ac.uk")
    results = requestor.submit_interpro_run(protein_sequence)
    
    return results

def standardise_ddd_de_novos(de_novos):
    """standardise the columns in the DDD de novos, so that we can easily
    integrate this with the externally reported de novos
    """
    
    de_novos["person_id"] = de_novos["person_stable_id"]
    de_novos["hgnc"] = de_novos["symbol"]
    de_novos["ref_allele"] = de_novos["ref"]
    de_novos["alt_allele"] = de_novos["alt"]
    de_novos["start_pos"] = de_novos["pos"]
    de_novos["end_pos"] = de_novos["start_pos"] + de_novos["ref_allele"].str.len() - 1
    de_novos["sex"] = de_novos["sex"].map({"F": "female", "M": "male"})
    de_novos["type"] = de_novos["var_type"].map({"DENOVO-SNP": "snv", \
        "DENOVO-INDEL": "indel"})
    de_novos["study_code"] = "ddd_unpublished"
    de_novos["publication_doi"] = numpy.nan
    de_novos["study_phenotype"] = "developmental_disorders"
    
    # only select the columns to match the extenal de novos
    de_novos = de_novos[["person_id", "sex", "chrom", "start_pos", "end_pos",
        "ref_allele", "alt_allele", "hgnc", "consequence", "study_code",
        "publication_doi", "study_phenotype", "type"]]
    
    return de_novos

def load_de_novos(de_novo_path, validation_path, external_path, diagnosed_path):
    """ loads the
    """
    
    de_novos = open_de_novos(de_novo_path, validation_path)
    de_novos = standardise_ddd_de_novos(de_novos)
    
    external = pandas.read_table(external_path, compression="gzip")
    de_novos = de_novos.append(external, ignore_index=True)
    
    diagnosed = pandas.read_table(diagnosed_path)
    de_novos = de_novos[~de_novos["person_id"].isin(diagnosed["person_id"])]
    
    return de_novos

def variants_table_to_dictionary(variants):
    """ convert the variants to a dictionary.
    
    This is so we can easily iterate through the variants, and use the
    attributes as necessary. We need to use a tuple of position and person_id
    as the key, since the position alone might not be unique.
    """
    
    var_dict = {}
    for idx, var in variants.iterrows():
        key = (var["start_pos"], var["person_id"])
        var_dict[key] = {}
        var_dict[key]["start_pos"] = var["start_pos"]
        var_dict[key]["ref_allele"] = var["ref_allele"]
        var_dict[key]["alt_allele"] = var["alt_allele"]
        var_dict[key]["consequence"] = var["consequence"]
        
        var_dict[key]["source"] = "external"
        if var["study_code"] == "ddd_unpublished":
            var_dict[key]["source"] = "internal"
    
    return var_dict

def main():
    args = get_options()
    
    de_novos = load_de_novos(args.de_novos, args.validations,
        args.external_sites, args.diagnosed)
    
    genes = pandas.read_table(args.results)
    genes = genes[genes["p_min"] < THRESHOLD]
    
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    
    ensembl = EnsemblRequest(cache_folder="cache", genome_build="grch37")
    for hgnc in genes["hgnc"]:
        
        variants = de_novos[de_novos["hgnc"] == hgnc]
        
        variants = variants_table_to_dictionary(variants)
        start_positions = [ v["start_pos"] for k, v in variants.items() ]
        
        transcripts = load_gene(ensembl, hgnc, start_positions)
        transcript = transcripts[0]
        protein_sequence = ensembl.get_protein_seq_for_transcript(transcript.get_name())
        
        filename = "{0}_gene_plot.pdf".format(hgnc)
        if args.output_dir is not None:
            filename = os.path.join(args.output_dir, filename)
        
        plotter = DiagramPlotter(transcript, hgnc, variants, filename)
        plotter.plot_gene()
        plotter.plot_transcript()
        plotter.plot_domains(get_domains(protein_sequence), protein_sequence)
        plotter.export_figure()

if __name__ == '__main__':
    main()
