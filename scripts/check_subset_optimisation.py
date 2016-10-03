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

from __future__ import division, print_function

import os
import argparse
import subprocess
import tempfile
import time
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from scipy.stats import wilcoxon
from numpy import log10, median
import pandas

from ddd_4k.combine_enrichment_and_clustering import get_gene_results
from mupit.open_ddd_data import open_known_genes
from ddd_4k.convert_doi import open_url
from ddd_4k.constants import CONSTRAINTS_URL

RATES_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.mutation_rates.2015-11-24.txt"
FILTERED_DE_NOVOS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-11-24.txt"
VALIDATIONS_PATH = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-11-24.txt"
FAMILIES_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt"
TRIOS_PATH = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt"
DDG2P_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
TEMP_DIR = "/lustre/scratch113/projects/ddd/users/jm33/temp"
DENOVONEAR_RATES = "denovonear/data/forSanger_1KG_mutation_rate_table.txt"
INITIAL_ENRICH = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.ddd_only.enrichment_results.2015-11-24.txt"
INITIAL_CLUSTER = "/lustre/scratch113/projects/ddd/users/jm33/results/de_novos.ddd_4k.with_diagnosed.ddd_only.clustering_results.txt"
EXTERNAL_COHORTS = "publishedDeNovos/data-raw/cohorts.tsv"

r_binary = "/software/R-3.2.2/bin/Rscript"
mupit = "mupit/scripts/ddd_analysis.R"
denovonear_batch = "denovonear/scripts/run_batch.py"
denovonear = "denovonear/scripts/clustering.py"

def get_pli_scores():
    headers = {}
    url = "{}/export?format=csv".format(pLI_URL)
    response, status_code, headers = open_url(url, headers)
    pLI = pandas.read_csv(StringIO(response))
    
    return pLI

def rate_limit(wait_time=60):
    """ temporarily sleep if we have tried to check the jobs too recently
    
    Args:
        wait_time: time to wait between checks
    """
    
    if "prev_time" not in globals():
        global prev_time
        prev_time = time.time() - wait_time
    
    current_time = time.time()
    delta = current_time - prev_time
    
    if delta < wait_time:
        time.sleep(wait_time - delta)
    
    prev_time = current_time

def get_jobs():
    """ get a set of job IDs which are currently on the cluster
    
    Returns:
        set of job names for currently running jobs
    """
    
    rate_limit()
    command = ["bjobs", "-o", "\"JOBID", "USER", "STAT", "QUEUE", "JOB_NAME", \
        "delimiter=';'\""]
    command = " ".join(command)
    jobs = subprocess.check_output(command, shell=True, stderr=open(os.devnull))
    
    # if there aren't any currently running or pending jobs, then the output
    if jobs == "":
        return set([])
    
    jobs = jobs.decode().strip().split("\n")
    
    current_jobs = set([])
    for line in jobs:
        if line.startswith("JOBID"): # ignore the header line
            continue
        
        line = line.split(";")
        job_name = line[4]
        current_jobs.add(job_name)
    
    return current_jobs

def run_mupit(subsets, enrich_path, clust_path):
    if type(subsets) == list:
        subsets = ",".join(subsets)
    
    command = [r_binary, mupit,
        "--rates", RATES_PATH,
        "--de-novos", FILTERED_DE_NOVOS_PATH,
        "--validations", VALIDATIONS_PATH,
        "--families", FAMILIES_PATH,
        "--trios", TRIOS_PATH,
        "--ddg2p", DDG2P_PATH,
        "--meta-analysis",
        "--meta-subset", subsets,
        "--out-enrichment", enrich_path,
        "--out-clustering", clust_path]
    
    call = subprocess.check_call(command, stdout=open(os.devnull, "w"),
        stderr=open(os.devnull, "w"))

def run_denovonear(infile, outfile):
    command = ["python", denovonear_batch,
        "--script", denovonear,
        "--temp-dir", TEMP_DIR,
        "--in", infile,
        "--rates", DENOVONEAR_RATES,
        "--out", outfile]
    
    call = subprocess.check_call(command)
    
    time.sleep(2)
    
    # wait until all the denovonear jobs have completed
    while len([ x for x in get_jobs() if "denovonear" in x ]) > 1:
        time.sleep(1)

def get_gene_sets(table, dominant):
    """ get the HGNC symbols for the genes which are genomewide or suggestive
    
    Args:
        table: pandas dataframe of test results, including columns for WHETHER
            the genes are genomewide significant, or suggestive significant
            (FDR < 1%)
    
    Returns:
        dictionary of sets of genes with genomewide or suggestive significance.
    """
    
    known = table[table["hgnc"].isin(dominant)]
    gwide = set(known["hgnc"][known["genomewide"]])
    sugg = set(known["hgnc"][known["suggestive"]])
    
    gene_sets = {"genomewide": gwide, "suggestive": sugg}
    
    return gene_sets

def check_subset(name, previous, cohorts, baseline, dominant):
    """
    """
    
    subset = [name] + previous
    
    # define the temporary files
    enrich = tempfile.NamedTemporaryFile(dir=TEMP_DIR)
    clust_in = tempfile.NamedTemporaryFile(dir=TEMP_DIR)
    cluster = tempfile.NamedTemporaryFile(dir=TEMP_DIR)
    
    # run the enrichment and clustering tests
    run_mupit(subset, enrich.name, clust_in.name)
    run_denovonear(clust_in.name, cluster.name)
    
    cohort_n = cohorts[["unique_male", "unique_female"]][cohorts["study_phenotype"] == name].sum().sum()

    current = get_gene_results(enrich.name, cluster.name)
    baseline_genes = get_gene_sets(baseline, dominant)
    current_genes = get_gene_sets(current, dominant)
    
    extra = current_genes["genomewide"] - baseline_genes["genomewide"]
    removed = baseline_genes["genomewide"] - current_genes["genomewide"]
    delta = len(extra)/cohort_n - len(removed)/cohort_n
    
    values = {"extra": extra, "removed": removed, "results": current, \
        "n": cohort_n, "delta": delta}
    
    return values

def main():
    # find the most constrained genes in the genome by pLI score, ~1800 genes
    pLI = pandas.read_table(CONSTRAINTS_URL)
    constrained = pLI["gene"][pLI["pLI"] >= 0.99]
    
    ddg2p = open_known_genes(DDG2P_PATH)
    dominant = ddg2p["gene"][ddg2p["mode"].isin(["Monoallelic", "X-linked dominant"])]
    
    cohorts = pandas.read_table(EXTERNAL_COHORTS)
    baseline = get_gene_results(INITIAL_ENRICH, INITIAL_CLUSTER)
    
    extra = []
    while len(extra) < len(cohorts["study_phenotype"].unique()):
        previous = [ x[0] for x in extra ]
        remaining = set(cohorts["study_phenotype"].unique()) - set(previous)
        
        # run through all the remaining cohorts, test the difference they each make
        temp = {}
        for name in remaining:
            temp[name] = check_subset(name, previous, cohorts, baseline, dominant)
            
            line = "{}\t{}\t{}\t{}\t{}".format(":".join(previous + [name]), \
                temp[name]["n"], ",".join(temp[name]["extra"]), ",".join(temp[name]["removed"]), \
                temp[name]["delta"])
            print(line)
        
        # find the subset that gave the greatest improvement
        max_delta = max([ temp[x]["delta"] for x in temp ])
        for (key, value) in temp.items():
            if value["delta"] == max_delta:
                break
        
        baseline = value["results"]
        extra.append((key, value))
    
    print("\nfinal_list\n")
    for key, value in extra:
        line = "{}\t{}\t{}\t{}\t{}".format(key, \
            value["n"], ",".join(value["extra"]), ",".join(value["removed"]), value["delta"])
        print(line)

if __name__ == '__main__':
    main()
