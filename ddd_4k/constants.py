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

import os

user_dir = "/lustre/scratch113/projects/ddd/users/jm33/"
datafreeze = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/"

DENOVO_PATH = "{}/de_novos.ddd_4k.ddd_only.2015-11-24.txt".format(user_dir)
VALIDATIONS = "{}/de_novos.validation_results.2015-11-24.txt".format(user_dir)
DIAGNOSED = "{}/ddd_4k.diagnosed.2015-11-24.txt".format(user_dir)
KNOWN_GENES = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/ddg2p_freeze_Jul15_corrected2_with_gencode_v19_coordinates_fixed.txt"
PHENOTYPES = "{}/phenotypes_and_patient_info.txt".format(datafreeze)
SANGER_IDS = "{}/person_sanger_decipher.txt".format(datafreeze)
FAMILIES = "{}/family_relationships.txt".format(datafreeze)
TRIOS = "{}/trios.txt".format(datafreeze)
DATATYPES = "{}/person_datatypes.txt".format(datafreeze)
PREVIOUS_VALIDATIONS = '/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv'

# define parameters for the multiple-testing corrected p-value threshold.
ALPHA = 0.01
NUM_GENES = 18500
META_ANALYSED_GENES = 6000
NUM_TESTS = NUM_GENES * 2 + META_ANALYSED_GENES * 2
THRESHOLD = ALPHA/NUM_TESTS

CONSTRAINTS_URL = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"

SEIZURE_GENES = ["ADSL", "ALDH7A1", "ALG13", "AMT", "ARFGEF2", "ARHGEF9",
    "ARX", "ATP1A3", "ATRX", "CACNA1A", "CACNB4", "CBL", "CDKL5", "CHD2",
    "CHRNA2", "CHRNA4", "CHRNB2", "CNTNAP2", "CREBBP", "CSNK1G1", "CSTB",
    "DCX", "DNM1", "DOCK7", "EHMT1", "EP300", "EPM2A", "FASN", "FLNA",
    "FOXG1", "GABBR2", "GABRA1", "GABRB3", "GABRD", "GABRG2", "GATAD2B",
    "GCSH", "GLRA1", "GLRB", "GRIN2A", "GRIN2B", "HCN1", "KCNA1", "KCNB1",
    "KCNJ10", "KCNMA1", "KCNQ2", "KCNQ3", "KCNT1", "KCTD7", "KIAA1279",
    "LGI1", "MAGI2", "MBD5", "MECP2", "MEF2C", "MOCS1", "MOCS2", "NHLRC1",
    "NRXN1", "PCDH19", "PCDH19", "PIGA", "PIGQ", "PLCB1", "PLCB1", "PNKP",
    "POLG", "PRICKLE1", "PRPH", "PRRT2", "QARS", "RYR3", "SCN1A", "SCN1B",
    "SCN2A", "SCN3A", "SCN4A", "SCN8A", "SCN9A", "SLC12A5", "SLC13A5",
    "SLC16A2", "SLC25A22", "SLC2A1", "SLC35A2", "SLC6A5", "SLC9A6",
    "SMARCA2", "SPTAN1", "ST3GAL5", "STXBP1", "SUOX", "SYNGAP1", "TBC1D24",
    "TCF4", "TSC1", "TSC2", "UBE2A", "UBE3A", "WDR45", "ZEB2"]
