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

DENOVO_PATH = "{}/de_novos.ddd_4k.ddd_only.2015-09-02.txt".format(user_dir)
VALIDATIONS = "{}/de_novos.validation_results.2015-09-22.txt".format(user_dir)
KNOWN_GENES = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/ddg2p_freeze_Jul15_corrected2_with_gencode_v19_coordinates_fixed.txt"
PHENOTYPES = "{}/phenotypes_and_patient_info.txt".format(datafreeze)
SANGER_IDS = "{}/person_sanger_decipher.txt".format(datafreeze)
FAMILIES = "{}/family_relationships.txt".format(datafreeze)
TRIOS = "{}/trios.txt".format(datafreeze)
DATATYPES = "{}/person_datatypes.txt".format(datafreeze)
