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

import argparse

import pandas

from ddd_4k.constants import DENOVO_PATH, VALIDATIONS

OUTPUT = "/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.high_quality.txt"

def get_options():
    """ parse the command line arguments
    """
    
    parser = argparse.ArgumentParser(description="script to check whether"
        "candidate de novo CNVs overlap candidate novel genes")
    parser.add_argument("--de-novos", default=DENOVO_PATH, \
        help="Path to table of de novo variants.")
    parser.add_argument("--validations", default=VALIDATIONS, \
        help="Path to table of validation results.")
    parser.add_argument("--output", default="high_quality_de_novos.txt", \
        help="Path to write table of high quality de novo variants to.")
    
    args = parser.parse_args()
    
    return args

def main():
    args = get_options()
    
    de_novos = pandas.read_table(args.denovos, sep="\t")
    validations = pandas.read_table(args.validations, sep="\t")

    de_novos = de_novos.merge(validations[["person_id", "chrom", "start_pos", "status"],
        how="left",
        on_left=["person_stable_id", "chrom", "pos"],
        on_right=["person_id", "chrom", "start_pos"])

    # set the de_novos that passed validation to pp_dnm scores of 1, i.e. high
    # posterior probability that these candidates are truly de novo mutations
    de_novos["pp_dnm"][(de_novos["status"] == "de_novo") & (~de_novos["status"].isnull())] = 1.000

    de_novos = de_novos[(de_novos["pp_dnm"] > 0.9) & ~(de_novos$["pp_dnm"].isnull())]

    de_novos.to_csv(args.output, sep="\t", index=False)

if __name__ == '__main__':
    main()
