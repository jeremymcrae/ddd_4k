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

import re

def find_hpo_matches(hpo_phenotypes, required):
    """ for a list of HPO term list, find if any are in a set of required terms.
    
    This function should work on lists of terms from the known genes (entries
    separated by " ;"), as well as terms from the phenotypes (entries separated
    by "|").
    
    Args:
        hpo_phenotypes: pandas series of HPO strings per gene
        required: set of HPO terms that we wish to check for.
    
    Returns:
        list of True/False entries indexed to the hpo phenotype entries
    """
    
    matches = []
    for hpo in hpo_phenotypes:
        if type(hpo) != str:
            matches.append(False)
        else:
            hpo = set(re.split(" ;|\|", hpo))
            matches.append(len(hpo & required) > 0)
    
    return matches
