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

import json

from denovonear.ensembl_requester import EnsemblRequest

class EnsemblVariant(EnsemblRequest):
    """ this class extends the EnsemblRequest class from denovonear, in order to
    extract HGVS codes for variants from the ensembl REST API.
    """
    
    def parse_alternate_consequences(self, variant, gene, vep_cq):
        """ at the moment we're looking for the HGVS nomenclature for the
        canonical transcript, where that transcript needs to match our most
        severe consequence.
        
        TODO: I might need to shift this to using HGVS nomenclature for the
        longest protein_coding transcript, where that transcript consequence
        matches our most severe consrquence prediction.
        
        Args:
            variant: dictionary of possible consequences or variant, obtained
                from the Ensembl REST API.
            gene: HGNC symbol for the variant that we need to match against.
            vep_cq: VEP consequence that we need to match against.
        
        Returns:
            dictionary entry for the canonical transcript for the correct gene
            and VEP consequence. Returns None if not found.
        """
        
        for cq in variant["transcript_consequences"]:
            if "gene_symbol" not in cq:
                continue
            if cq["gene_symbol"] != gene:
                continue
            if cq["biotype"] not in ["protein_coding", "polymorphic_pseudogene"]:
                continue
            if vep_cq not in cq["consequence_terms"]:
                continue
            
            if "canonical" in cq:
                return cq
    
    def get_variant_hgvs(self, chrom, pos, ref, alt, gene, vep_cq):
        """ extract the HGVS predictions for a variant
        """
        
        self.request_attempts = 0
        end_pos = pos + len(ref) - 1
        
        headers = {"Content-Type": "application/json"}
        coordinate = "{}:{}:{}".format(chrom, pos, end_pos)
        ext = "/vep/human/region/{}/{}?hgvs=1;canonical=1".format(coordinate, alt)
        
        r = self.ensembl_request(ext, coordinate, headers)
        variant = json.loads(r)[0]
        
        return self.parse_alternate_consequences(variant, gene, vep_cq)
