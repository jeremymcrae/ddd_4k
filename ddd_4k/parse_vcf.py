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

def get_vcf_header(vcf):
    """ get the full VCF header as a list of lines
    
    Args:
        vcf: handle to vcf file
    
    Returns:
        list of header lines
    """
    
    current_pos = vcf.tell()
    vcf.seek(0)
    
    header = []
    for line in vcf:
        if not line.startswith("#"):
            break
        
        header.append(line)
    
    vcf.seek(current_pos)
    
    return(header)

def exclude_vcf_header(vcf):
    """ move the file handle to just beyond the VCF header lines
    
    Args:
        vcf: file handle to VCF
    """
    
    current_pos = vcf.tell()
    
    while vcf.readline().startswith("#"):
        current_pos = vcf.tell()
    
    vcf.seek(current_pos)

def parse_vcf_info(info):
    """ parse the info string from the VCF INFO field.
    
    Args:
        info: VCF INFO field for a single variant.
    
    Returns:
        dictionary of key, values from the info field. Info keys that are flags
        are given True values.
    """
    
    info_dict = {}
    for item in info.strip().split(";"):
        if "=" in item:
            item = item.split("=")
            key = item[0]
            value = item[1]
        else:
            key = item
            value = True
        
        info_dict[key] = value
    
    return info_dict

def parse_vcf_format(keys, values):
    """ get a dictionary of format values
    
    Args:
        keys: string of format keys from the format field
        values: string for the sample format values
    
    Returns:
        a dictionary of sample format values e.g.
        {"GT": "0/1", "DP": "10", "GQ": "20"}
    """
    
    return dict(zip(keys.split(":"), values.split(":")))
