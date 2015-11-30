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

import sys
import json

IS_PYTHON3 = sys.version_info[0] == 3
 
if IS_PYTHON3:
    import urllib.request as request
else:
    import urllib2 as request

def open_url(url, headers):
    """ open url with python libraries
    """
    
    req = request.Request(url, headers=headers)
    
    try:
        handler = request.urlopen(req)
    except Exception as e:
        handler = e
    
    status_code = handler.getcode()
    response = handler.read()
    if IS_PYTHON3:
        response = response.decode("utf-8")
    
    # parse the headers into a key, value dictionary
    headers = {}
    for key, value in zip(handler.headers.keys(), handler.headers.values()):
        headers[key.lower()] = value
    
    return response, status_code, headers

def doi_to_pubmed(doi, email):
    """ convert a doi to a pubmed ID, using the NCBI DOI converter REST API
    
    Args:
        doi: digitical object identifier for an article e.g. "10.1093/nar/gks1195"
        email: email address of the person using this tool (for the NCBI server
            records).
    
    Returns:
        The pubmed ID for the article.
    """
    
    ncbi = "http://www.ncbi.nlm.nih.gov"
    tool = "gene_reporter"
    doi_url = "{}/pmc/utils/idconv/v1.0/?tool={}&email={}".format(ncbi, tool, email)
    
    url = "{}{}&ids={}&format=json".format(doi_url, email, doi)
    
    headers = {}
    response, status_code, headers = open_url(url, headers)
    
    if status_code != 200:
        sys.exit("failed to access the NCBI doi converter API at {}, returned status code: {}, with response:".format(url, status_code, response))
    
    response = json.loads(response)
    records = response["records"]
    
    # What if there are multiple pubmed IDs for an article? I'll handle this
    # when it happens
    assert len(records) == 1
    
    return records[0]["pmid"]
