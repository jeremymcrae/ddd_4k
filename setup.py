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

from setuptools import setup

setup(
    name = "ddd_4k",
    version = "0.1.0",
    author = "Jeremy McRae",
    author_email = "jeremy.mcrae@sanger.ac.uk",
    description = ("Analysis of de novos in the DDD 4k trios."),
    license = "MIT",
    packages=["ddd_4k"],
    install_requires=['pandas >= 0.13.1',
                      'statsmodels >= 0.5.0',
                      'matplotlib >= 1.3.1',
                      'seaborn >= 0.6.0',
                      'psycopg2 >= 2.6.0',
                      'denovonear >= 0.1.1',
                      'inflect >= 0.2.5',
                      'mupit >= 0.3.0'
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
    ]
)
