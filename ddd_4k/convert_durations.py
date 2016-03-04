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

from datetime import timedelta

import pandas

def get_duration(value_and_age, cap_age=20):
    """ standardises a duration string into number of seconds
    
    Many children have the time from birth until they achieved several
    developmental milestones, such as smiled, sat up, walked and spoke their
    first words. These can be recorded as number of weeks, months, years, and
    can mix the units within each achievement. We need to standardise these to
    be able to analyse differences in groups for these achievements.
    
    Args:
        value: the string value for the duration. These can indicate a number of
            weeks, months, years etc. For example:
                weeks: "1 week", 6 weeks, "52 weeks or more"
                months: "7 months", "23 months"
                years: "2-2.5 years", "3-4 years", "5 years and over"
                missing values: "Unknown", "NA", "Not yet achieved"
        ages: decimal age of probands in years.
        cap_age: some probands have not yet achieved the milestone. Rather than
            losing them from the cohort, we give them their decimal age, so long
            as they are older than the normal age at which a milestone would be
            considered delayed. This age differs by phenotype, so the age for
            delayed social smiling would be > 3 months, but the age for delayed
            walking would be 18 months.
    
    Returns:
        the total number of seconds in the duration
    """
    
    value = value_and_age[0]
    age = value_and_age[1]
    
    weeks_per_month = 4.33
    weeks_per_year = 52.18
    
    if pandas.isnull(value):
        return None
    elif 'not yet achieved' in value.lower():
        if age > cap_age:
            value = "{} years".format(age)
        else:
            return None
    elif not any([ x in value for x in ["week", "month", "year"] ]):
        return None
    
    # get the numeric value from the string
    count = value.split(" ")[0]
    if "-" in value:
        count = count.split("-")[0]
    count = float(count)
    
    if "week" in value:
        duration = timedelta(weeks=count)
    elif "month" in value:
        duration = timedelta(weeks=count * weeks_per_month)
    elif "year" in value:
        duration = timedelta(weeks=count * weeks_per_year)
    
    return duration.total_seconds()
