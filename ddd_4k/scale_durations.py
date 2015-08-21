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

def autoscale_durations(durations):
    """ convert the durations (in seconds) to a better unit
    
    Args:
        durations: list/pandas Series of durations in seconds
    
    Returns:
        The time in seconds converted to the better represenation, along with
        the name of the unit ("hour", "day" etc).
    """
    
    # Define the number of seconds in each time interval (note that the number
    # of seconds in a month doesn't take different month lengths into account)
    # and seconds in a year is the average when including leap years.
    hour = 3600
    day = hour * 24
    week = day * 7
    month = day * 30.4375
    year = day * 365.25
    intervals = {"hour": hour, "day": day, "week": week, "month": month, "year": year}
    
    delta = max(durations) - min(durations)
    
    for unit in ["hour", "day", "week", "month", "year"]:
        if 30 > (delta / intervals[unit]) > 3:
            break
    
    durations = [ x/intervals[unit] for x in durations ]
    
    return durations, unit
