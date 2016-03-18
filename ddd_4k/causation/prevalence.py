
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

from __future__ import division

import itertools

import pandas
import numpy
from scipy.stats import gaussian_kde

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, gridspec
import seaborn

seaborn.set_context("notebook", font_scale=2)
seaborn.set_style("white", {"ytick.major.size": 10, "xtick.major.size": 10})

def get_birth_prevalence(cohort_n, excess, cnv_yield=0.1, missing_variants=119,
    enrichment=120.0):
    """ calculate the prevalence of live births with developmental disorders
    caused by de novo mutations.
    
    Args:
        cohort_n: size of the studied cohort
        excess: dictionary with numbers of excess mutations for each of
            synonymous, missense and loss-of-function mutations
        cnv_yield: estimated diagnostic yield from de novo CNVs in probands with
            developmental disorders, obtained from published estimates
            (doi:10.1038/ng.909, doi:10.1097/GIM.0b013e318194ee8f and
            doi:10.1016/j.ajhg.2010.04.006)
        missing_variants: number of variants missing from the cohort due to ease
             of clinical recognisability.
        enrichment: factor by whcih the cohort is enriched for loss-of-function
            de novo mutations in known haploinsufficient genes with low clinical
            recognisability.
    
    Returns:
        estimate of the prevalence of live births with developmental disorders
        caused by de novo mutations
    """
    
    cohort_yield = (excess["loss-of-function"]["excess"] + excess["missense"]["excess"])/cohort_n
    missing_yield = missing_variants/cohort_n + cnv_yield
    
    # adjust the yield in the DDD for the proportions that exist outside the
    # DDD, since these expand the DDD cohort size.
    cohort_yield = cohort_yield/(1 + missing_yield)
    
    return (cohort_yield + missing_yield)/enrichment

def plot_prevalence_by_age(prevalence, phenotypes, diagnosed, uk_ages,
    dad_rate=1.53, mom_rate=0.86):
    """ plot the prevalence of developmental disorders from de novo mutations
    
    This function generates three plots in one figure. The first
    
    Args:
        prevalence: estimated birth prevalence for children with developmental
            disorders caused by de novo mutations
        phenotypes: pandas DataFrame of phenotypic data for probands within the
            trios in the cohort.
        mutations_per_year: estimate for the number of additional mutations
            acquired due to older fathers, as additional mutations per year.
    """
    
    fig = pyplot.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[4, 1])
    
    prevalence_axes = pyplot.subplot(gs[0, 0])
    dad_axes = pyplot.subplot(gs[1, 0], sharex=prevalence_axes)
    mom_axes = pyplot.subplot(gs[0, 1], sharey=prevalence_axes)
    
    low_age = 20
    high_age = 45
    
    plot_birth_prevalences(prevalence_axes, prevalence, phenotypes,
        dad_rate, mom_rate, lower_age=low_age, upper_age=high_age)
    plot_paternal_ages(dad_axes, phenotypes, uk_ages, low_age, high_age)
    plot_maternal_ages(mom_axes, phenotypes, uk_ages, low_age, high_age)
    
    e = prevalence_axes.set_xlim((low_age - 2, high_age + 2))
    e = prevalence_axes.set_ylim((low_age - 2, high_age + 2))
    e = prevalence_axes.invert_yaxis()
    
    e = dad_axes.legend(fontsize='medium')
    e = dad_axes.invert_yaxis()
    
    e = mom_axes.legend(fontsize='medium')
    
    fig.savefig("results/prevalence_by_age.smaller_increment.pdf", format="pdf",
        bbox_inches='tight', pad_inches=0, transparent=True)
    pyplot.close()

def plot_birth_prevalences(ax, prevalence, phenotypes, dad_rate, mom_rate,
    lower_age=20, upper_age=40):
    """plot the birth prevalence of developmental disorders from de novo mutations
    
    We plot a point estimate of birth prevalence of children with developmental
    disorders caused by de novo mutations, at the median age of our cohort. We
    extrapolate the birth prevalence at lower and higher ages from the number of
    additional mutations acquired from fathers at older ages (per year). This
    requires knowning how many mutations on average a father would have at the
    median age, and how many of those would cause a developmental disorder. We
    can estimate this from the birth prevalence divided by the median number of
    mutations.
         
    Args:
        ax: matplotlib pyplot axes object for the plot.
        prevalence: estimated birth prevalence for children with developmental
            disorders caused by de novo mutations
        phenotypes: pandas DataFrame of phenotypic data for probands within the
            trios in the cohort.
        dad_rate: estimate for the number of additional mutations
            acquired in older fathers per year.
        mom_rate: estimate for the number of additional mutations
            acquired in older mothers per year.
        lower_age: lower age limit to estimate birth prevalence at.
        upper_age: upper age limit to estimate birth prevalance at.
    """
    
    median_dad_age = numpy.median(phenotypes["fathers_age"])
    median_mom_age = numpy.median(phenotypes["mothers_age"])
    
    # estimate the median number of mutations each child obtains at the
    # median paternal and maternal ages.
    median_mutations = median_dad_age * dad_rate + median_mom_age * mom_rate
    
    increment = 3.0
    ages = range(lower_age, upper_age + int(increment), int(increment))
    ages = [ x + increment/2 for x in ages ]
    prevalences = []
    for mom_age, dad_age in itertools.product(ages, repeat=2):
        mutations = dad_age * dad_rate + mom_age * mom_rate
        aged_prevalence = mutations/median_mutations * prevalence * 100
        prevalences.append(aged_prevalence)
    
    # reshape the list to a n x n array, then insert into a dataframe
    prevalences = numpy.reshape(prevalences, (len(ages), len(ages)))
    ages = [ x - increment/2 for x in ages ]
    
    green_blue1 = get_colormap()
    mesh = ax.pcolormesh(numpy.array(ages), numpy.array(ages), prevalences,
        cmap=green_blue1)
    annotate_heatmap(ax, mesh, ages, increment)
    e = pyplot.colorbar(mesh, ticks=[0.30, 0.35, 0.40, 0.45, 0.50, 0.55])
    
    # scale_heatmap_by_parental_proportion(ax, prevalence, phenotypes, dad_rate,
    #     mom_rate, lower_age=20, upper_age=40)
    
    e = ax.spines['bottom'].set_visible(False)
    e = ax.spines['right'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('top')
    e = ax.yaxis.set_ticks_position('left')
    
    e = ax.set_xlabel("Paternal age (years)")
    e = ax.set_ylabel("Maternal age (years)")
    e = ax.xaxis.set_label_position('top')

def get_colormap():
    """ define a custom colormap for shading the prevalence plot
    
    This colormap goes green -> yellow -> red. The yellow is at the midpoint.
    
    Returns:
        matplotlib linear colormap.
    """
    
    cdict = {'red':   ((0.0,  0.118, 0.118),
                       (0.5,  0.902, 0.902),
                       (1.0,  0.902, 0.902)),

             'green': ((0.0, 0.588, 0.588),
                       (0.5, 0.902, 0.902),
                       (1.0, 0.118, 0.118)),

             'blue':  ((0.0, 0.392, 0.392),
                       (0.5, 0.235, 0.235),
                       (1.0, 0.118, 0.118))}

    return matplotlib.colors.LinearSegmentedColormap('GreenBlue1', cdict)

def relative_luminance(color):
    """Calculate the relative luminance of a color according to W3C standards
    Parameters
    ----------
    color : matplotlib color or sequence of matplotlib colors
        Hex code, rgb-tuple, or html color name.
    Returns
    -------
    luminance : float(s) between 0 and 1
    """
    rgb = matplotlib.colors.colorConverter.to_rgba_array(color)[:, :3]
    rgb = numpy.where(rgb <= .03928, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    lum = rgb.dot([.2126, .7152, .0722])
    try:
        return lum.item()
    except ValueError:
        return lum

def annotate_heatmap(ax, mesh, ages, increment):
    """Add textual labels with the value in each cell."""
    
    ages = [ x + increment/2 for x in ages[:-1] ]
    mesh.update_scalarmappable()
    xpos, ypos = numpy.meshgrid(ages, ages)
    for x, y, val, color in zip(xpos.flat, ypos.flat,
                                mesh.get_array(), mesh.get_facecolors()):
        if val is not numpy.ma.masked:
            l = relative_luminance(color)
            text_color = ".0" if l > .408 else "w"
            val = "{:.2f}".format(val)
            text_kwargs = dict(color=text_color, ha="center", va="center")
            ax.text(x, y, val, **text_kwargs)

def scale_heatmap_by_parental_proportion(ax, prevalence, phenotypes, dad_rate, mom_rate,
        lower_age=20, upper_age=40):
    """ heatmap of parental age by birth prevalence of developmental disorders
    from de novo mutations, where the heatmap points size is scaled by the
    proportion of the population within each age bracket. This is an alternate
    to the standard heatmap.
    
    Args:
        ax: matplotlib pyplot axes object for the plot.
        prevalence: estimated birth prevalence for children with developmental
            disorders caused by de novo mutations
        phenotypes: pandas DataFrame of phenotypic data for probands within the
            trios in the cohort.
        dad_rate: estimate for the number of additional mutations
            acquired in older fathers per year.
        mom_rate: estimate for the number of additional mutations
            acquired in older mothers per year.
        lower_age: lower age limit to estimate birth prevalence at.
        upper_age: upper age limit to estimate birth prevalance at.
    """
    
    median_dad_age = numpy.median(phenotypes["fathers_age"])
    median_mom_age = numpy.median(phenotypes["mothers_age"])
    
    # estimate the median number of mutations each child obtains at the
    # median paternal and maternal ages.
    median_mutations = median_dad_age * dad_rate + median_mom_age * mom_rate
    
    increment = 5.0
    ages = range(lower_age, upper_age + int(increment), int(increment))
    ages = [ x + increment/2 for x in ages ]
    data = pandas.DataFrame(columns=["dad_age", 'mom_age', 'prevalence', 'freq'])
    for mom_age, dad_age in itertools.product(ages, repeat=2):
        mutations = dad_age * dad_rate + mom_age * mom_rate
        aged_prevalence = mutations/median_mutations * prevalence * 100
        
        half = increment/2
        in_range = ((dad_age - half < phenotypes["fathers_age"]) &
                        (phenotypes["fathers_age"] < dad_age + half)) & \
            ((mom_age - half < phenotypes["mothers_age"]) &
                (phenotypes["mothers_age"] < mom_age + half))
        freq = sum(in_range)/len(phenotypes)
        data = data.append({'mom_age': mom_age, 'dad_age': dad_age,
            'prevalence': aged_prevalence, 'freq': freq}, ignore_index=True)
    
    data['freq'] = 200.0 + 1000 * \
        (data['freq'] - min(data['freq']))/(max(data['freq']) - min(data['freq']))
    data['prevalence'] = 0.9 - 0.9 * \
        (data['prevalence'] - min(data['prevalence']))/(max(data['prevalence']) - min(data['prevalence']))
    
    data['freq'] = data['freq'].astype(numpy.float64)
    data['prevalence'] = data['prevalence'].astype(str)
    
    e = ax.scatter(data['dad_age'], data['mom_age'], s=data['freq'],
        c=data['prevalence'])

def plot_paternal_ages(ax, phenotypes, uk_ages, lower_age=20, upper_age=40):
    """ plot the age distribution for fathers in the DDD cohort, along with
    the proportion with a child with a diagnostic de novo across age bins.
    
    Args:
        ax: matplotlib pyplot axes object for the plot.
        phenotypes: pandas DataFrame of phenotypic data for probands within the
            trios in the cohort, including a column for "fathers_age".
        lower_age: lower age limit to estimate birth prevalence at.
        upper_age: upper age limit to estimate birth prevalance at.
    """
    
    plot_uk_age_distribution(ax, uk_ages, 'fathers_count', lower_age, upper_age)
    
    ages = phenotypes["fathers_age"]
    ages = ages[~ages.isnull()]
    
    density = gaussian_kde(ages)
    x = numpy.arange(lower_age, upper_age, 0.1)
    e = ax.plot(x, density(x), label="DDD")
    e = ax.fill_between(x, density(x), alpha=0.5)
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['top'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    
    e = ax.yaxis.set_ticks_position('left')
    e = ax.set_ylabel("Density")
    e = ax.xaxis.set_ticks_position('none')

def plot_maternal_ages(ax, phenotypes, uk_ages, lower_age=20, upper_age=40):
    """ plot the age distribution for fathers in the DDD cohort, along with
    the proportion with a child with a diagnostic de novo across age bins.
    
    Args:
        ax: matplotlib pyplot axes object for the plot.
        phenotypes: pandas DataFrame of phenotypic data for probands within the
            trios in the cohort, including a column for "fathers_age".
        lower_age: lower age limit to estimate birth prevalence at.
        upper_age: upper age limit to estimate birth prevalance at.
    """
    
    plot_uk_age_distribution(ax, uk_ages, 'mothers_count', lower_age, upper_age)
    
    ages = phenotypes["mothers_age"]
    ages = ages[~ages.isnull()]
    
    density = gaussian_kde(ages)
    x = numpy.arange(lower_age, upper_age, 0.1)
    e = ax.plot(density(x), x, label="DDD")
    e = ax.fill_betweenx(x, density(x), alpha=0.5)
    
    e = ax.spines['right'].set_visible(False)
    e = ax.spines['bottom'].set_visible(False)
    e = ax.spines['left'].set_visible(False)
    
    e = ax.xaxis.set_ticks_position('top')
    e = ax.yaxis.set_ticks_position('none')
    e = ax.set_xlabel("Density")
    e = ax.xaxis.set_label_position('top')

def plot_uk_age_distribution(ax, uk_ages, column, lower_age=20, upper_age=40):
    """ plot the age distribution of fathers in the UK, to compare to our cohort
    
    Args:
        ax: matplotlib pyplot axes object for the plot.
        uk_ages: pandas DataFrame of ages and numbers of parents at each age for
            mothers and fathers.
        column: string for 'fathers_age' or 'mothers_age'
        lower_age: lower age limit to estimate birth prevalence at.
        upper_age: upper age limit to estimate birth prevalance at.
    """
    
    # we have a DataFrame containing ages and counts. We need to convert that to
    # a list of ages, where each age is repeated as many times as the count, in
    # order to use in a kernel density plot. Rather than have the ages jumping
    # from integer to integer, we add a bit of noise, so that the people in each
    # age are uniformly spread between that age and the succeeding age, e.g from
    # 20.0 to 20.99999. This avoids lumpiness in the kernel density estimates.
    # Another way to get around the lumpiness would be to adjust the kernel
    # density covariance factor, such as:
    #    # gaussian_kde.covariance_factor = lambda x: 0.095
    ages = []
    for k, x in uk_ages.iterrows():
        ages += list(x["age"] + numpy.random.uniform(size=x[column]))
    
    density = gaussian_kde(ages)
    x = numpy.arange(lower_age, upper_age, 0.1)
    y = density(x)
    
    # swap axes for maternal age distribution
    if column == 'mothers_count':
        y, x = x, y
    
    e = ax.plot(x, y, color="gray", label="UK")
    
    if column == 'mothers_count':
        e = ax.fill_betweenx(y, x, color='gray', alpha=0.5)
    else:
        e = ax.fill_between(x, y, color='gray', alpha=0.5)
    
