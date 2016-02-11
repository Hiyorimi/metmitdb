#!/usr/bin/python
"""
==================================================================
Phylogenetic representativeness: a new method for evaluating taxon
sampling in evolutionary studies
==================================================================

CODENAME:     PhyRe
DESCRIPTION:

Copyright (c) 2009 Ronald R. Ferrucci, Federico Plazzi, and Marco Passamonti.
Rewritten and enhanced by Malev K, 2015

According to license limitations:

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

"""

from re import * #NOQM
import sys
from optparse import OptionParser
from random import sample
from collections import defaultdict, Counter


def PathLength(population, taxon, count_of_the_taxons_entries_on_each_level):
    """Calculates Path length for given master list. Results are used later."""

    number_of_unique_taxons_entries_on_each_level = {}
    taxon_names = defaultdict()
    unique_taxon_names = defaultdict()
    taxon_level_counter = defaultdict()

    for t in taxon:
        taxon_names[t] = []

    for entry in population:
        for t in taxon:
            taxon_names[t].append(population[entry][t])

    for t in taxon:
        unique_taxon_names[t] = frozenset(taxon_names[t])
        taxon_level_counter[t] = Counter(taxon_names[t])

    for t in taxon:
        count_of_the_taxons_entries_on_each_level[t] = {}

        if taxon.index(t) == 0:
            for i in unique_taxon_names[t]:
                # count number of taxon entries for each unique taxon in group
                count_of_the_taxons_entries_on_each_level[t][i] = taxon_level_counter[t][i]
        else:
            for i in unique_taxon_names[t]:
                if i not in unique_taxon_names[taxon[taxon.index(t) - 1]]:
                    # if this taxon doesn't belong to previous level, we get count of it
                    # it might get the name of the parent, if we have missing data
                    count_of_the_taxons_entries_on_each_level[t][i] = taxon_level_counter[t][i]

        number_of_unique_taxons_entries_on_each_level[t] = len(count_of_the_taxons_entries_on_each_level[t])

    print "Finished getting number of unique entries"

    n = [float(len(count_of_the_taxons_entries_on_each_level[t])) for t in taxon]
    # we insert N for the common ancestor of all entries
    n.insert(0, 1.0)

    unscaled_weights = []
    for i in range((len(n) - 1)):
        j = i + 1

        if n[i] > n[j]:
            c = 1
        else:
            c = (1 - n[i] / n[j])

        unscaled_weights.append(c)

    s = sum(unscaled_weights)
    adjusted_coefficient = [i * 100 / s for i in unscaled_weights]

    coef = {}
    pathLengths = {}
    for i in range(len(taxon)):
        t = taxon[i]
        coef[t] = sum(adjusted_coefficient[i:])
        pathLengths[t] = adjusted_coefficient[i]

    return coef, number_of_unique_taxons_entries_on_each_level, pathLengths


def ATDmean(data, sample, taxon, coef):
    """Calculates mean average taxonomic distinctness for given data"""

    N = len(sample)

    Taxon = {}
    taxonN = {}
    AvTD = 0
    n = 0
    # Taxon are counts of taxa at each level, taxonN are numbers of pairwise differences
    # at each level, with n being the accumlation of pairwise differences at that level. the difference
    # between n and TaxonN is the number of species that are in different taxa in that level
    # but not in upper levels

    for t in taxon:
        Taxon[t] = {}
        x = [data[i][t] for i in sample]
        for i in set(x):
            Taxon[t][i] = x.count(i)

    for t in taxon:
        taxonN[t] = sum([Taxon[t][i] * Taxon[t][j] for i in Taxon[t] for j in Taxon[t] if i != j])
        n = taxonN[t] - n
        AvTD = AvTD + (n * coef[t])
        n = taxonN[t]

    # print sample
    AvTD /= (N * (N - 1))

    return AvTD, taxonN, Taxon


def ATDvariance(taxonN, sample, atd, taxon, coef):
    """Calculates variance of average taxonomic distinctness for given data"""

    vtd = []

    vtd = 0
    N = 0
    n = 0

    for t in taxon:
        n = taxonN[t] - n
        vtd = vtd + n * coef[t] ** 2
        n = taxonN[t]

    N = len(sample)
    n = N * (N - 1)

    vtd = (vtd - ((atd * n) ** 2) / n) / n

    return vtd


def euler(data, atd, TaxonN, taxon, Taxon, coef):
    """Calculates von Euler's index of imbalance for given data"""

    sample = data.keys()

    n = len(sample)
    TDmin = 0
    N = 0
    for t in taxon:
        k = len(Taxon[t])
        TDmin += coef[t] * (((k - 1) * (n - k + 1) * 2 + (k - 1) * (k - 2)) - N)
        N += ((k - 1) * (n - k + 1) * 2 + (k - 1) * (k - 2)) - N

    TDmin /= (n * (n - 1))

    taxon.reverse()
    TaxMax = defaultdict()
    taxonN = defaultdict()

    for t in taxon:
        TaxMax[t] = []
        if taxon.index(t) == 0:
            TaxMax[t] = []
            for i in range(len(Taxon[t])):
                TaxMax[t].append([])
            for i in range(len(Taxon[t])):
                TaxMax[t][i] = [sample[j] for j in range(i, n, len(Taxon[t]))]
        else:
            TaxMax[t] = []
            for i in range(len(Taxon[t])):
                TaxMax[t].append([])
                s = taxon[taxon.index(t) - 1]

                Tax = [TaxMax[s][j] for j in range(i, len(Taxon[s]), len(Taxon[t]))]

                for j in Tax:
                    TaxMax[t][i] += j
        TaxMax[t].reverse()

    taxon.reverse()
    TDmax = 0
    n = 0
    N = len(sample)

    for t in taxon:
        taxonN[t] = sum(
                [len(TaxMax[t][i]) * len(TaxMax[t][j]) for i in range(len(TaxMax[t])) for j in range(len(TaxMax[t])) if
                 i != j])
        n = taxonN[t] - n
        TDmax += n * coef[t]
        n = taxonN[t]
        # for i in TaxMax[t]:
        #    print t, len(i)

    TDmax /= (N * (N - 1))

    EI = (TDmax - atd) / (TDmax - TDmin)

    return {'EI': EI, 'TDmin': TDmin, 'TDmax': TDmax}


def get_sample_subset_from_sample_file(samplefile, population):
    """Function, which parses parameters

     Parameters
    ----------
    samplefile : string
        Sample file filename
    population : dict
        Dict with master list

    Returns
    -------
    sample : dict
        subsample from master list
    """

    sample = {}

    for i in open(samplefile):

        if match('Taxon:', i):
            continue
        elif match('Coefficients:', i):
            continue

        x = i.split()

        species = x[0]
        try:
            sample[species] = population[species]
        except KeyError:
            print "Key Erorr: ", species
            return sample

    return sample


def printResults(taxon, taxonN, popN, pathLengths, results):
    """Prints result of performed analysis in the preset format

    Parameters
    ----------
    taxonN : dict
    taxonN : dict
    popN : dict
    pathLengths : dict
        Dictionary with taxon Path length data. Taxon name serves as key
    results : list
        List of computed coefficients

    Returns
    -------
    analysis_result : dict
        Dictionary with Average taxonomic distinctness, Variation in taxonomic distinctness,
        Minimum and maximum taxonomic distinctness and von Euler's index of imbalance

    """

    analysis_result = "Output from Average Taxonomic Distinctness\n\n"
    analysis_result += "Number of taxa and path lengths for each taxonomic level:\n"

    for t in taxon:
        analysis_result += '%-10s\t%d\t%.4f\n' % (t, popN[t], pathLengths[t])
        n = taxonN[t]

    analysis_result += "\n"

    for f in results:
        analysis_result += "---------------------------------------------------\n\n"
        analysis_result += "Results for sample: %s\n\n" % (f)
        analysis_result += "Dimension for this sample is %s\n\n" % (results[f]['n'])
        analysis_result += "Number of taxa and pairwise comparisons  at each taxon level:\n"

        n = 0
        for t in taxon:
            N = results[f]['N'][t] - n
            analysis_result += '%-10s\t%i\t%i\n' % (t, len(results[f]['taxon'][t]), N)
            n = results[f]['N'][t]

        analysis_result += """\nNumber of pairwise comparisons is for pairs that differ \
at each level excluding comparisons that differ at upper levels\n"""

        analysis_result += "Average taxonomic distinctness      = %.4f\n" % results[f]['atd']
        analysis_result += "Variation in taxonomic distinctness = %.4f\n" % results[f]['vtd']
        analysis_result += "Minimum taxonomic distinctness      = %.4f\n" % results[f]['euler']['TDmin']
        analysis_result += "Maximum taxonomic distinctness      = %.4f\n" % results[f]['euler']['TDmax']
        analysis_result += "von Euler's index of imbalance      = %.4f\n" % results[f]['euler']['EI']

    analysis_result += "\n\n---------------------------------------------------\n"

    return analysis_result


def print_funnel_data(p, d1, d2, population, taxon, coef):
    """Prints computed funnel plot data

    Parameters
    ----------
    permutations_number : int
        Permutations number
    d1 : int
    d2 : int
        d1 and d2 are range for number of species for funnel plot
    population : dict
        Dictionary with population information
    taxon : dict
        Dictionary with taxon data
    coef : list
        List of computed coefficients

    Returns
    -------
    analysis_result : str
        Results as a ready-to-print report

    """

    pop = population.keys()

    dims = []
    ciarray = []
    x = []
    carray = []

    analysis_result = """Confidence limits for average taxonomic distinctness and variation in taxonomic distinctness
    limits are lower 95% limit for AvTD and upper 95% limit for VarTD\n"""
    analysis_result += "Number of permutations for confidence limits = %s dimension\n\n" % (p)

    analysis_result += "AvTD05%   AvTDmean  AvTD95%   AvTDup    VarTDlow   VarTD05%   VarTDmean  VarTD95%\n\n"

    for d in range(d1, d2 + 1):

        x.append(d)
        AvTDci = []
        VarTDci = []
        for j in range(p):
            rsamp = sample(pop, d)

            atd, taxonN, Taxon = ATDmean(population, rsamp, taxon, coef)
            AvTDci.append(atd)
            vtd = ATDvariance(taxonN, rsamp, atd, taxon, coef)
            VarTDci.append(vtd)

        AvTDci.sort()
        VarTDci.sort()

        AvTD = AvTDci[int(.05 * p)], sum(AvTDci) / p, AvTDci[int(.95 * p)], max(AvTDci)
        VarTD = min(VarTDci), VarTDci[int(.05 * p)], sum(VarTDci) / p, VarTDci[int(.95 * p)]

        dims.append(d)
        ciarray.append(AvTD[0])
        carray.append(AvTD[1])

        analysis_result += '%i        %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n' \
                           % (d, AvTD[0], AvTD[1], AvTD[2], AvTD[3],
                              VarTD[0], VarTD[1], VarTD[2], VarTD[3])

    return analysis_result


def process_population_file_line(input_line, taxon, Index, missing='y'):
    """Processes line from population file,
    pushes read values into species_taxonomy dict and
    then returns them along with the species name

    Parameters
    ----------
    input_line : str
        Line from a population file, split by tabs or spaces
    taxon : dict
        Dictionary with taxons
    Index : dict
        Dict containing index information
    missing : str
        y if missing taxons are met in file

    Returns
    -------
    species_name : str
        Name of species (or least taxon)
    species_taxonomy : dict
        Dictionary with taxon as a key and taxon value as value

    """

    input_line.strip()
    splitted_line = input_line.split()

    species_name = splitted_line[0]
    species_taxonomy = {}

    if missing == 'y':
        mtax = ''
        for t in taxon:
            if splitted_line[Index[t]] == '/':
                temp = mtax
            else:
                temp = splitted_line[Index[t]]
                mtax = splitted_line[Index[t]]

            species_taxonomy[t] = temp

    else:
        for t in taxon:
            species_taxonomy[t] = splitted_line[Index[t]]

    return species_name, species_taxonomy


def phy_re_analysis(options, args):
    """
    Script should be launched as:

    python PhyRe.py [samplefile] [masterlistfile] s1 s2 [options]

    Parameters
    ----------

    p : int
        permutations for confidence intervals
    d1 : int
        d1 and d2 are range for number of species for funnel plot
    d2 : int
        d1 and d2 are range for number of species for funnel plot


    Returns
    -------
    first_file.out : file
        Results from analyses of the sample. By default,
        the output file has the same name of the sample file
        with extension .OUT
    second_file.out : file
        Results from random subsamples of the master list.
        The funnel output file has the same name with suffix
        "_funnel" and extension .OUT. I


    Notes
    ----------
    Described in Phylogenetic representativeness: a new
    method for evaluating taxon sampling in evolutionary studies [1]

     References
    ----------
    .. [1] http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-209
    .. [2] http://www.mozoolab.net/downloads/manual.pdf
    """

    samplefile = args['samplefile']
    popfile = args['popfile']
    d1 = args['d1']
    d2 = args['d2']

    output_as_string = False

    if options['m']:
        missing = options['m']
    else:
        missing = 'n'

    if options['o']:
        out = options['o']
    else:
        out = samplefile.split('.')[0]

    if options['p']:
        p = options['p']
    else:
        p = 1000

    if options['c']:
        ci = options['c']
    else:
        ci = 'y'

    if options['b']:
        batch = options['b']
    else:
        batch = 'n'

    if options['l']:
        pathlengths = options['l']
    else:
        pathlengths = 'n'

    if options['s']:
        output_as_string = True

    sample = defaultdict()
    # Population - dictionary with population file information
    population = defaultdict()

    if batch == 'y':
        Files = []
    else:
        Files = [samplefile]

    Index = {}
    Taxon = defaultdict()
    coef = {}
    taxon = []
    pathLengths = defaultdict()

    for i in open(samplefile):

        if batch == 'y':
            j = i.strip()
            Files.append(j)
        else:
            break

    duplicates = []

    with open(popfile) as fp:
        population_file_entries = fp.readlines()

    # We reed two first lines to check if we have System information
    # If we have it, we rebuild the array of lines
    lines_removal_counter = 0

    for i in population_file_entries[:2]:
        # If we encounter string, starting with "Taxon" we get
        # information about taxons
        if match('Taxon:', i):
            lines_removal_counter += 1
            x = i.split()
            x.remove('Taxon:')

            for i in x:
                taxon.append(i)
                j = x.index(i)
                # Index list is used to get the indexation of
                # taxon during string parsing
                Index[i] = j + 1
            continue

        elif match('Coefficients:', i):
            lines_removal_counter += 1
            x = i.split()
            x.remove('Coefficients:')
            x = map(eval, x)

            for t in taxon:
                i = taxon.index(t)
                coef[t] = sum(x[i:])
                pathLengths[t] = x[i]

            continue

    population_file_entries = population_file_entries[lines_removal_counter:]

    # opening population file and getting information
    for i in population_file_entries:

        # here starts entry processing
        (species_name, species_taxonomy) = \
            process_population_file_line(i, taxon, Index, missing)

        if species_name in population:
            duplicates.append(species_name)
        else:
            population[species_name] = species_taxonomy

    sample = population.copy()

    if len(duplicates) > 0:
        print "Population master list contains %s duplicates" \
              % (len(duplicates))

    if pathlengths == 'n':
        coef, popN, pathLengths = PathLength(population, taxon, Taxon)
    if pathlengths == 'y':
        XXX, popN, YYY = PathLength(population, taxon, Taxon)

    print "Finished path length calculation"

    results = {}

    """Opening all sample files (or the sample file)
     and getting information about it"""
    for f in Files:
        sample = get_sample_subset_from_sample_file(f, population)
        f = f.split('.')
        f = f[0]

        results[f] = {}

        samp = sample.keys()

        atd, taxonN, Taxon = ATDmean(sample, samp, taxon, coef)
        average_taxonomic_distinctness_variance = \
            ATDvariance(taxonN, samp, atd, taxon, coef)
        euler_results = euler(sample, atd, taxonN, taxon, Taxon, coef)

        results[f]['atd'] = atd
        results[f]['vtd'] = average_taxonomic_distinctness_variance
        results[f]['euler'] = euler_results
        results[f]['N'] = taxonN
        results[f]['n'] = len(sample)
        results[f]['taxon'] = Taxon

    phy_re_result = printResults(taxon, taxonN, popN, pathLengths, results)

    funnel_data = ''
    if ci == 'y':
        funnel_data = print_funnel_data(p, d1, d2, population, taxon, coef)

    if output_as_string:
        return phy_re_result, funnel_data
    else:
        with open(out + '.out', 'w') as fp:
            fp.write(phy_re_result)
        with open(out.split('_')[0] + '_funnel.out', 'w') as fp:
            fp.write(funnel_data)


if __name__ == "__main__":
    samplefile = sys.argv[1]
    del sys.argv[1]
    popfile = sys.argv[1]
    del sys.argv[1]

    p = 1000
    d1 = 10
    d2 = 70
    ci = 'y'
    b = 'n'
    l = 'n'
    batch = b
    pathlengths = l
    missing = 'n'

    parser = OptionParser()

    d1 = int(sys.argv[1])
    del sys.argv[1]
    d2 = int(sys.argv[1])
    del sys.argv[1]

    parser.add_option('-o')
    parser.add_option('-p', type='int')
    parser.add_option('-c')
    parser.add_option('-b')
    parser.add_option('-l')
    parser.add_option('-m')

    (options, _) = parser.parse_args()

    args = dict()

    args['samplefile'] = samplefile
    args['popfile'] = popfile
    args['d1'] = d1
    args['d2'] = d2

    options = vars(options)
    options['s'] = None

    phy_re_analysis(options, args)
