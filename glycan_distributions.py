#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from glycopeptidepy import PeptideSequence
from glycopeptidepy.io import fasta
from csv import DictReader
from requests import get
from os.path import basename
from collections import OrderedDict
from ast import literal_eval
from math import sqrt
from math import exp
from os import makedirs
from os.path import exists
import numpy as np
import re
import json
import csv
import argparse

#global variables
monosaccharide_order = ('Fuc', 'Hex', 'HexNAc', 'Neu5Ac')
glycan_order = ('HexNAc', 'Hex', 'Fuc', 'Neu5Ac')

# classes


class SiteSpecificGlycan(object):
    """An object that holds a glycan composition for a particular site,
    with additional information on protein accession id, glycan total signal,
    replicates in which the glycan is found, and the total number of
    observations of the glycan."""

    def __init__(self, glycan, glycosite, accession_id, total_signal, replicate_ids, ms2_score,
                 num_observations=1, averaged=False):
        """Creates a SiteSpecificGlycan object with inputs as attributes.
        num_observations defaults to 1. If glycan parsed to OrderedDict with
        keys in order Fuc,Hex,HexNAc,Neu5Ac. No return value."""

        # set all attributes by glycan
        self.glycosite = int(glycosite)
        self.accession_id = accession_id
        self.total_signal = float(total_signal)
        self.replicate_ids = set(replicate_ids)
        self.num_observations = int(num_observations)
        self.averaged = averaged
        self.glycan = glycan

    def increment_observations(self, increment=1):
        """Increment the number of observations of this glycan by increment.
        Returns the new number of observations. Increment size defaults to 1."""
        self.num_observations += increment
        return self.num_observations

    def add_replicate(self, replicate):
        """Add replicate to set existing replicates. Returns the new set of replicates."""
        self.replicate_ids.add(replicate)
        return self.replicate_ids

    def update_signal(self, signal):
        """Increase total_signal by the provided amount. Returns the new total_signal value."""
        self.total_signal += signal
        return self.total_signal

    def average(self):
        """'Averages' glycan by dividing its total signal by its number of observations.
        Can only be used if self.averaged is False. Changes self.averaged to True and
        returns the newly averaged total_signal."""
        assert self.averaged == False, 'Cannot average glycan ' + \
            self.glycan + '; has already been averaged'

        self.total_signal = self.total_signal / self.num_observations
        self.averaged = True
        return self.total_signal

    def __repr__(self):
        """Returns the string representation of this object (not pretty)."""
        return 'SiteSpecificGlycan ({}, {})'.format(self.glycan, self.total_signal)


class PairedGlycan(object):
    """Stores paired values for glycan observed in an expected and experimental sample.
    Information includes glycan compostion in OrderedDict form, site numbers, and total signals."""

    def __init__(self, glycan, expected_site, experimental_site, expected_signal, experimental_signal):
        """Creates a PairedGlycan object with inputs as attributes. Glycan must be in OrderedDict form."""
        assert type(glycan) == OrderedDict, 'glycan is ' + \
            str(type(glycan)) + ', not OrderedDict'
        self.glycan = glycan
        self.expected_site = expected_site
        self.experimental_site = experimental_signal
        self.expected_signal = expected_signal
        self.experimental_signal = experimental_signal

    def __repr__(self):
        """"""
        return 'PairedGlycan({} {} {})'.format(self.glycan, self.expected_signal, self.experimental_signal)

# functions


def cosine_similarity(v1, v2):
    """Returns the cosine similarity values for two vectors, which is simply the
    cosine value of the angle between the vectors, where
    cos(x) = (v1 * v2) / (||v1|| ||v2||). The vectors must be of equal length and
    cannot equal the zero vector."""

    # check v1 and v2 are valid
    assert len(v1) == len(v2), 'Vectors must must be the same length; len(v1) = ' +\
        str(len(v1)) + ' != ' + str(len(v2)) + ' len(v2)'
    assert set(v1) != set([0]), 'v1 = ' + str(v1) + ' which is the zero vector'
    assert set(v2) != set([0]), 'v2 = ' + str(v2) + ' which is the zero vector'

    # get lengths and dot product
    lenx = sqrt(sum([i**2 for i in v1]))
    leny = sqrt(sum([i**2 for i in v2]))
    dot = sum([v1[i] * v2[i] for i in range(len(v1))])

    # return cosine similarity
    return dot / (lenx * leny)


def kullback_leibler_divergence(p, q):
    """Measures how one distribution (q) diverges from a second distribution (p).
    p and q must be probability distributions. Think of p as the control distribution
    and q as the experimental distribution."""
    # make sure vectors are normalized
    assert len(p) == len(q)
    assert round(
        sum(p), 4) == 1, 'p is not a probability distribution; sum is ' + round(sum(p), 4)
    assert round(
        sum(q), 4) == 1, 'q is not a probabiltiy distribution; sum is ' + round(sum(q), 4)

    # set all zeroes to 1e-3 and then re-normalize vectors
    p = np.array(p)
    p[p == 0] = 1e-3
    p = p / p.sum()

    q = np.array(q)
    q[q == 0] = 1e-3
    q = q / q.sum()

    assert round(
        sum(p), 4) == 1, 'p is not a probability distribution; sum is ' + round(sum(p), 4)
    assert round(
        sum(q), 4) == 1, 'q is not a probabiltiy distribution; sum is ' + round(sum(q), 4)

    logdifs = np.log(p) - np.log(q)
    return np.dot(p, logdifs)


def normalized_kullback_leibler_divergence(p, q):
    """Calculates kl divergence for probability distributions p and q before
    returning normalized value"""
    div = kullback_leibler_divergence(p, q)
    return 1 - exp(-div)


def get_sequoninfodict_from_files_accid(csvfilelist, accessionid=None, seq='', score_cutoff=30, replicate_cutoff=2, sequon_length=4):
    """Takes list of DictReader CSV file objects produced from GlycReSoft
    glycopeptide-identification. Requires either a protein sequence or UniProt accession
    ID; if both provided, defaults to sequence. Returns dictionary mapping glycosylation
    sites of protein associated with accessionid to sequon and list of averaged
    SiteSpecificGlycan objects. Glycans with an ms2 score less than score_cutoff and
    observed in fewer replicates than replicate_cutoff will not be included in the
    returned siteinfodict. Sequon length determined by sequon_length. Dictionary returned
    is of the form:
    {
        SITE: {
            "sequon": "XXXX",
            "glycans": [ glycan1, glycan2, ... ]
        },
        ...
    }"""

    # check number of csvs is not less than the number of required replicates
    assert replicate_cutoff <= len(csvfilelist), \
        'Replicate cutoff = ' + \
        str(replicate_cutoff) + ' > ' + \
        str(len(csvfilelist)) + ' = number of csvs'

    # check either accessionid or seq provided
    assert accessionid is not None or seq is not None, 'accessionid or sequence required'

    # if sequence not provided, get sequence from accessionid
    if len(seq) == 0:
        seq = get_seq(accessionid)
    pepseq = PeptideSequence(seq)

    # init siteinfodict
    site_to_glycans = {x: [] for x in pepseq.n_glycan_sequon_sites}

    # use csv index as replicate ID
    # add initial glycans to siteinfodict
    for replicate, content in enumerate(csvfilelist):
        site_to_glycans = update_siteinfodict(site_to_glycans, content, replicate, accessionid,
                                              pepseq, score_cutoff=score_cutoff)

    # prune off all glycans with replicate count < replicate_cutoff
    # averaging remaining glycans
    site_to_glycans = prune_siteinfodict(
        site_to_glycans, replicate_cutoff=replicate_cutoff)

    # build dict relating site to sequon and glycan list, then return
    sequoninfodict = {seq[x:x + sequon_length]: {'site': x,
                                                 'glycans': y} for x, y in site_to_glycans.items()}
    return sequoninfodict


def write_distribution_csvs(sequoninfodict, outdir):
    """Writes CSV files for the glycan distribution of each site with coverage to the
    provided output directory. No return value."""

    for sequon, sdict in sequoninfodict.items():
        site = sdict['site']
        glycans = sdict['glycans']
        if len(glycans) == 0:
            continue
        gcomp_to_signal = {x.glycan: x.total_signal for x in glycans}

        # turn glycans into ordered dicts and order
        gdict_signal = []
        for gcomp, signal in gcomp_to_signal.items():
            gdict = str_to_dict(gcomp.replace(';', ','))
            gdict = OrderedDict(sorted(gdict.items(),
                                       key=lambda x: tuple(x[0] == y for y in monosaccharide_order), reverse=True))
            gdict_signal.append((gdict, signal))

        # order glycans by complexity; then order scores
        ordered_glycans_scores = sorted(gdict_signal, key=lambda x: tuple(
            x[0][y] for y in glycan_order if y in x[0]))
        ordered_scores = [x[1] for x in ordered_glycans_scores]

        # get normalized scores
        ordered_glycan_strings = [get_glycan_string(
            x[0]) for x in ordered_glycans_scores]
        norm_ordered_scores = [x / sum(ordered_scores) for x in ordered_scores]

        # get csv
        fname = '{}/{}_{}_glycan_distibution.csv'.format(outdir, sequon, site)
        header = ['Glycan composition', 'Probability']
        with open(fname, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)
            for glycan, score in zip(ordered_glycan_strings, norm_ordered_scores):
                row = [glycan, round(score, 3)]
                writer.writerow(row)
    return


def get_distribution_axes(sequoninfodict):
    """Builds glycan distribution plots for sequons with coverage. Returns a dictionary of the form
    {("SEQUON", SITE): AXIS, ...}"""
    sequon_site_to_axis = {}

    for sequon, sdict in sequoninfodict.items():
        site = sdict['site']
        glycans = sdict['glycans']
        if len(glycans) == 0:
            continue
        gcomp_to_signal = {x.glycan: x.total_signal for x in glycans}

        # turn glycans into ordered dicts and order
        gdict_signal = []
        for gcomp, signal in gcomp_to_signal.items():
            gdict = str_to_dict(gcomp.replace(';', ','))
            gdict = OrderedDict(sorted(gdict.items(), key=lambda x: tuple(
                x[0] == y for y in monosaccharide_order), reverse=True))
            gdict_signal.append((gdict, signal))

        # order glycans by complexity; then order scores
        ordered_glycans_scores = sorted(gdict_signal, key=lambda x: tuple(
            x[0][y] for y in glycan_order if y in x[0]))
        ordered_scores = [x[1] for x in ordered_glycans_scores]

        # get normalized scores and colors
        colors = [get_glycan_color(x[0]) for x in ordered_glycans_scores]
        ordered_glycan_strings = [get_glycan_string(
            x[0]) for x in ordered_glycans_scores]
        norm_ordered_scores = [x / sum(ordered_scores) for x in ordered_scores]

        # get plot
        title = '{} ({})'.format(sequon, site)
        axis = plot_bar_chart(
            norm_ordered_scores, ordered_glycan_strings, '', 'Probability', title, colors=colors)
        sequon_site_to_axis[(sequon, site)] = axis

    return sequon_site_to_axis


def get_seq(accession_id, return_id=False):
    """Queries Uniprot for the provided accession ID and returns the sequence. If the accession ID is not found,
    a values of None will be returned."""

    # build and query for url
    url = 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=' + accession_id +\
          '&limit=1&force=no&format=fasta'
    resp = get(url)
    assert resp.status_code == 200, \
        'query for accession id ' + accession_id + \
        ' returned status code ' + str(resp.status_code)
    # get sequence from returned text and return
    # seq will be all lines but first, since first line will be fasta carrot
    # line
    seq = None
    if resp.text != '':
        id_seq = parse_fasta_text(resp.text)
    else:
        raise ValueError("Could not get sequence for %s" % accession_id)

    if return_id is False:
        id_seq = id_seq[1]
    return id_seq


def parse_fasta_text(fasta_text):
    """Returns sequence and ID from text of single fasta entry. ID is the
    first word after > character in header"""
    assert fasta_text.count(
        '>') == 1, 'Only one sequence allowed; ' + str(fasta_text.count('>')) + ' found'
    id = fasta_text.split('\n')[0].split(' ')[0].replace('>', '')
    seq = ''.join(fasta_text.split('\n')[1:])
    return id, seq


def update_siteinfodict(siteinfodict, csvcontent, replicate_id, accession_id, peptidesequence, score_cutoff=30):
    """Updates siteinfodict with csvcontent and returns updated siteinfodict. siteinfodict is
    mapping of glycosylation sites to list of SiteSpecificGlycan object found at that site.
    csvcontent is DictReader content of GlycReSoft glycopeptide-identification csv output.
    replicate_id is the ID of the csv. accession_id is the ID of the protein to which the
    siteinfo dict corresponds. peptidesequence is a PeptideSequence object of the AA sequence
    associated with accession_id. score_cutoff defaults to 30."""

    required_keys = set(['protein_name', 'glycopeptide',
                         'total_signal', 'peptide_start', 'peptide_end', 'ms2_score'])

    for i in csvcontent:
        # check required keys are present
        assert required_keys.issubset(
            i.keys()), 'row keys are ' + str(i.keys())

        # check accid is correct and ms2_score >= score_cutoff
        accid = i['protein_name'].split('|')[1]
        if accid != accession_id or float(i['ms2_score']) < score_cutoff:
            continue

        # get glycopeptide attributes
        glycopeptide = i['glycopeptide']
        composition = '{' + glycopeptide.split('{')[1]
        signal = float(i['total_signal'])
        start = int(i['peptide_start'])
        end = int(i['peptide_end'])
        ms2_score = float(i['ms2_score'])

        # get site
        assert glycopeptide.count('(N-Glycosylation)') == 1
        site = get_matching_nglycan_site(
            glycopeptide, start, end, peptidesequence)
        assert site in siteinfodict.keys(), 'site ' + str(site) + \
            ' not in infodict keys ' + str(siteinfodict.keys())

        # try to find glycan by comp
        # if found, update comp obj
        found = False
        for j in siteinfodict[site]:
            if j.glycan == composition:
                j.add_replicate(replicate_id)
                j.increment_observations()
                j.update_signal(signal)
                found = True
                break

        # if found, go to next line
        # otherwise, add SiteSpecificGlycan
        if found:
            continue
        glycan = SiteSpecificGlycan(composition, site, accid, signal, [
                                    replicate_id], ms2_score)
        siteinfodict[site].append(glycan)

    # return updated dict
    return siteinfodict


def get_glycosite_from_glycopeptide_sequence(glycopeptide_sequence, start_position):
    """Calculates position of glycosylation site on the protein based on glycopeptide
    string and glycopeptide start position. Returns calculated position."""

    splitseq = re.split('\(N-Glycosylation\)', glycopeptide_sequence)
    assert len(splitseq) == 2, 'Multiple n-glycosylations on glycopeptide ' + \
        glycopeptide_sequence

    # take first part of AA sequence, remove other modifications, and rejoin
    # again
    subseq = splitseq[0]
    splitseq = re.split('\([^\(\)]*\)', subseq)
    subseq = ''.join(splitseq)

    # strip off "-" at n-term if present
    if subseq[0] == '-':
        subseq = subseq[1:]

    return start_position + len(subseq) - 1


def get_matching_nglycan_site(gpep_sequence, gpep_start, gpep_end, peptidesequence):
    """Calculates position of glycosylation site by gpep_sequence and gpep_start,
    and check that calculated site is valid based on PeptideSequence sequon sites.
    Returns position of glycosylation site."""

    sites_from_protein = [x for x in peptidesequence.n_glycan_sequon_sites if x >= gpep_start
                          and x <= gpep_end]
    site_from_peptide = get_glycosite_from_glycopeptide_sequence(
        gpep_sequence, gpep_start)

    # check sites are valid before returning
    assert site_from_peptide in sites_from_protein, 'Calculated site = ' + str(site_from_peptide) +\
        ' != sites from pepseq = ' + \
        str(sites_from_protein) + ' for glycopeptide ' + gpep_sequence
    return site_from_peptide


def prune_siteinfodict(siteinfodict, replicate_cutoff=2):
    """Averages glycans in siteinfodict and removes glycans with fewer replicates than
    replicate_cutoff. Returns pruned siteinfodict."""

    for site, glycanlist in siteinfodict.items():
        pruned_glycanlist = []
        for g in glycanlist:
            # if glycan observed in too few replicates, do not use
            if len(g.replicate_ids) < replicate_cutoff:
                continue

            # average and append glycan to pruned list
            calc_average = g.total_signal / g.num_observations
            g.average()
            assert g.total_signal == calc_average, 'Signal for glycan after averaging = ' +\
                str(round(g.total_signal, 3)) + \
                ' != calculated average ' + str(round(calc_average, 3))
            pruned_glycanlist.append(g)

        # update siteinfodict
        siteinfodict[site] = pruned_glycanlist

    return siteinfodict


def get_matched_sequon_infodict(expected_sequoninfodict, experimental_sequoninfodict):
    """Matches glycosylation sites by sequon and relates sequon to site pair
    and list of PairedGlycan Objects ordered by glycan composition. Returns
    matched sequondict of form
    {
        "SEQUON":{
            "sitepair": (EXPECTED_SITE, EXPERIMENTAL_SITE),
            "glycans": [ glycan, ... ]
        },
        ...
    }"""
    matched_dict = {}
    for sequon, expected_sequondict in expected_sequoninfodict.items():
        # skip this sequon if no match or no glycans
        # otherwise, get sequondict for experimental
        if sequon not in experimental_sequoninfodict or \
                len(expected_sequondict['glycans']) == 0 or \
                len(experimental_sequoninfodict[sequon]['glycans']) == 0:
            continue
        experimental_sequondict = experimental_sequoninfodict[sequon]

        # init entry in matched_dict
        matched_dict[sequon] = {'sitepair': (expected_sequondict['site'], experimental_sequondict['site']),
                                'glycans': []}

        # get glycan sets union
        expected_glycans = expected_sequondict['glycans']
        experimental_glycans = experimental_sequondict['glycans']
        glycan_comps = set(
            [x.glycan for x in expected_glycans + experimental_glycans])

        # for each glycan, turn into OrderedDict and build PairedGlycan object
        # before adding to glist
        glist = []
        for gcomp in glycan_comps:
            # get signals for expected and experimental
            expected_signal = get_signal_by_composition(
                expected_glycans, gcomp)
            experimental_signal = get_signal_by_composition(
                experimental_glycans, gcomp)

            # transform gcomp into OrderedDict
            gdict = str_to_dict(gcomp.replace(';', ','))
            gdict = OrderedDict(sorted(gdict.items(),
                                       key=lambda x: tuple(
                                           x[0] == y for y in monosaccharide_order),
                                       reverse=True))

            # make PairedGlycan object and append to matched_dict glycan list
            paired = PairedGlycan(gdict, expected_sequondict['site'], experimental_sequondict['site'],
                                  expected_signal, experimental_signal)
            glist.append(paired)

        # order glist and add to matched_dict
        ordered_glist = sorted(glist, key=lambda x: tuple(
            x.glycan[y] for y in glycan_order if y in x.glycan))
        matched_dict[sequon]['glycans'] = ordered_glist

    # return matched_dict
    return matched_dict


def get_signal_by_composition(glycanlist, composition, not_found=0):
    """Searches list of SiteSpecificGlycan objects for object with specified
    composition (string). Returns total_signal attribute once found. If not
    found, returns not_found value."""
    signal = not_found
    for glycan in glycanlist:
        if glycan.glycan == composition:
            signal = glycan.total_signal
            break
    return signal


def str_to_dict(s):
    """Takes sring of form {a:1, b:2, c:3} and returns it as a dictionary.
    Keys of returned dictionary will be strings. Values cannot be strings."""
    # add double quotes to end of keys
    s = re.split(':', s)
    s = '":'.join(s)

    # add double quotes to front of keys
    s = re.split('{| ', s)
    s = '{' + '"'.join(s)

    return literal_eval(s)


def get_comparison_calculation_dict(matched_dictionary, expected_distribution_name, experimental_distribution_name):
    """Takes a dictionary of sequons mapping to site coordinates and a list of PairGlycans. Plots glycan
    distibutions and calculates cosine similarity and normalized Kullback-Leibler divergence. Returns
    dictionary of format
    {
        "SEQUON": {
            "sitepair": (EXPECT_SITE, EXPER_SITE),
            "cosine": COSINE,
            "nkld": NKLD,
            "comparision": COMPARISON DISTRIBUTION PLOT
        },
        ...
    }"""
    sequon_calculation_axes = {}

    # build axes and calculations from matched dict sequons
    for sequon, sequondict in matched_dictionary.items():
        expected_site, experimental_site = sequondict['sitepair']
        glycans = sequondict['glycans']  # already ordered by composition

        # get and score vectors
        expected_signals = [x.expected_signal for x in glycans]
        experimental_signals = [x.experimental_signal for x in glycans]

        # normalize score vectors
        expected_norm_signals = [x / sum(expected_signals)
                                 for x in expected_signals]
        experimental_norm_signals = [
            x / sum(experimental_signals) for x in experimental_signals]
        assert round(sum(expected_norm_signals), 4) == 1
        assert round(sum(experimental_norm_signals), 4) == 1

        # calculate cosine similarity and normalized Kullback-Leibler
        # divergence
        cosine = cosine_similarity(
            expected_norm_signals, experimental_norm_signals)
        nkld = normalized_kullback_leibler_divergence(
            expected_norm_signals, experimental_norm_signals)

        # get comparison axis
        comparison_title = '{} Glycan Composition Comparision\nCosine similarity: {}\nnKLD: {}'.format(
            sequon, round(cosine, 3), round(nkld, 3))
        glycan_strings = [get_glycan_string(x.glycan) for x in glycans]
        comparison_axis = get_comparison_distribution_axis(expected_norm_signals, experimental_norm_signals,
                                                           expected_distribution_name, experimental_distribution_name,
                                                           comparison_title, glycan_strings)

        # save in dict
        sequon_calculation_axes[sequon] = {
            'sitepair': (expected_site, experimental_site),
            'cosine': cosine,
            'nkld': nkld,
            'comparison': comparison_axis
        }

    # return completed dict
    return sequon_calculation_axes


def get_single_distribution_axis(probabilities, pairedglycans, charttitle, xtitle='Glycans', ytitle='Probability', width=0.35, color=True):
    """Given a list of probabilities and corresponding PairedGlycan objects,
    returns plot of probability distribution"""
    assert round(sum(probabilities),
                 4) == 1, 'Not a probability distribution; sum is ' + str(sum(probabilities))
    assert len(probabilities) == len(pairedglycans), 'Num probabilities = ' + str(len(probabilities)) +\
        ' != ' + str(len(pairedglycans)) + ' = num glycans'

    # remove zeroes, get colors if true, and get glycan string
    filtered_probabilities = []
    colors = []
    glycanstrings = []
    for i, j in zip(probabilities, pairedglycans):
        # do not add if probability is zero
        if i == 0:
            continue
        filtered_probabilities.append(i)
        colors.append(get_glycan_color(j.glycan))
        glycanstrings.append(get_glycan_string(j.glycan))

    # get axis
    if color:
        axis = plot_bar_chart(filtered_probabilities, glycanstrings, xtitle, ytitle, charttitle,
                              width=width, colors=colors)
    else:
        axis = plot_bar_chart(filtered_probabilities, glycanstrings, xtitle, ytitle, charttitle,
                              width=width)
    return axis


def get_glycan_color(glycan, low='crimson', medium='purple', high='navy'):
    """Takes OrderedDict glycan composition and returns color based on number of HexNAc
    (approximation of complexity). Low, medium, and high glycans will be colored red,
    purple, and blue respectively, unless otherwise specified."""
    assert type(glycan) == dict or type(glycan) == OrderedDict, 'glycan type is ' + str(type(glycan)) +\
        'which is not valid'
    assert 'HexNAc' in glycan.keys(), 'HexNAc not in glycan keys: ' + \
        str(glycan.keys())
    count = glycan.get('HexNAc')

    color = None
    if count <= 2:
        color = low
    elif count <= 4:
        color = medium
    else:
        color = high
    return color


def get_glycan_string(glycan, abbreviated=True):
    """Builds glycan string for given glycan. Glycan must be dict or OrderedDict. If abbreviated is true,
    returned string will be of form 'FUC,HEX,HEXNAC,NEU5AC'. Otherwise, returned
    string will be 'Fuc:FUC, Hex:HEX, HexNAc:HEXNAC, Neu5Ac:NEU5AC'."""

    counts = [str(glycan.get(x, 0)) for x in monosaccharide_order]

    s = ''
    if abbreviated:
        s = ','.join(counts)
    else:
        s = ', '.join(
            [x + ':' + y for x, y in zip(monosaccharide_order, counts)])

    return s


def plot_bar_chart(values, xlabels, xtitle, ytitle, charttitle, width=0.35, colors=None):
    # init
    ind = np.arange(len(values))
    fig, axis = plt.subplots()
    if colors is None:
        rects = axis.bar(ind, values, width=width)
    else:
        rects = axis.bar(ind, values, width=width, color=colors)

    # remove top and right borders
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    # set axis labels and title
    axis.set_ylabel(ytitle)
    axis.set_xlabel(xtitle)
    axis.set_title(charttitle)

    # add x-axis ticks and labels
    axis.set_xticks(ind)
    axis.set_xticklabels(xlabels, ha='center', rotation=90)

    return axis


def get_comparison_distribution_axis(distribution1, distribution2, name1, name2, charttitle, xlabels, xtitle='Glycans', ytitle='Probability', width=0.35, color1='indianred', color2='steelblue'):
    """Plots two probability distributions and returns the axis. Given names will be used in a key."""
    assert len(distribution1) == len(distribution2), 'Distribution sizes not equal: ' + str(len(distribution1)) +\
        ' != ' + str(len(distribution2))

    ind = np.arange(len(distribution1))
    fig, axis = plt.subplots()
    rects1 = axis.bar(ind, distribution1, width=width, color=color1)
    rects2 = axis.bar(ind + width, distribution2, width=width, color=color2)

    # remove top and right borders
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    # set axis labels and titles
    axis.set_ylabel(ytitle)
    axis.set_xlabel(xtitle)
    axis.set_title(charttitle)

    # add ticks and labels
    axis.set_xticks(ind + width / 2)
    axis.set_xticklabels(xlabels, rotation=90)

    # add key
    axis.legend((rects1[0], rects2[0]), (name1, name2))

    return axis


def save_plot(filename, axis):
    """Writes axis to filename. Simple function; used for one-line convenience.
    Does not return value."""
    with open(filename, 'w') as f:
        axis.figure.savefig(filename, bbox_inches='tight')
    return

if __name__ == '__main__':

    # arguments
    parser = argparse.ArgumentParser(description='Takes sets of csvs produced with glycresoft for expected and experimental samples and UniProt accession IDs for the expected and experimental proteins to be examined. Creates probabiltiy distributions of the glycans at each glycosylation site and compares the distributions using cosine similarity and normalized Kullback-Leibler divergence.')
    parser.add_argument('-xc', '--expect-csv', dest='expect_csvs', action='append',
                        required=True, help='CSV files for the expected distribution')
    parser.add_argument('-pc', '--exper-csv', dest='exper_csvs', action='append',
                        required=True, help='CSV files for the experimental distribution')
    parser.add_argument('-xa', '--exect-accid', dest='expect_accid',
                        required=True, help='accession ID of expected protein')
    parser.add_argument('-pa', '--exper_accid', dest='exper_accid',
                        required=True, help='accession ID of experimental protein')
    parser.add_argument('-o', '--outdir', dest='outdir',
                        required=True, help='output directory for plots and CSV')
    parser.add_argument('-xn', '--expect-name', dest='expect_name',
                        help='expected distribution key label; defaults to -xa')
    parser.add_argument('-pn', '--exper-name', dest='exper_name',
                        help='experimental distribution key label; defaults to -pa')
    parser.add_argument('-r', '--rcutoff', dest='replicate_cutoff', default=1,
                        type=int, help='replicate observation number required; defaults to 1')
    parser.add_argument('-s', '--scutoff', dest='score_cutoff',
                        default=30, type=int, help='MS2 score required; defaults to 30')
    parser.add_argument('-l', '--length', dest='sequon_length', default=4,
                        type=int, help='length of sequons to match; defaults to 4')
    parser.add_argument("-f", "--fasta", dest="protein_fasta", required=False, )
    args = parser.parse_args()

    # check number of csvs is greater than replicate cutoff
    if args.replicate_cutoff > len(args.expect_csvs):
        print 'Replicate cutoff number cannot be greater than number of expected CSVs'
    elif args.replicate_cutoff > len(args.exper_csvs):
        print 'Replicate cutoff number cannot be greater than number of experimental CSVs'

    # if names not set, set to accession IDs
    if args.expect_name is None:
        args.expect_name = args.expect_accid
    if args.exper_name is None:
        args.exper_name = args.exper_accid

    protein_dict = {}
    if args.protein_fasta is not None:
        for prot in fasta.ProteinFastaFileParser(args.protein_fasta):
            assert prot
            key = None
            try:
                key = prot.name.accession
            except (AttributeError, KeyError) as e:
                key = prot.name[0]
            protein_dict[key] = str(prot)
        expect_seq = protein_dict[args.expect_accid]
        exper_seq = protein_dict[args.exper_accid]
    else:
        expect_seq = exper_seq = ""

    # get file objects
    print 'Reading in files'
    expected_csvfiles = [DictReader(open(x, 'r')) for x in args.expect_csvs]
    experimental_csvfiles = [DictReader(open(x, 'r')) for x in args.exper_csvs]

    # get distibutions for each set
    print 'Getting glycan distributions for {} and {}'.format(args.expect_name, args.exper_name)
    expected_sequon_info = get_sequoninfodict_from_files_accid(
        expected_csvfiles, args.expect_accid, seq=expect_seq, score_cutoff=args.score_cutoff,
        replicate_cutoff=args.replicate_cutoff, sequon_length=args.sequon_length)
    expected_sequon_to_axis = get_distribution_axes(expected_sequon_info)

    experimental_sequon_info = get_sequoninfodict_from_files_accid(
        experimental_csvfiles, args.exper_accid, seq=exper_seq, score_cutoff=args.score_cutoff,
        replicate_cutoff=args.replicate_cutoff, sequon_length=args.sequon_length)
    experimental_sequon_to_axis = get_distribution_axes(
        experimental_sequon_info)

    # match by sequon
    print 'Matching by sequon'
    matched_dict = get_matched_sequon_infodict(
        expected_sequon_info, experimental_sequon_info)

    # get plots, nKLD, and cosine similarity
    print 'Getting axes and similarity measurements'
    calculation_axes = get_comparison_calculation_dict(
        matched_dict, args.expect_name, args.exper_name)

    # write csv
    if not exists(args.outdir):
        makedirs(args.outdir)
    header = ['Sequon', args.expect_name, args.exper_name, 'Cosine', 'nKLD']
    csv_output = '{}/cos_nkld.csv'.format(args.outdir)
    print 'Writing csv to {}'.format(csv_output)
    with open(csv_output, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for sequon, sdict in calculation_axes.items():
            row = [sequon, sdict['sitepair'][0], sdict['sitepair'][1],
                   str(round(sdict['cosine'], 3)), str(round(sdict['nkld'], 3))
                   ]
            writer.writerow(row)

    # write combined distribution plots
    # go in order expected, experimental, comparison
    comparison_subdir = '{}/comparison_distributions'.format(args.outdir)
    if not exists(comparison_subdir):
        makedirs(comparison_subdir)
    for sequon, sdict in calculation_axes.items():
        comparison_sites = '{}_{}'.format(
            sdict['sitepair'][0], sdict['sitepair'][1])
        fname = '{}/{}_{}.svg'.format(comparison_subdir,
                                      sequon, comparison_sites)
        save_plot(fname, sdict['comparison'])

    # save expected plots
    expected_subdir = '{}/{}_distributions'.format(
        args.outdir, args.expect_accid)
    if not exists(expected_subdir):
        makedirs(expected_subdir)
    for (sequon, site), axis in expected_sequon_to_axis.items():
        fname = '{}/{}_{}.svg'.format(expected_subdir, sequon, site)
        save_plot(fname, axis)

    # save experimental plots
    experimental_subdir = '{}/{}_distributions'.format(
        args.outdir, args.exper_accid)
    if not exists(experimental_subdir):
        makedirs(experimental_subdir)
    for (sequon, site), axis in experimental_sequon_to_axis.items():
        fname = '{}/{}_{}.svg'.format(experimental_subdir, sequon, site)
        save_plot(fname, axis)
