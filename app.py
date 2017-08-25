#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt

#from flask import Flask, g, jsonify, render_template, Response, request
from flask import Flask, render_template, Response, request
from csv import DictReader
from io import BytesIO
import glycan_distributions
import urllib

def render_plot(figure, format=None, **kwargs):
    if isinstance(figure, matplotlib.axes.Axes):
        figure = figure.get_figure()

    kwargs["format"] = format
    kwargs['bbox_inches'] = 'tight'
    kwargs['patchless'] = True
    if "height" in kwargs:
        figure.set_figheight(kwargs["height"])
    if "width" in kwargs:
        figure.set_figwidth(kwargs['width'])
    if kwargs.get("bbox_inches") != 'tight' or kwargs.get("patchless"):
        figure.patch.set_alpha(0)
        figure.axes[0].patch.set_alpha(0)
    data_buffer = BytesIO()
    figure.savefig(data_buffer, **kwargs)
    plt.close(figure)
    return data_buffer

def svg_plot(figure, **kwargs):
    data_buffer = render_plot(figure, format='svg', **kwargs)
    return etree.tostring(etree.fromstring(data_buffer.getvalue()))

def png_plot(figure, **kwargs):
    data_buffer = render_plot(figure, format='png', **kwargs)
    return "<img src='data:image/png;base64,%s'>" % urllib.quote(data_buffer.getvalue().encode("base64"))

def process_response(request):
    """Processes response object for submission page and checks conditions are met. 
    Returns dict with all variables needed to proceed and files loaded as DictReader objects"""
    rd = {}

    # check csvs
    rd['expect_list'] = [ DictReader(x) for x in request.files.getlist('expect-csvs') ]
    rd['exper_list'] = [ DictReader(x) for x in request.files.getlist('exper-csvs') ]
    assert len(rd['expect_list']) > 0 and len(rd['exper_list']) > 0, \
    'Expected and experimental distributions require at least one CSV each.'

    # check dataset names
    rd['expect_name'] = request.values.get('expect-name')
    rd['exper_name'] = request.values.get('exper-name')
    assert rd['expect_name'] is not None and rd['exper_name'] is not None, \
    'An expected and experimental dataset name must be provided'

    # get accIDs
    rd['expect_accid'] = request.values.get('expect-accid')
    rd['exper_accid'] = request.values.get('exper-accid')

    # get fastas
    rd['expect_fasta'] = request.files.get('expect-fasta').read()
    rd['exper_fasta'] = request.files.get('exper-fasta').read()

    # get seq if fasta and check that both datasets have either accID or valid fasta
    if len(rd['expect_fasta']) > 0:
        # get sequence
        accid, seq = glycan_distributions.parse_fasta_text(rd['expect_fasta'])
        rd['expect_fasta'] = seq
        assert len(rd['expect_fasta']) > 0, \
        'Valid fasta file required for expected dataset'
        # set accid if not given in form
        if rd['expect_accid'] == '':
            rd['expect_accid'] = accid
    # check accid is set
    assert rd['expect_accid'] is not None and len(rd['expect_accid']) > 0, \
    'Accession ID not given or fasta improperly formatted for expected dataset'

    if len(rd['exper_fasta']) > 0:
        # get sequence
        accid, seq = glycan_distributions.parse_fasta_text(rd['exper_fasta'])
        rd['exper_fasta'] = seq
        assert len(rd['exper_fasta']) > 0, \
        'Valid fasta file required for experimental dataset'
        # set accid if not given in form
        if rd['exper_accid'] == '':
            rd['exper_accid'] = accid
    # check accid is set
    assert rd['expect_accid'] is not None and len(rd['expect_accid']) > 0, \
    'Accession ID not given or fasta improperly formatted for expected dataset'

    if len(rd['exper_fasta']) > 0:
        # get sequence
        accid, seq = glycan_distributions.parse_fasta_text(rd['exper_fasta'])
        rd['exper_fasta'] = seq
        assert len(rd['exper_fasta']) > 0, \
        'Valid fasta file required for experimental dataset'
        # set accid if not given in form
        if rd['exper_accid'] == '':
            rd['exper_accid'] = accid
    # check accid is set
    assert rd['exper_accid'] is not None and len(rd['exper_accid']) > 0, \
    'Accession ID not given or fasta improperly formatted for experimental dataset'

    # constants
    rd['rcutoff'] = request.values.get('rcutoff', 1, type=int)
    rd['scutoff'] = request.values.get('scutoff', 30, type=int)
    rd['sequonlen'] = request.values.get('sequonlen', 4, type=int)

    return rd

def get_table(calculation_dict, expected_name, experimental_name):
    """Builts table out of dictionary containing sequons, sitevalues, CS, and nKLD. 
    Returns list of lists for rows of table, with first list being the header"""
    header = [ ['Sequon', expected_name+' site', experimental_name+' site', 'CS', 'nKLD'] ]

    # build rows
    l = []
    for sequon, sdict in calculation_dict.items():
        expect_site = sdict['sitepair'][0]
        exper_site = sdict['sitepair'][1]
        cos = round(sdict['cosine'], 3)
        nkld = round(sdict['nkld'], 3)
        row = [sequon, expect_site, exper_site, cos, nkld]
        l.append(row)

    # sort rows based on expect_site
    l = sorted(l, key=lambda x: x[1])
    l = header+l
    return l

# Flask stuff below
app = Flask(__name__)

@app.route("/")
def index():
    return render_template("home.html")

@app.route("/csv")
def csv():
    return render_template("csv_form.html")

@app.route("/csv-form", methods=["POST"])
def csv_form():
    # process request object and get back values
    processed_response = process_response(request)
    replicate_cutoff = processed_response['rcutoff']
    score_cutoff = processed_response['scutoff']
    sequon_length = processed_response['sequonlen']

    # expected
    expected_files = processed_response['expect_list']
    expected_name = processed_response['expect_name']
    expected_accid = processed_response['expect_accid']
    expected_seq = processed_response['expect_fasta']
    expected_sequon_info = glycan_distributions.get_sequoninfodict_from_files_accid(
        expected_files, accessionid=expected_accid, seq=expected_seq, score_cutoff=score_cutoff, 
        replicate_cutoff=replicate_cutoff, sequon_length=sequon_length)
    expected_sequon_to_axis = glycan_distributions.get_distribution_axes(expected_sequon_info)

    # experimental
    experimental_files = processed_response['exper_list']
    experimental_name = processed_response['exper_name']
    experimental_accid = processed_response['exper_accid']
    experimental_seq = processed_response['exper_fasta']
    experimental_sequon_info = glycan_distributions.get_sequoninfodict_from_files_accid(
        experimental_files, accessionid=experimental_accid, seq=experimental_seq, score_cutoff=score_cutoff, 
        replicate_cutoff=replicate_cutoff, sequon_length=sequon_length)
    experimental_sequon_to_axis = glycan_distributions.get_distribution_axes(experimental_sequon_info)


    # match and get calculation axes
    matched_dict = glycan_distributions.get_matched_sequon_infodict(
        expected_sequon_info, experimental_sequon_info)
    calculation_axes = glycan_distributions.get_comparison_calculation_dict(
        matched_dict, expected_name, experimental_name)
    
    # plotting
    plot_comparisons = [ png_plot(x['comparison']) for x in calculation_axes.values() ]
    plot_expected = [ png_plot(x) for x in expected_sequon_to_axis.values() ]
    plot_experimental = [ png_plot(x) for x in experimental_sequon_to_axis.values() ]

    # build table and CSV string
    value_table = get_table(calculation_axes, expected_name, experimental_name)
    csvstr = ''
    for row in value_table:
        csvstr += ','.join([str(x) for x in row])+'\n'

    # create new template and return
    return render_template("plot_form.html", value_table_header=value_table[0],
        value_table_rows=value_table[1:], plot_comparisons=plot_comparisons, 
        plot_expected=plot_expected, plot_experimental=plot_experimental,
        csvstr=csvstr)

@app.route("/csv-file", methods=["POST"])
def csv_file():
    csv_str = request.values.get('csvstr')
    return Response(csv_str, mimetype='text/csv', 
        headers={'Content-disposition':'attachment; filename=similarities.csv'})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=9898, use_reloader=True, threaded=True)
