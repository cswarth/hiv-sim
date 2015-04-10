#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""With code liberally cribbed from Connor McCoy's 'sidbweb' application.

"""
from flask import Flask
from flask import render_template, abort, Response, redirect, url_for, request, g, jsonify
from flask import Flask, make_response
import flask

import re
import os
import sys


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import filters
import process
from util import *

app = Flask(__name__, template_folder='templates/')

filters.register(app)	# register jinja filters in the app


def degapify(s):
    return re.sub('[-_]', '', s)

class FastaIterator():
    """
    leaks file descriptors and I'm too lazy to fix...
    probably the answer is to push filehandle management outside this class so the user can use 'with' contexts to manage their lifetime?
    """
    
    def __init__(self, fname, seqname):
        self.handle = open(fname)
        self.records = SeqIO.parse(self.handle, 'fasta')
        if seqname is not None:
            pattern = re.compile(seqname)
            def namefilter(rec):
                return(pattern.search(rec.id) is not None)
            self.records = ifilter(namefilter, records)

    def __iter__(self):
        return self.records;

    def next(self):
        return self.records.next()

def parse_log(fp):
    """
    parse a beast trait output file.
    this is expected to be a a state output file from beast holding inferred ancestral sequences.
    """
    for line in fp:
        if line[0] == '#':
            continue
        if not line[0].isdigit():
            continue
        line = line.rstrip()
        tokens = line.split()
        tokens = map(lambda t: t[1:-1], tokens)
        yield tokens
        
def str2Seq(id, s):
    """Convert string to SeqRecord insance
    """
    return SeqRecord(Seq(s, IUPAC.IUPACUnambiguousDNA),
              id=id, name="anonymous",
              description="dynamically created sequence form string")
    
def beastseq(path):
    """
    calculate a consensus sequence of ancestral sequences from the posterior.
    Use a 10% burnin period.
    """
    from Bio.Align.Generic import Alignment
    from Bio.Alphabet import IUPAC, Gapped
    from Bio.Align import AlignInfo

    align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
    beastfile = os.path.join(path, 'ancestralSequences.log')
    with open(beastfile) as fh:
        i = 10
        for seq in parse_log(fh):
            # parse_log() returns a list of tokens.
            # first token is mcmc iteration
            # sequences begin at the second token and continue to end of line.
            # we are only looking at the first sequence in the line.
            align.add_sequence("", seq[1])
            i -= 1
            if i == 0:
                break
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.dumb_consensus()
    return SeqRecord(consensus, id='consensus')

                
def prankseq(path):
    prank = os.path.join(path, 'prank.best.anc.fas')
    iseqs = FastaIterator(prank, None)
    pattern = re.compile(r'^#(\d+)#$')

    # identify the root node.  It is the one with the highest 
    rootnode = None
    rootnum = 0
    for seq in iseqs:
        n = pattern.match(seq.id)
        if n is not None:
            n = int(n.group(1))
            if n > rootnum:
                rootnode = seq
                rootnum = n
    return rootnode


@app.route('/runs/<transmit>/<tsi_donor>/<tsi_acceptor>/')
def transmission_detail(transmit, tsi_donor, tsi_acceptor):
    """
    provide a basic landing page that summarizes the information ina single directory of the simulation hierarchy.
    """
    vars = {}

    founderfile = os.path.join("../sims/runs/replicate_0", transmit, 'founder.fa')
    founder = ''
    with open(founderfile, "rU") as handle:
        founder = SeqIO.parse(handle, "fasta").next()

    prank = prankseq(os.path.join("../sims/runs/replicate_0", transmit, tsi_donor, tsi_acceptor))
    beast = beastseq(os.path.join("../sims/runs/replicate_0", transmit, tsi_donor, tsi_acceptor))

    vars['transmit'] = transmit
    vars['tsi_donor'] = tsi_donor
    vars['tsi_acceptor'] = tsi_acceptor

    vars['founder'] = str(founder.seq)
    vars['beast'] = str(beast.seq)
    vars['prank'] = str(prank.seq)
    vars['tree'] = 'sometree'
    vars['directory'] = '/'.join([transmit, tsi_donor, tsi_acceptor])

    prank.id = 'prank'
    beast.id = 'beast'
    founder.id = 'founder'
    vars['sequences'] = process.multi_align(beast, prank, founder)
    vars['highlight'] = highlight(vars['sequences'])
    
    return render_template('transmission_detail.html', **vars)

# http://code.activestate.com/recipes/577879-create-a-nested-dictionary-from-oswalk/
def get_directory_structure(rootdir):
    """
    Creates a nested dictionary that represents the folder structure of rootdir
    """
    dir = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        subdir = dict() # dict.fromkeys(files)
        parent = reduce(dict.get, folders[:-1], dir)
        parent[folders[-1]] = subdir
    return dir

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.resources import CDN
from bokeh.templates import RESOURCES
from bokeh.util.string import encode_utf8


colors = {
    'Black': '#000000',
    'Red':   '#FF0000',
    'Green': '#00FF00',
    'Blue':  '#0000FF',
}


def getitem(obj, item, default):
    if item not in obj:
        return default
    else:
        return obj[item]


@app.route("/test/")
def polynomial():
    """ Very simple embedding of a polynomial chart"""
    # Grab the inputs arguments from the URL
    # This is automated by the button
    args = flask.request.args

    # Get all the form arguments in the url with defaults
    color = colors[getitem(args, 'color', 'Black')]
    _from = int(getitem(args, '_from', 0))
    to = int(getitem(args, 'to', 10))

    # Create a polynomial line graph
    x = list(range(_from, to + 1))
    fig = figure(title="Polynomial")
    fig.line(x, [i ** 2 for i in x], color=color, line_width=2)

    # Configure resources to include BokehJS inline in the document.
    # For more details see:
    #   http://bokeh.pydata.org/en/latest/docs/reference/resources_embedding.html#module-bokeh.resources
    plot_resources = RESOURCES.render(
        js_raw=INLINE.js_raw,
        css_raw=INLINE.css_raw,
        js_files=INLINE.js_files,
        css_files=INLINE.css_files,
    )

    # For more details see:
    #   http://bokeh.pydata.org/en/latest/docs/user_guide/embedding.html#components
    script, div = components(fig, INLINE)
    html = flask.render_template(
        'embed.html',
        plot_script=script, plot_div=div, plot_resources=plot_resources,
        color=color, _from=_from, to=to
    )
    return encode_utf8(html)

from bokeh.plotting import *
from bokeh.models import HoverTool, TapTool, OpenURL
from bokeh.models.glyphs import Rect

import bokeh
import numpy
import pandas as pd
from collections import OrderedDict
import brewer2mpl

measures = {
    'prank':'PRANK',
    'beast':'BEAST',
    'ratio':'Prank / beast ratio',
    'control':'Control'
}
    
@app.route("/")
def index():
    data = pd.io.parsers.read_csv('../sims/distances.tsv', sep='\t')

    # split the 'dir' column into three separate columns
    s = data['dir'].str.split('/').apply(pd.Series)
    s.drop([0], axis=1, inplace=True)
    s.rename(columns={1:'replicate',2:'xmit',3:'dtsi',4:'rtsi'},inplace=True)
    print(s.head())
    s['replicate'] = s.replicate.str.extract('replicate_(\d+)').astype(int)
    print(s.head())
    
    data = data.join(s).drop(data.columns[0], axis=1)
    print(data.head())
    data = data.groupby(['xmit', 'dtsi', 'rtsi'])['prank','beast','control'].mean()
    data.reset_index(inplace=True)
    print(data.head())
    #    prank  beast  control  replicate  xmit   dtsi   rtsi
    # 0  13695  13655     1000          0  1000      0   1600
    # 1  13535  13585     1000          0  1000      0  25000
    # 2  13535  13745     1000          0  1000      0    400
    # 3  13535  13440     1000          0  1000      0   6300
    # 4  13345  12900    25000          0  1000  24000   1600

    # data.groupby([‘col1’, ‘col2’])[‘col3’].mean()

    # save the possible values for the tramission event
    xmit_events = sorted(data.xmit.unique())
    
    # Grab the inputs arguments from the URL
    args = flask.request.args

    # Get all the form arguments in the url with defaults
    measure = getitem(args, 'measure', 'beast')
    transmission = getitem(args, 'event', xmit_events[0])
    print("transmission = {}".format(transmission))
    data = data[data.xmit == transmission ]

    # use colorbrewer to come up with a nice selction of colors to use
    blues = brewer2mpl.get_map('Blues', 'sequential', 9).hex_colors
    reds = brewer2mpl.get_map('Reds', 'sequential', 9).hex_colors

    colors = reds
    data['ratio'] = data.prank/data.beast

    # map the range of scores into evenly divided bins across the color range
    v = data[measure].tolist()
    bins = numpy.linspace(*minmax(v)+(len(colors),))
    digitized = numpy.digitize(v, bins)
    print(set(digitized))
    
    # calculate the rank of the beast scores
    temp = data.beast.argsort()
    rbeast = numpy.empty(len(temp), int)
    rbeast[temp] = numpy.arange(len(temp))

    # calculate the rank of the prank scores
    temp = data.prank.argsort()
    rprank = numpy.empty(len(temp), int)
    rprank[temp] = numpy.arange(len(temp))

    print(len(data.index))
    # A "ColumnDataSource" is like a dict, it maps names to columns of data.
    # These names are not special we can call the columns whatever we like.
    source = ColumnDataSource(
        data=dict(
            donor = data.dtsi,
            recipient = data.rtsi,
            color = [colors[i-1] for i in digitized],
            dtime_transmission = data.dtsi,
            dtime_infection = data.rtsi,
            transmission = data.xmit,
            prank = data.prank, 
            beast=data.beast,
            n=len(data),
            rprank=rprank,
            rbeast=rbeast
        )
    )

    # We need a list of the categorical coordinates
    compare_mixed = lambda x: int(x) if is_number(x) else float("inf")
    x_range = sorted(data.dtsi.unique(), key=compare_mixed)
    #x_range = sorted(list(set(data['dtsi'].tolist())), key=compare_mixed)

    y_range = sorted(data.rtsi.unique(), key= lambda x: int(x))
    #y_range = sorted(list(set(data['rtsi'].tolist())), key= lambda x: int(x))

    tooltips = OrderedDict([
        ("tranmission / donor / recipient", "@transmission / @donor / @recipient"),
        ("PRANK","@prank - @rprank"),
        ("BEAST","@beast - @rbeast")
    ])
    url = "/runs/@transmission/@donor/@recipient/"

    tap = TapTool(action=OpenURL(url=url))
    hover = HoverTool(tooltips=tooltips)

    fig = figure(title="{} score".format(measure),
                 x_range=x_range, y_range=y_range,
                 x_axis_location="below", plot_width=800, plot_height=600,
                 toolbar_location="left", tools=[tap,hover])

    # Prevent the TapTool from highlighting the selected tile
    # https://groups.google.com/a/continuum.io/d/msg/bokeh/ytxc1fQ_6nE/nYIOtzvi1TgJ
    rect = Rect(x='donor', y='recipient', width=.97, height=.97, fill_color='color', line_color=None)
    fig.add_glyph(source, rect, nonselection_glyph=rect)

    # fig.rect(source=source, 
    #         x='donor', y='recipient', width=.97, height=.97, fill_color='color', line_color=None)
    
    fig.xaxis.axis_label="Donor generations since transmission"
    fig.yaxis.axis_label="Recipient generations since infection"

    # # and add the legend just next to the data
    # x, y = int(x_range[0]), int(y_range[0])
    # for i,color in enumerate(reds):
    #     fig.rect([x], [y], width=100, height=100, color=color, 
    #              x_range=[0,1], 
    #              y_range=[0,1])
    #     #fig.text([x], [y], text=area, angle=0, text_font_size="8pt", text_align="center", text_baseline="middle")
    #     y = y + 1000


    
    # Configure resources to include BokehJS inline in the document.
    # For more details see:
    #   http://bokeh.pydata.org/en/latest/docs/reference/resources_embedding.html#module-bokeh.resources
    plot_resources = RESOURCES.render(
        js_raw=CDN.js_raw,
        css_raw=CDN.css_raw,
        js_files=CDN.js_files,
        css_files=CDN.css_files,
    )

    script, div = components(fig, INLINE)
    template_vars = dict(
        plot_script=script, plot_div=div, plot_resources=plot_resources,
        xmit_events = xmit_events, measures=measures, selected_measure=measures[measure])

    html = flask.render_template('embed_bokeh.html', **template_vars)
    return encode_utf8(html)
    

@app.route('/figtree/<transmit>/<tsi_donor>/<tsi_acceptor>/mcc.svg')
def mcc_tree_svg(transmit, tsi_donor, tsi_acceptor):
    tree = process.tree_svg(os.path.join("../sims/runs/replicate_0", transmit, tsi_donor, tsi_acceptor, 'mcc.tree'), compress=False)
    resp = make_response(tree)
    resp.headers['Content-Type'] = 'image/svg+xml'
    # resp.headers['Content-Encoding'] = 'gzip'

    return resp
 
if __name__ == "__main__":
    # add files to be watched for changes
    # http://stackoverflow.com/a/9511655/1135316
    extra_files = ["/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/util.py",
                   "/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/process.py",
                   "/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/filters.py"]
    app.run(host="0.0.0.0", debug=True, extra_files=extra_files)
