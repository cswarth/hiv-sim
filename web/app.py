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

from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.resources import CDN
from bokeh.templates import RESOURCES
from bokeh.util.string import encode_utf8
from bokeh.models import ColumnDataSource, FactorRange, TapTool, HoverTool, OpenURL


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
        founder.id = 'founder'

    prank = prankseq(os.path.join("../sims/runs/replicate_0", transmit, tsi_donor, tsi_acceptor))
    prank.id = 'prank'
    beast = beastseq(os.path.join("../sims/runs/replicate_0", transmit, tsi_donor, tsi_acceptor))
    beast.id = 'beast'

    (nreplicates, data) = load_data()
    data.set_index(['xmit','dtsi','rtsi'], inplace=True)
    irow = data.loc[transmit, tsi_donor, tsi_acceptor]

    vars['transmit'] = transmit
    vars['tsi_donor'] = tsi_donor
    vars['tsi_acceptor'] = tsi_acceptor

    vars['founder'] = str(founder.seq)
    vars['beast'] = str(beast.seq)
    vars['prank'] = str(prank.seq)
    vars['tree'] = 'sometree'
    vars['directory'] = '/'.join([transmit, tsi_donor, tsi_acceptor])

    for id in list(data):
        vars[id] = irow[id]
    vars['count'] = data.shape[0]
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


colors = {
    'Black': '#000000',
    'Red':   '#FF0000',
    'Green': '#00FF00',
    'Blue':  '#0000FF',
}


def getitem(obj, item, default):
    if item not in obj or len(obj[item])==0:
        return default
    else:
        return obj[item]

def mklegend(colors, values, title="Legend"):
    """Create a legend for a tile plot

    I have not been able to determine if bokeh support legends for tile
    plots, so this is a poor-man's legend
    implemented as a separate plot.  This routine is incestuous with
    map_colors() because it has to know how that routine maps values to
    colors.

    :param colors: the list of colors (in order) used to color-code values in the tile plot
    :param values: the complete list of values that are mapped to colors
    :returns: 
    :rtype:
    """
    
    # make the labels for the tick marks on the legend
    # if all the labels are greater than 1, round to integers
    # otherwise round to three decimal places
    x = np.linspace(*minmax(values)+(len(colors),))
    labels = [np.mean([x[n],x[n-1]]) for n in range(1,len(x))]
    if all(map(lambda x: x > 1, labels)):
        labels = map(lambda x: int(round(x)),labels)
    else:
        labels = map(lambda x: round(x,3),labels)

    legend_ydr = FactorRange(factors=map(str,labels))
    legend_xdr = FactorRange(factors=[1])

    legend = figure(title=title, x_range=legend_xdr, y_range=legend_ydr,
                    plot_width=100, plot_height=250,
                    min_border=5, min_border_left=75, toolbar_location=None, tools=[],
                   x_axis_type=None)
    legend.title_text_font_style = "bold"
    legend.title_text_color = "black"
    legend.title_text_font_size="9pt"

    legend.rect(y=range(1,len(colors)), 
                x=[1] * len(colors), 
                width=1, height=1, 
                fill_color=colors,
                line_color=None)
    return legend

    

    
from bokeh.plotting import *
from bokeh.models import HoverTool, TapTool, OpenURL
from bokeh.models.glyphs import Rect

import bokeh
import numpy as np
import pandas as pd
from collections import OrderedDict
import brewer2mpl

titles = {
    'p_score':'Score of PRANK inferred founder',
    'b_score':'Score of Beast inferred founder',
    'ratio':'Prank / Beast ratio',

    'b_identity':'Beast proportion identical sites',
    'b_gaps':'Proportion gap sites - Beast founder',
    'b_len':'Beast founder length',

    'p_identity':'PRANK proportion identical sites',
    'p_gaps':'Proportion gap sites - PRANK founder',
    'p_len':'PRANK founder length',

    'control':'Color palette - testing and calibration'
}
    
measures = {
    'p_score':'PRANK Score',
    'b_score':'BEAST Score',
    'ratio':'Prank / beast ratio',
    #----
    'p_identity':'Prank % identical',
    'p_gaps':'Prank % gaps',
    'p_len':'Prank founder length',
    #----
    'b_identity':'Beast % identical',
    'b_gaps':'Beast % gaps',
    'b_len':'Beast founder length',
    #----
    'control':'Control'
}

@app.route("/single/")
def index():
    (nreplicates, data) = load_data()

    # save the possible values for the transmission event.
    # we will use these to populate the choices on the web page.
    xmit_events = sorted(data.xmit.unique())
    
    # Grab the inputs arguments from the URL
    args = flask.request.args

    # Get all the form arguments in the url with defaults
    measure = getitem(args, 'measure', 'b_score')
    transmission = getitem(args, 'event', xmit_events[0])
    if len(transmission) == 0:
        transmission = '300'

    data = data[data.xmit == transmission ]

    # use colorbrewer to come up with a nice selection of colors to use
    colors = brewer2mpl.get_map('Reds', 'sequential', 9, reverse=True).hex_colors
    fill_color = map_to_colors(data[measure], colors)

    # A "ColumnDataSource" is like a dict, it maps names to columns of data.
    # These names are not special we can call the columns whatever we like.
    source = ColumnDataSource(
        data=dict(
            donor = data.dtsi,
            recipient = data.rtsi,
            color = fill_color,
            dtime_transmission = data.dtsi,
            dtime_infection = data.rtsi,
            transmission = data.xmit,
            prank = data.p_score, 
            beast=data.b_score,
            ratio=data.ratio,
            n=len(data),
            rprank=data.p_rank,
            rbeast=data.b_rank
        )
    )

    # We need a list of the categorical coordinates
    compare_mixed = lambda x: int(x) if is_number(x) else float("inf")
    x_range = sorted(data.dtsi.unique(), key=compare_mixed)
    #x_range = sorted(list(set(data['dtsi'].tolist())), key=compare_mixed)

    y_range = sorted(data.rtsi.unique(), key= lambda x: int(x))
    #y_range = sorted(list(set(data['rtsi'].tolist())), key= lambda x: int(x))

    tooltips = OrderedDict([
        ("transmission / donor / recipient", "@transmission / @donor / @recipient"),
        ("PRANK","@prank - @rprank"),
        ("BEAST","@beast - @rbeast"),
        ("prank/beast ratio","@ratio")
    ])
    url = "/runs/@transmission/@donor/@recipient/"

    tap = TapTool(action=OpenURL(url=url))
    hover = HoverTool(tooltips=tooltips)

    fig = figure(title=titles[measure].format(min=min(data[measure]), max=max(data[measure])),
                 x_range=x_range, y_range=y_range,
                 x_axis_location="below", plot_width=800, plot_height=600,
                 min_border_right=0,
                 toolbar_location=None, tools=[tap,hover])

    # Prevent the TapTool from highlighting the selected tile
    # https://groups.google.com/a/continuum.io/d/msg/bokeh/ytxc1fQ_6nE/nYIOtzvi1TgJ
    rect = Rect(x='donor', y='recipient', width=.97, height=.97, fill_color='color', line_color=None)
    fig.add_glyph(source, rect, nonselection_glyph=rect)

    # fig.rect(source=source, 
    #         x='donor', y='recipient', width=.97, height=.97, fill_color='color', line_color=None)
    
    fig.xaxis.axis_label="Donor generations since transmission"
    fig.yaxis.axis_label="Recipient generations since infection"

    # draw a second figure as a legend.
    legend = mklegend(colors, data[measure])
    
    # Configure resources to include BokehJS inline in the document.
    # For more details see:
    #   http://bokeh.pydata.org/en/latest/docs/reference/resources_embedding.html#module-bokeh.resources
    plot_resources = RESOURCES.render(
        js_raw=CDN.js_raw,
        css_raw=CDN.css_raw,
        js_files=CDN.js_files,
        css_files=CDN.css_files,
    )

    script, div = components(bokeh.plotting.vplot(bokeh.plotting.hplot(fig, legend)), INLINE)
    template_vars = dict(
        plot_script=script, plot_div=div, plot_resources=plot_resources,
        xmit_events = xmit_events, measures=measures, selected_measure=measures[measure],
        nreplicates=nreplicates)

    html = flask.render_template('embed_bokeh.html', **template_vars)
    return encode_utf8(html)



# map the range of scores into evenly divided bins across the color range
def map_to_colors(values, colors):
    bins = np.linspace(*minmax(values)+(len(colors),))
    values = np.digitize(values, bins)
    colors = [colors[i-1] for i in values]
    return colors

def rank_values(values):
    temp = values.argsort()
    rank = np.empty(len(temp), dtype=np.int32)
    rank[temp] = np.arange(len(temp), dtype=np.int32)
    return rank
    
def load_data():
    data = pd.io.parsers.read_csv('../sims/distances.csv', comment='#')

    # split the 'root' column into separate columns
    # see http://pandas.pydata.org/pandas-docs/stable/generated/pandas.core.strings.StringMethods.extract.html
    s = data.root.str.extract('/replicate_(?P<replicate>\d+)/(?P<xmit>\d+)/(?P<dtsi>.+)/(?P<rtsi>\d+)')
    nreplicates = len(set(s.replicate))

    # assert that the extract worked correctly and there are no NANs in the resulting data frame.
    assert(not s.isnull().any(1).any())

    # join new columns and drop old 'root' column
    data = data.join(s).drop(['replicate', 'root'], axis=1)

    # take the average across replicates
    data = data.groupby(['xmit', 'dtsi', 'rtsi']).mean()
    data.reset_index(inplace=True)

    # calculate the rank of the beast scores
    data['b_rank'] = rank_values(data.b_score)
    data['p_rank'] = rank_values(data.p_score)

    data['ratio'] = data.p_score/data.b_score

    # FYI - if the control column in the distance table is bogus, this is how we can compute our own.
    # data.control = data.rtsi.astype('int') + map(lambda x: 0 if x=='nodonor' else int(x), data.dtsi)
    return (nreplicates, data)


@app.route('/')
def multitile():
    (nreplicates, data) = load_data()
    
    # save the possible values for the transmission event.
    # we will use these to populate the choices on the web page.
    xmit_events = sorted(data.xmit.unique())
    
    # Grab the inputs arguments from the URL
    args = flask.request.args

    # Get all the form arguments in the url with defaults
    measure = getitem(args, 'measure', 'b_score')
    transmission = getitem(args, 'event', xmit_events[0])
    if len(transmission) == 0:
        transmission = '300'

    data = data[data.xmit == transmission ]

    # use colorbrewer to come up with a nice selction of colors to use
    blues = brewer2mpl.get_map('Blues', 'sequential', 9, reverse=True).hex_colors
    reds = brewer2mpl.get_map('Reds', 'sequential', 9, reverse=True).hex_colors

    # Get all the form arguments in the url with defaults
    # display beast scores by default if no other measure is selected
    outer_measure = getitem(args, 'outer_measure', 'b_score')
    inner_measure = getitem(args, 'inner_measure', 'b_score')
    inner_scale = outer_scale = reds
    if outer_measure == inner_measure:
        outer_scale = inner_scale
    else:
        outer_scale = blues
    inner_colors = map_to_colors(data[inner_measure], inner_scale)
    outer_colors = map_to_colors(data[outer_measure], outer_scale) 

    # A "ColumnDataSource" is like a dict, it maps names to columns of data.
    # These names are not special we can call the columns whatever we like.
    source = ColumnDataSource(
        data=dict(
            donor = data.dtsi,
            recipient = data.rtsi,
            inner_color = inner_colors,
            outer_color = outer_colors,
            dtime_transmission = data.dtsi,
            dtime_infection = data.rtsi,
            transmission = data.xmit,
            prank = data.p_score, 
            beast=data.b_score,
            ratio=data.ratio,
            n=len(data),
        )
    )

    # We need a list of the categorical coordinates
    compare_mixed = lambda x: int(x) if is_number(x) else float("inf")
    x_range = sorted(data.dtsi.unique(), key=compare_mixed)
    y_range = sorted(data.rtsi.unique(), key= lambda x: int(x))

    tooltips = OrderedDict([
        ("PRANK","@prank"),
        ("BEAST","@beast"),
        ("prank/beast ratio","@ratio")
    ])
    url = "/runs/@transmission/@donor/@recipient/"

    tap = TapTool(action=OpenURL(url=url))
    hover = HoverTool(tooltips=tooltips)
    fig = figure(title=titles[measure].format(min=min(data[measure]), max=max(data[measure])),
                 x_range=x_range, y_range=y_range,
                 x_axis_location="below", plot_width=800, plot_height=600,
                 min_border_right=0,
                 toolbar_location=None, tools=[tap,hover])

    # Prevent the TapTool from highlighting the selected tile
    # https://groups.google.com/a/continuum.io/d/msg/bokeh/ytxc1fQ_6nE/nYIOtzvi1TgJ
    rect = Rect(x='donor', y='recipient', width=.97, height=.97, fill_color='outer_color', line_color=None)
    fig.add_glyph(source, rect, nonselection_glyph=rect)
    rect = Rect(x='donor', y='recipient', width=.60, height=.60, fill_color='inner_color', line_color=None)
    fig.add_glyph(source, rect, nonselection_glyph=rect)

    # fig.rect(source=source, 
    #         x='donor', y='recipient', width=.97, height=.97, fill_color='color', line_color=None)
    
    fig.xaxis.axis_label="Donor generations since transmission"
    fig.yaxis.axis_label="Recipient generations since infection"


    # draw a second figure as a legend.
    inner_legend = mklegend(inner_scale, data[inner_measure], title=measures[inner_measure])
    outer_legend = mklegend(outer_scale, data[outer_measure], title=measures[outer_measure])

    # Configure resources to include BokehJS inline in the document.
    # For more details see:
    #   http://bokeh.pydata.org/en/latest/docs/reference/resources_embedding.html#module-bokeh.resources
    plot_resources = RESOURCES.render(
        js_raw=CDN.js_raw,
        css_raw=CDN.css_raw,
        js_files=CDN.js_files,
        css_files=CDN.css_files,
    )

    if outer_measure != inner_measure:
        legends = [inner_legend,outer_legend]
    else:
        legends = [inner_legend]
    script, div = components(bokeh.plotting.vplot(bokeh.plotting.hplot(fig, bokeh.plotting.vplot(*legends))), INLINE)

    template_vars = dict(
        plot_script=script, plot_div=div, plot_resources=plot_resources,
        xmit_events = xmit_events,
        measures=measures,
        selected_measure=measures[measure],
        outer_measures=measures,
        inner_measures=measures,
        nreplicates=nreplicates)

    html = flask.render_template('multitile.html', **template_vars)
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
    extra_files = [os.path.join(os.getcwd(), "util.py"),
                   os.path.join(os.getcwd(), "process.py"),
                   os.path.join(os.getcwd(), "filters.py")]
    app.run(host="0.0.0.0", port=5000, debug=True, extra_files=extra_files)
