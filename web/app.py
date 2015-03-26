#!/usr/bin/env python
"""With code liberally cribbed from Connor McCoy's 'sidbweb' application.

"""
from flask import Flask
from flask import render_template, abort, Response, redirect, url_for, request, g, jsonify
from flask import Flask, make_response

import re
import os
import sys


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import filters
import process
import util

app = Flask(__name__, template_folder='templates/')

filters.register(app)	# register jinja filters in the app

@app.route('/')
def index():
    return app.root_path

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


@app.route('/runs/<transmit>/<tsi_donor>/<tsi_acceptor>/<clockmodel>/')
def transmission_detail(transmit, tsi_donor, tsi_acceptor, clockmodel):
    """
    provide a basic landing page that summarizes the information ina single directory of the simulation hierarchy.
    """
    vars = {}

    founderfile = os.path.join("../sims/runs", transmit, 'founder.fa')
    founder = ''
    with open(founderfile, "rU") as handle:
        founder = SeqIO.parse(handle, "fasta").next()

    prank = prankseq(os.path.join("../sims/runs", transmit, tsi_donor, tsi_acceptor, clockmodel))
    beast = beastseq(os.path.join("../sims/runs", transmit, tsi_donor, tsi_acceptor, clockmodel))

    vars['transmit'] = transmit
    vars['tsi_donor'] = tsi_donor
    vars['tsi_acceptor'] = tsi_acceptor
    vars['clockmodel'] = clockmodel

    vars['founder'] = str(founder.seq)
    vars['beast'] = str(beast.seq)
    vars['prank'] = str(prank.seq)
    vars['tree'] = 'sometree'
    vars['directory'] = '/'.join([transmit, tsi_donor, tsi_acceptor, clockmodel])

    prank.id = 'prank'
    beast.id = 'beast'
    founder.id = 'founder'
    vars['sequences'] = process.multi_align(beast, prank, founder)
    vars['highlight'] = util.highlight(vars['sequences'])
    
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



@app.route("/simple.png")
def simple():
    import datetime
    import StringIO
    import random
 
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.dates import DateFormatter
 
    fig=Figure()
    ax=fig.add_subplot(111)
    x=[]
    y=[]
    now=datetime.datetime.now()
    delta=datetime.timedelta(days=1)
    for i in range(10):
        x.append(now)
        now+=delta
        y.append(random.randint(0, 1000))
    ax.plot_date(x, y, '-')
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()
    canvas=FigureCanvas(fig)
    png_output = StringIO.StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response

@app.route('/runs/<transmit>/<tsi_donor>/<tsi_acceptor>/<clockmodel>/mcc.svg')
def mcc_tree_svg(transmit, tsi_donor, tsi_acceptor, clockmodel):
    tree = process.tree_svg(os.path.join("../sims/runs", transmit, tsi_donor, tsi_acceptor, clockmodel, 'mcc.tree'), compress=True)
    resp = make_response(tree)
    resp.headers['Content-Type'] = 'image/svg+xml'
    resp.headers['Content-Encoding'] = 'gzip'

    return resp
 
if __name__ == "__main__":
    # add files to be watched for changes
    # http://stackoverflow.com/a/9511655/1135316
    extra_files = ["/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/util.py", "/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/process.py", "/shared/silo_researcher/Matsen_F/MatsenGrp/working/cwarth/hiv-sim/web/filters.py"]
    app.run(host="0.0.0.0", debug=True, extra_files=extra_files)
