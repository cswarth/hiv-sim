#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    make a new BEAST  config file by inserting FASTA sequences from the RV217 trail.
    The sequence names int the RV217 data look like this:
	>RV217_PDB|1M|01WG|NFLG|2011/11/10

    The 5th component is the sample data and should eb converted to an appropriate tipdate for BEAST.

    usage:

         mkbeast.py -p foo template.xml sequences.fasta  >beast_in.xml

    This will take the fasta sequences and insert them into template.xml
    to produce an XML file that is suitable to pass to BEAST.  

    Once the BEAST config file is generated, you would run,

         beast beast_in.xml

    This will produce various output files, among which is foo.trees.
    That file gets fed to the 'annotatetrees' program and the output of that gets
    visualized with 'figtree'.
'''

from __future__ import print_function

import jinja2

import re
from Bio import SeqIO
from collections import defaultdict
import getpass
from datetime import datetime, date, timedelta
import calendar
import sys
import argparse
import os.path
import operator

# patients dict stores instance of class Patient keyed by patient
# id string.  we create a new one when we run into a patient id
# string that we haven't seem before.  patient id string is
# derived from the first part of the name on an alignment,
# e.g. 'patient1_100_1' yields patient id 'patient1'
#
patients = defaultdict()

def dayofyear(value):
    t = value.timetuple()
    return "{}".format(t.tm_year+(t.tm_yday-1)/(366.0 if calendar.isleap(t.tm_year) else 365.0))

def dateformat(value, format='%d/%m/%Y'):
    return value.strftime(format)

def deltayears(value):
    return value.total_seconds()/(60*60*24*365)

def render(params, template, fp):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"))
    env.filters['dayofyear'] = dayofyear
    env.filters['deltayears'] = deltayears
    env.filters['dateformat'] = dateformat
    
    # Alias str.format to strformat in template
    env.filters['strformat'] = str.format
    template = env.get_template(os.path.abspath(template))

    # define additional values available in the template
    vars = dict(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                  user=getpass.getuser(),
                  command=" ".join(sys.argv),
                  workdir=os.getcwd(),
                  now=datetime.now(),
                  timedelta=timedelta)
    # merge in the values passed in to this routine.
    vars.update(params)
    
    template.stream(**vars).dump(fp)

# Convert a string like '1M', '3M', or '6M' to a timedelta object corresponding to number of months indicated in string.
#
def str2timedelta(s):
    m = re.match(r'(\d+)M$', s)
    if m:
        delay = 30 * int(m.group(1))
    else:
        delay = 30 * 12 
    return(timedelta(days=delay))

    
def processFasta(datafile):
    '''
    Read sequences from a FASTA file (datafile) and create a nested data structure thet organizaes the sequences by patient and sample date.
    '''
    patient = defaultdict(dict)
    
    # define a regex to extract the generation number from the fasta id string
    # we use this to provide tip dates to BEAST.
    patientId = None
    with open(datafile, "rU") as handle:
    
        # for each fasta sequence in the data file, create a taxon node and a sequence node.
        for record in SeqIO.parse(handle, "fasta") :
            # extract the patient id and generation from the fasta name.
            fields = record.id.split('|')
            patientId = fields[0]
            timePoint = fields[1] if len(fields) > 0 else "0"
            sampleDate = fields[4] if len(fields) > 3 else "0"
            taxon = record

            collectiondate = patient[sampleDate]
            if not collectiondate:
                collectiondate['taxa'] = []
                collectiondate['date'] = datetime.strptime(sampleDate, '%Y/%m/%d')
                collectiondate['delta'] = str2timedelta(timePoint)

            collectiondate['taxa'].append(taxon)
    
    if patientId is not None:
        return(patientId, patient)
    else:
        raise Exception('Empty fasta file - {}'.format(datafile))



    

def build_parser():
    """
    Build the command-line argument parser.
    """
    def commaSplitter(str):
        """
        Argparse a comm-seperated list
        """
        # leave this here as a reminder of what I should do to make the argument parsing more robust

        # if sqrt != int(sqrt):
        #      msg = "%r is not a perfect square" % string
        #      raise argparse.ArgumentTypeError(msg)
        # return value
        return str.split(',')

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--template', '-t', help='templated BEAST XML config',
            required=True, dest='template')
    parser.add_argument('--fasta', help='produce a FASTA file (default: produce XML file)',
            action='store_true', default=False, dest='createFasta')
    parser.add_argument('--prefix', help='Specify a prefix for all output log filename',
            default="", dest='prefix')
    parser.add_argument('--outgroup', help='FASTA files for outgroup sequences',
            default=None, dest='outgroup', type=argparse.FileType('r'))
    parser.add_argument('datafiles', nargs='+', help='FASTA input', type=existing_file)

    return parser


def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()
    
    patients = dict([processFasta(datafile) for datafile in a.datafiles])

    samples = [sample for p in patients.values() for sample in p.values()]
    dates = [s['date'] for s in samples]

    earliest_timepoint = min([ s['date'] - s['delta'] for s in samples])
    latest_timepoint = max(dates)
    
    outgroup = list(SeqIO.parse(a.outgroup, "fasta")) if a.outgroup else []

    params = dict(patients=patients,
                  outgroup=outgroup,
                  earliest_timepoint=earliest_timepoint,
                  latest_timepoint=latest_timepoint)

    render(params, a.template, sys.stdout)



    
if __name__ == "__main__":
   main(sys.argv[1:])
   


