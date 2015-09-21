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

import sys
import argparse
import os.path

# patients dict stores instance of class Patient keyed by patient
# id string.  we create a new one when we run into a patient id
# string that we haven't seem before.  patient id string is
# derived from the first part of the name on an alignment,
# e.g. 'patient1_100_1' yields patient id 'patient1'
#
patients = defaultdict()

startDate = datetime.strptime("6/13/1994", "%m/%d/%Y").date()

def render(patients, outgroup, template, fp):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"))
    # Alias str.format to strformat in template
    env.filters['strformat'] = str.format
    template = env.get_template(os.path.abspath(template))
    template.stream(
            patients=patients,
            outgroup=outgroup,
            date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            user=getpass.getuser(),
            command=" ".join(sys.argv),
            workdir=os.getcwd()).dump(fp)


    
def processFasta(datafile, generations=None):
    '''
    Read sequences from a FASTA file (datafile) and create a nested data structure thet organizaes the sequences by patient and sample date.
    if 'generations' is set (list of strings or list of numbers?), only select sequences from the indicated generations.
    '''
    patient = defaultdict(dict)
    
    # define a regex to extract the generation number from the fasta id string
    # we use this to provide tip dates to BEAST.
    patientId = ""
    with open(datafile, "rU") as handle:
    
        # for each fasta sequence in the data file, create a taxon node and a sequence node.
        for record in SeqIO.parse(handle, "fasta") :
            # extract the patient id and generation from the fasta name.
            fields = record.id.split('|')
            patientId = fields[0]
            sampleDate = fields[4] if len(fields) > 3 else "0"
            taxon = record

            collectiondate = patient[sampleDate]
            if not collectiondate:
                collectiondate['taxa'] = []
                collectiondate['date'] = sampleDate

            collectiondate['taxa'].append(taxon)
    return(patientId, patient)



    

def build_parser():
    """
    Build the command-line argument parser.
    """
    def commaSplitter(str):
        """
        Argparse a comm-seperated list
        """
        # leave this here as a reminder of what I shold do to make the argument parsing more robust

        # if sqrt != int(sqrt):
        #      msg = "%r is not a perfect square" % string
        #      raise argparse.ArgumentTypeError(msg)
        # return value
        return str.split(',')

    def existing_file(arg):
        """
        Argparse type for an existing file
        """
        fname,generation = arg.partition(":")[::2]
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        if generation == '':
            generation = None
        return [fname,generation]

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
    
    patients = dict([processFasta(*datafile) for datafile in a.datafiles])

    outgroup = list(SeqIO.parse(a.outgroup, "fasta")) if a.outgroup else []

    render(patients, outgroup, a.template, sys.stdout)



    
if __name__ == "__main__":
   main(sys.argv[1:])
   


