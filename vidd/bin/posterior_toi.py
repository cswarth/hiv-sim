#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
extract an estimate of time of infection from a mcmc posterior distribution.

expects a log file as output by Beast.
Typical usage:

$ posterior_Toi.py outout/beastout.log samples/CAPRISA002_PDC_GP120_1M_aln.fa

'''
from __future__ import print_function

import re
from Bio import SeqIO
from collections import defaultdict
from datetime import datetime, date, timedelta
import calendar
import sys
import argparse
import os.path
import pandas as pd

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

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--burnin', '-b', help='burning percentage, expressed as whole integer [default %default]', default=10)
    parser.add_argument('logfile', nargs=1, help='Beast log file', type=existing_file)
    parser.add_argument('fasta', nargs=1, help='FASTA input', type=existing_file)

    return parser


def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()
    
    tbl = pd.read_table(a.logfile[0], skiprows=2)
    rootHeight = tbl['treeModel.rootHeight']
    n = len(rootHeight)
    burnin = int(n*(a.burnin / 100.0))
    toi = tbl['treeModel.rootHeight'][burnin:].mean()


    patients = dict([processFasta(datafile) for datafile in a.fasta])
    samples = [sample for p in patients.values() for sample in p.values()]
    dates = [s['date'] for s in samples]
    latest_timepoint = max(dates)
    doi = latest_timepoint - timedelta(days=365.0*toi)
    print('Estimated date of infection: {}'.format(datetime.strftime(doi, '%Y/%m/%d')))
    #print(tbl)


    
if __name__ == "__main__":
   main(sys.argv[1:])
   
