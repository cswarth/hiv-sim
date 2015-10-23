#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
split the sequences in a fasta file into multiple files by date and region.

We have sample files with sequence names like, 

>CAPRISA002_PDC|6M|64C2C3|gp120|2006/09/28

These samples files sometimes have sequences from multiple timepoints
mixed in them, like '1m6m' meaning the file has sequences from both the
1 month and 6 month timepoints.  Furthermore, sample files may also
contain sequences from multiple regions, like 'NFLG', 'LH', and 'RH'.
The goal is to split the files so sequences from each
timepoint and region end up in separate files.

Typical usage:

$ find ../../../sample_data -name '*fasta' | grep -v _aln | xargs ./split_samples.py --outdir samples
'''
from __future__ import print_function

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

region_dd = lambda: defaultdict(list)
sample_dd = lambda: defaultdict(region_dd)
patient_dd = lambda: defaultdict(sample_dd)


def processFasta(datafile, patients):
    '''
    Read sequences from a FASTA file (datafile) and create a nested data structure thet organizaes the sequences by patient and sample date.
    '''
    # extract gene  from the filename,
    # as the gene name in the sequence id is unreliable.
    # rv217_pdc_gp120_1m6m_aln.fasta
    # caprisa002_pda_gp120_1m6m_aln.fasta
    fields = os.path.basename(datafile).split('_')
    fgene = fields[2].upper()	# file gene name
    
    patientId = None
    with open(datafile, "rU") as handle:
    
        # for each fasta sequence in the data file, create a taxon node and a sequence node.
        for record in SeqIO.parse(handle, "fasta") :
            # extract the patient id and generation from the fasta name.
            # >RV217_PDA|Other|01LH|LH|2012/07/25
            # >RV217_PDB|1M|01WG|NFLG|2011/11/10
            # >CAPRISA002_PDB|1M|01C2C3|gp120|2006/03/16


            fields = record.id.split('|')
            sgene = fields[3]	# sequence gene name
            # Prefer the gene name from the sequence id unless it is NFLG,
            # otherwise use the gene from the file name.
            gene = fgene if fgene != 'NFLG' else sgene
            taxon = record
            patientId = fields[0]
            timePoint = fields[1] if len(fields) > 0 else "0"
            sampleDate = fields[4] if len(fields) > 3 else "0"

            patient = patients[patientId]
            
            sample = patient[sampleDate]
            if not sample:
                sample['regions'] = region_dd()
                sample['date'] = datetime.strptime(sampleDate, '%Y/%m/%d')
                sample['timepoint'] = timePoint
            region = sample['regions'][gene]
            region.append(taxon)
    
    if patientId is None:
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

    parser.add_argument('--fasta', help='produce a FASTA file (default: produce XML file)',
            action='store_true', default=False, dest='createFasta')
    parser.add_argument('--dryrun', '-n', help='print the names of fils that would be created.',
             action='store_true', default=False, dest='dryrun')
    parser.add_argument('--outdir', help='output directory',
            default='.', dest='outdir')
    parser.add_argument('datafiles', nargs='+', help='FASTA input', type=existing_file)

    return parser


def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()
    patients = patient_dd()

    for datafile in a.datafiles:
        processFasta(datafile, patients)

    samples = [sample for p in patients.values() for sample in p.values()]
    dates = [s['date'] for s in samples]


    for pid, patient in patients.items():
        for sid,sample in patient.items():
            for gene,region in sample['regions'].items():
                if not a.dryrun:
                    print("{}  {} {}  {}".format(pid, gene, sid, sample['timepoint'])) 
                filename = '{}_{}_{}.fa'.format(pid, gene, sample['timepoint'])
                filename = os.path.join(a.outdir, filename)
                if a.dryrun:
                    print(filename)
                else:
                    with open(filename, "w") as fh:
                        SeqIO.write(region, fh, "fasta")

    
if __name__ == "__main__":
   main(sys.argv[1:])
   
