#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
	Create a SANTA config file by combining sequences sampled from a FASTA file with a templated XML file.
    SANTA simulates the evolution of a population of gene sequences forwards through time. It models the
    underlying biological processes as discrete components; replication (including recombination), mutation, fitness and selection.
    See https://code.google.com/p/santa-sim/

    It is expected that most often, a single sequence will be sampled to create a founder sequence for the simulation.
    However if you want to use more than one sequence to found your simulation, you can do so with the 'count' positional argument.
'''

from __future__ import print_function

import jinja2

from Bio import SeqIO
from collections import defaultdict
import getpass
from datetime import datetime, date, timedelta

import sys
import argparse
import os.path
import re

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

    parser.add_argument('-p', '--prefix', help='dont really know what this does...',
            action='store', default='patient', dest='prefix')
    parser.add_argument('-d', '--date', help='dont really know what this does...',
            action='store', default='', dest='sampledate')
    parser.add_argument('template', type=argparse.FileType('r'), help='BEAST config template file')
    parser.add_argument('fasta', type=argparse.FileType('r'), help='file of sequences (in FASTA format)')

    return parser



def main():
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()

    sourcefile = a.fasta.name
    
    # find the sequences from the chosen <generation> and
    cumulativeGeneration = 0

    records = list(SeqIO.parse(a.fasta, "fasta"))

    # extract generation and cumulative generation from each sequence.
    generation = set([int(r.id.split('|')[0].split('_')[1]) for r in records])
    cumulative = set([int(r.id.split('|')[3]) for r in records])
    assert(len(generation) == 1)
    assert(len(cumulative) == 1)
    cumulative = cumulative.pop()
    generation = generation.pop()
    
    # compose a label to indicate where the starting population from the sample
    cumulativeGeneration = cumulative + generation
    label= "{}_%g_%s|{}|{}|{}".format(a.prefix, os.path.basename(sourcefile), generation, cumulativeGeneration)

    with sys.stdout as fp:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"))
        # Alias str.format to strformat in template
        env.filters['strformat'] = str.format
        template = env.get_template(os.path.abspath(a.template.name))
        template.stream(data=records,
                label=label,
                date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                user=getpass.getuser(),
                command=" ".join(sys.argv),
                workdir=os.getcwd()).dump(fp)

    
if __name__ == "__main__":
   main()
   


