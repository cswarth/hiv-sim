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

import re
from Bio import SeqIO
from collections import defaultdict
import getpass
from datetime import datetime, date, timedelta

import sys
import argparse
import os.path

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
    parser.add_argument('-g', '--generation', help='generation to sample from',
            action='store', default=None, dest='generation')
    parser.add_argument('template', type=argparse.FileType('r'), help='BEAST config template file')
    parser.add_argument('fastaFile', type=argparse.FileType('r'), help='file of sequences (in FASTA format)')
    parser.add_argument('count', nargs='?', type=int, default=-1, help='number of sample sequences to extract')

    return parser



def main():
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()

    sourcefile = a.fastaFile.name
    
    # find the sequences from the chosen <generation> and
    # randomly choose <count> of them to pass on to the next generation
    data = []
    cumulativeGeneration = 0
    with a.fastaFile as fp:
        for record in SeqIO.parse(fp, "fasta") :
            # import pdb; pdb.set_trace()

            fields = record.id.split('|')
            subfields = fields[0].split('_')
            patient = subfields[0]
            generation = subfields[1] if len(subfields) >= 2 else None
            
            if not a.generation or generation == a.generation:
                if a.count != 0:
                    # Randomly choose a sequence.  For now
                    # assume they are already randomy sorted and choose from the
                    # order in which they are presented.
                    #
                    # need to assert the all the lengths are the same....
                    # if not, we will have to align all the sequences to one another...
                    ll = len(record.seq)
                    data.append(record)
                    a.count -= 1
                    cumulativeGeneration = int(fields[3]) if len(fields) >=3 else 0

            if a.count == 0:
                # break out when we have collected all the sequences we need
                break



    # compose a label to indicate where the starting population from the sample
    a.generation = 0 if not a.generation else a.generation
    cumulativeGeneration += int(a.generation)
    
    label= "{}_%g_%s|{}|{}|{}".format(a.prefix, os.path.basename(sourcefile), a.generation, cumulativeGeneration)

    with sys.stdout as fp:
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"))
        # Alias str.format to strformat in template
        env.filters['strformat'] = str.format
        template = env.get_template(os.path.abspath(a.template.name))
        template.stream(data=data,
                label=label,
                date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                user=getpass.getuser(),
                command=" ".join(sys.argv),
                workdir=os.getcwd()).dump(fp)



    
if __name__ == "__main__":
   main()
   


