#!/usr/bin/env python

# read in a set of sequences as output from beast an caclulate consensus sequences.
#

from __future__ import print_function

from lxml import etree
import re
import math
import StringIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

from collections import defaultdict
from datetime import datetime, date, timedelta

import sys
import argparse
import os.path
import numpy as np

IUPAC_rev = {('A',): 'A',
             ('A', 'C'): 'M',
             ('A', 'C', 'G'): 'V',
             ('A', 'C', 'T'): 'H',
             ('A', 'G'): 'R',
             ('A', 'G', 'T'): 'D',
             ('A', 'C', 'G', 'T'): 'N',
             ('A', 'T'): 'W',
             ('C',): 'C',
             ('C', 'G'): 'S',
             ('C', 'G', 'T'): 'B',
             ('C', 'T'): 'Y',
             ('G',): 'G',
             ('G', 'T'): 'K',
             ('T',): 'T'}


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--burnin', type=int, help='skip this many entries as a burn-in period (default: no burn-in period)',
            action='store', default=0, dest='burnin')
    parser.add_argument('--minpct', type=float, help='minimum percentage at each site necessary to call a base at that site (default: 0.70)',
            action='store', default=0.70, dest='minpct')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    return parser

def is_int(s):
    '''
    return True if string 's' is an integer.
    Used in this script to test if the first field on a line is an int.
    '''
    
    try:
        int(s)
        return True
    except ValueError:
        return False

# inspired by Module Bio.Align.AlignInfo.information_content() 
def tabulate( seqList, start = 0, end = None ):
    """ make a list of dictionaries, each holding nucleotide frequences for a single site in the multiple alignment.
    """
    
    # if no end was specified, then we default to the end of the sequence 
    if end is None: 
        end = len(seqList.alignment._records[0].seq) 
   
    if start < 0 or end > len(seqList.alignment._records[0].seq): 
        raise ValueError("Start (%s) and end (%s) are not in the range %s to %s" 
                        % (start, end, 0, len(seqList.alignment._records[0].seq))) 

    all_letters = seqList._get_all_letters()
    chars_to_ignore = []
    dictList = []
    for residue_num in range(start, end): 
        dictList.append(seqList._get_letter_freqs(residue_num, 
                                            seqList.alignment._records, 
                                            all_letters, chars_to_ignore))

    return dictList

def consensus( tabdict, plu=.33, gap='-', errorchar='X', use_ambi=True ):
    """Given a dictionary representing character frequencies
    at a single position, returns the most common char at
    that position subject to the rules below.

    plu            plurality for calling a consensus character

    use_ambi - uses IUPAC ambiguity codes if possible in place of errorchar
    """

    # remove any base that has a zero count.
    # convert fractional frequencies to log2 scale
    tabdict = dict([(k,math.log(v,2)) for k,v in tabdict.iteritems() if v > 0])

    # if our frequency table no remaining entries, error
    if len(tabdict) == 0:
        return errorchar

    # if our frequency table has only one entry, we are done.
    if len(tabdict) == 1:
        return tabdict.keys()[0]

    # order nucleotides by increasing frequency
    sortedkeys = sorted(tabdict.keys(), key=lambda k: tabdict[k], reverse=True)

    # Collect the nucleotides that constitute a majority of samples at a site.
    # We want to emit ambiguity codes that summarize a position
    # while avoiding truly low-frequency bases.  An alignment motif would convey
    # more information, and perhaps that's what one should use, but a consensus
    # sequence must express each site in only a single letter so some loss of
    # information is inevitable.
    #
    # Our cheesy algorithm is to include all bases whose frequency is
    # at least half that of the next highest frequency nucleotide.
    nuc = [sortedkeys[0]]

    for i,v in enumerate(sortedkeys):
        if i == 0:
            continue
            
        # calculate log fold change between successive frequencies.
        logf = tabdict[sortedkeys[i - 1]]- tabdict[v]

        if logf <= 1:
            nuc.append(v)
        else:
            # otherwise stop looking.
            break
    
    if len(nuc) > 1:
        if use_ambi:
            return IUPAC_rev.get(tuple(sorted(nuc)), errorchar)
        else:
            return errorchar
    else:
        return sortedkeys[0]



    
def main(args=sys.argv[1:]):
    '''
    Calculate consensus sequences from a set of ancestral sequences output by BEAST

    The sequences are assumed to be created by beast sampling trait at particular nodes
    of the trees sampled from the posterior distribution.

    <state> "sequence1" "sequence2" "sequence3" ...

    state: the number of the iteration of the monte carlo simulation
    sequenceN: one point sample from the posterior distribution of sequences at a particular node in the tree.
    e.g.     200	"CCCAACGGATAATCGTAGCC" "CCGGACGGATAATCGTAGCC" "CCCGACGGATAATCGTAGCC"

    A consensus sequence will be calculated for each sequence on a line, with each line adding one sample to each consensus

    '''
    parser = build_parser()
    arguments = parser.parse_args(args)
    burnin = arguments.burnin
    
    # Keep a list of accumulators, one for each sequence that we run into.
    alignments = []	# collection of numpy matricies that accumulate results for each sequence
    labels = []			# sequence labels, found on the the first non-comment line in the file

    for line in arguments.infile:
        if line[0] == '#':
            continue
        fields = line.rstrip().split()
        
        # if the first field is not a comment and non-numeric,
        # assume it is the sequence labels.
        if  not is_int(fields[0]):
            if labels:
                print("warning: saw more than one line that looks like labels.", file=sys.stderr)
            labels = line.rstrip().split()[1:]
            continue

        # allow some burnin period before we start collecting data.
        if burnin > 0:
            burnin -= 1
            continue

        # sanity check - warn the user if all sequence appear identical.
        # this may mean a misconfigured BEAST config file (missing <mrca> tag in <ancestralTrait>?)
        if len(fields) > 2 and fields[2:] == fields[1:-1]:
            print("warning: all sequences are identical for state {}".format(fields[0]), file=sys.stderr)

        # iterate over each sequence on a line...
        for i,seq in enumerate(fields[1:]):
            if len(alignments) <= i:
                alignments.append(MultipleSeqAlignment([]))

            alignment = alignments[i]
            seq = seq.strip('"').rstrip('"')
            alignment.append(SeqRecord(Seq(seq, generic_dna), id=fields[0]))

    for i,v in enumerate(alignments):
        summary_align = AlignInfo.SummaryInfo(v)
        tabulated = tabulate( summary_align )
        consensus_str = ''.join([consensus(d, use_ambi=True) for d in tabulated]).upper()
        print(">{}\n{}".format(labels[i], consensus_str))


            
if __name__ == "__main__":
   main(sys.argv[1:])
   


