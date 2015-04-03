#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Uses 'needle' program from Emboss package
(http://emboss.sourceforge.net/download/) to measure distance between
actual founder sequence and inferred founder sequence. the distace is
measured as the needleman-munsch score for a pairwise alignment.  We
don't hae fast N-W alignment software accessible in Python so we
actually call the 'needle'' program and parse the results to extract the
score we need.

usage:
    cd ~/src/matsen/hiv-sim/sims
    ../bin/distance.py 
'''
from __future__ import print_function
import argparse
import os.path
from  itertools import ifilter
import  itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import re
import subprocess
import tempfile
import fnmatch
import numpy as np
import logging

def parse_args():
    ''' returns command-line arguments parsed by argparse.
        >>> args = parse_args()
    '''

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
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False)
    return parser.parse_args()


def needle_score(seq1, seq2, verbose=False, keep=False):
    """Calculate needlman-wunsch score for aligning two sequences.

    Takes two Bio.Seq objects.  The second sequence may be just a python string, in which case an anonymous 
    Bio.Seq object will be create automatically.  
    (really?  that seems silly to put here.  have the caller create a Seq object if they need to... - csw)
    """
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq1', delete = not keep) as fh1, \
         ntf(prefix='seq2', delete = not keep) as fh2, \
         ntf(prefix='align_out') as outfile:
        SeqIO.write(seq1, fh1, 'fasta')
        fh1.flush()
        SeqIO.write(seq2, fh2, 'fasta')
        fh2.flush()

        cmd = ['needle', '-gapopen', '0',
               '-gapextend', '0',
               '-outfile',  outfile.name,
               fh1.name, fh2.name]
        if verbose:
            print(' '.join(cmd))
        subprocess.check_call(cmd, stderr=outfile)
        result = outfile.read()
        pattern = re.compile(r'# Score: (.*)')
        score = pattern.search(result)
        if score is not None:
            return float(score.group(1))
        return 0


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
    count = 2
    for line in fp:
        if line[0] == '#':
            continue
        if not line[0].isdigit():
            continue
        line = line.rstrip()
        tokens = line.split()
        count -= 1
        if count == 0:
            break
        yield tokens
        

def process_founder(root):
    """Calculate the distance from founder sequence to each founder sequences inferred by PRANK.

    There will be multiple simulated transmission hierarchies below the founder sequence file.
    For each inferred founder file, identify the root noe and calculate the score betweenthe root and the founder.
    Higher scores are better.
    """
    with open(os.path.join(root, "founder.fa")) as fh:
        founder = SeqIO.read(fh, "fasta")

    for root, dirnames, filenames in os.walk(root):
        if 'donor.fasta' in filenames:
            pscore = calculate_prank_score(founder, root)
            bscore = calculate_beast_score(founder, root)
            cscore = calculate_control_score(founder, root)
            print("{}\t{}\t{}\t{}".format(root, pscore, bscore,cscore))
            dirnames = []  # prune the tree below this point

def calculate_prank_score(founder, root):
    pscore = None
    # identify the root node.  It is the one with the highest numeric label
    pattern = re.compile(r'^#(\d+)#$')
    with open(os.path.join(root, 'prank.best.anc.fas')) as fh:
        rootnode = None
        rootnum = 0
        iseqs = SeqIO.parse(fh, "fasta")
        for seq in iseqs:
            n = pattern.match(seq.id)
            if n is not None:
                n = int(n.group(1))
                if n > rootnum:
                    rootnode = seq
                    rootnum = n
        pscore = needle_score(founder, rootnode)
    return pscore
    
def calculate_beast_score(founder, root):
    def str2Seq(s):
        """Our distance measure takes a Bio.Seq object,
        but we will be reading strings from the beast ancestror log.
        Define a routine to do the conversion for us.
        """
        return SeqRecord(Seq(s, IUPAC.IUPACUnambiguousDNA),
                  id="dynamic", name="anonymous",
                  description="dynamically created sequence form string")

    bscore = None
    beast = os.path.join(root, 'ancestralSequences.log')
    with open(beast) as fh:
        # parse_log() returns a list of tokens.
        # sequences begin at the second token in the list
        scores = [ needle_score(founder, str2Seq(seq[1])) for seq in parse_log(fh)]
        bscore = np.mean(scores)
    return bscore

def calculate_control_score(founder, root):
    # cscore is a control value that increases monotonically with donor and recipient distance.
    # it is used to validate plotting methods.
    # for now have it indicate the total distance between samples, so just the sum of donor and recipient dimes since infection.
    def is_number(s):
        """return true if argument looks like a python number, false otherwise."""
        try:
            float(s)
            return True
        except ValueError:
            return False
    tokens = root.split("/")
    donor = int(tokens[-3]) if is_number(tokens[-3]) else 0
    recipient = int(tokens[-2]) if is_number(tokens[-2]) else 0
    cscore = donor+recipient

    return cscore
    
def main():
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s:%(name)s: %(message)s')
    a = parse_args()

    if a.verbose:
        Logger.setLevel(logging.INFO)

    print("{}\t{}\t{}\t{}".format("dir", "prank", "beast", "control"))

    matches = []
    for root, dirnames, filenames in os.walk('runs'):
        if 'founder.fa' in filenames:
            logging.info('found "founder.fa" in {}'.format(root))
            process_founder(root)
            dirnames = []  # prune rest of tree.

      
    
        
if __name__ == '__main__':
    main()


