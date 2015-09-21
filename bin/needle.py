#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Calculate a one-to-many needleman-wunsch global alignment score. 

Typically used to calculate a score between a single founder sequence
and man inferred founder candidate sequences.  This script invokes the
`needle` program from the Emboss toolkit to do the actual alignment.

usage:

needle.py founder.fasta sequences.fasta
'''
from __future__ import print_function
import argparse
import os.path
from  itertools import ifilter
import  itertools
from Bio import SeqIO
import sys
import re
import subprocess
import tempfile

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
    parser.add_argument('founder', help='input fasta file',  type=existing_file)
    parser.add_argument('otherseqs', help='input fasta file', type=existing_file)
    parser.add_argument('--max', dest='max', action='store_true')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False, help="print more incremental progress information")
    parser.add_argument('-k', '--keep', help='keep temporary files, for debugging.', dest='keep', action='store_true', default=False)
    return parser.parse_args()


def needle_score(seq1, seq2, verbose=False, keep=False):
    """
    get needlman-wunsch score for aligning two sequences
    """
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq1', delete = not keep) as fh1, \
         ntf(prefix='seq2', delete = not keep) as fh2, \
         ntf(prefix='align_out') as outfile, \
         open(os.devnull) as dn:
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
        subprocess.check_call(cmd, stderr=dn)
        result = outfile.read()
        pattern = re.compile(r'# Score: (.*)')
        score = pattern.search(result)
        if score is not None:
            return float(score.group(1))
        return 0


class FastaIterator():
    '''
    leaks file descriptors and I'm too lazy to fix...
    probably the answer is to push filehandle management outside this class so the user can use 'with' contexts to manage their lifetime?
    '''
    
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

    
def main():
    a = parse_args()

    iseqs = FastaIterator(*a.otherseqs)
    if a.max:
        pattern = re.compile("^#(\d+)#$")
        root = None
        nroot = 0
        for seq in iseqs:
            n = pattern.match(seq.id)
            if n is not None:
                n = int(n.group(1))
                if n > nroot:
                    root = seq
                    nroot = n
        iseqs = iter([root])
        
    ifounder = FastaIterator(*a.founder)
    founder = ifounder.next()

    for seq in iseqs:
        score = needle_score(founder, seq, a.verbose, a.keep)
        print("{}-{}: {}".format(founder.id, seq.id, score))
            
    
        
if __name__ == '__main__':
    main()


