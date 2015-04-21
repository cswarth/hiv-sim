#!/usr/bin/env python
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
    ../bin/distances.py 
'''
from __future__ import print_function
import argparse
import os.path
from  itertools import ifilter
import  itertools
from Bio import SeqIO
from Bio import Phylo	# for reading Prank trees
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import re
import subprocess
import tempfile
import fnmatch
import pandas as pd
import numpy as np
import logging


logging.basicConfig(stream=sys.stdout)
log = logging.getLogger(__name__)

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
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False)
    return parser.parse_args()

def needle_records(fh):
    """separate output of needleman-wunsch alignment into individual records.

    We are only interested in the statistics asociated with each alignment, so discard all the 
    details of specific base alignment.  Wish there were a way to make needle be quiet and just
    output the parts we need.
    """
    recordsep = '#======================================='
    for line in fh:
        if line.startswith(recordsep):
            record = ''
            for line in fh:
                if line.startswith(recordsep):
                    yield record
                    break
                record += line
    
def needle_score(seq1, seq2):
    """Calculate needlman-wunsch score for aligning two sequences.
    """
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq1', delete=False) as fh1, \
         ntf(prefix='seq2', delete=False) as fh2, \
         ntf(prefix='align_out', delete=False) as outfile:
        SeqIO.write(seq1, fh1, 'fasta')
        fh1.flush()
        SeqIO.write(seq2, fh2, 'fasta')
        fh2.flush()

        # invoke Needleman-Wunsch global alignment of two sequences from Emboss toolkit
        # http://emboss.sourceforge.net/apps/release/6.3/emboss/apps/needle.html
        # expect to find this in /home/matsengrp/local/bin/needle
        # NB the -gapopen and -gapextend parameters are NOT optional!
        # needle claims default values for these, but will not run if they are not specified.
        cmd = ['needle', '-gapopen', '10.0', '-gapextend', '0.5',  
               '-outfile',  outfile.name,
               fh1.name, fh2.name]
        logging.info(' '.join(cmd))
        subprocess.check_call(cmd, stderr=outfile)

        score_pattern = re.compile(r'# Score: (.*)')
        gaps_pattern = re.compile(r'# Gaps:\s+(\d+)/(\d+)')
        ident_pattern = re.compile(r'# Identity:\s+(\d+)/(\d+)')

        for record in needle_records(outfile):
            score = score_pattern.search(record)
            score = float(score.group(1)) if score is not None else 0

            gaps = gaps_pattern.search(record)
            length = int(gaps.group(2)) if gaps is not None else 0
            gaps = float(gaps.group(1))/float(gaps.group(2)) if gaps is not None else 0
        
            identity = ident_pattern.search(record)
            identity = float(identity.group(1))/float(identity.group(2)) if identity is not None else 0
            yield (score, gaps, identity, length)


def calculate_needle_score(founder, seq_iter):

    df = pd.DataFrame(columns=('score', 'identity', 'gaps', 'len') )
    x = 0
    for scores in needle_score(founder, seq_iter):
        df.loc[x] = scores
        x += 1
    logging.debug(df.shape)
    logging.debug(df.mean())

    return df.mean()


def parse_log(fp, burnin=0.9):
    """
    parse a beast ancestral sequence output file.
    this is expected to be a state output file from beast holding inferred ancestral sequences.
    Each non-comment line  consists of a number, 
    """
    count = 0
    for line in fp:
        if line[0] == '#':
            continue
        if not line[0].isdigit():
            continue
        count += 1

    skip = round(count * burnin)
    # Reposition pointer at the beginning once again
    fp.seek(0, 0);
    if skip > 0:
        for line in fp:
            if skip <= 0:
                break
            if line[0] == '#':
                continue
            if not line[0].isdigit():
                continue
            skip -= 1

    for line in fp:
        if line[0] == '#':
            continue
        if not line[0].isdigit():
            continue
        line = line.rstrip()
        tokens = line.split()
        yield tokens
        

def beast_iter(root):
    """
    iterate over BEAST sequences inferred founder sequences

    This routine knows how to retrieve and iterate over sequences
    produced by Beast.  It yields a series of sequences corresponding to
    high posterior probability after a 10% burn-in period.

    :param root: string name of directory where sequence files can be found.
    :returns: iterator over seqeuences
    :rtype:
    """

    def str2Seq(s, id='dynamic'):
        """Our distance measure takes a Bio.Seq object,
        but we will be reading strings from the beast ancestror log.
        Define a routine to do the conversion for us.
        """
        return SeqRecord(Seq(s, IUPAC.IUPACUnambiguousDNA),
                  id=id, name="anonymous",
                  description="dynamically created sequence from ancestralSequences.log")
    
    beast = os.path.join(root, 'ancestralSequences.log')
    logging.debug("iterating beast {}".format(beast))

    with open(beast) as fh:
        # NB parse_log() will skip a burn-in period
        for seq in parse_log(fh):
            yield str2Seq(seq[1], id=seq[0])


def prank_iter(root):
    """
    iterate over Prank inferred founder sequence(s)

    This routine knows how to retrieve and iterate over sequences
    produced by Prank. It will yield a single sequence corresponding to
    the root of the guide tree.

    :param root: string name of directory where sequence files can be found.
    :returns: iterator over seqeuences
    :rtype:
    """

    # identify the sequence associated with the root node.
    tree = Phylo.read(os.path.join(root, 'prank.best.anc.dnd'),'newick')
    rootid = str(tree.root)
    record_dict = SeqIO.index(os.path.join(root, 'prank.best.anc.fas'), "fasta")
    yield record_dict[rootid] 

    
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
    donor = int(tokens[-2]) if is_number(tokens[-2]) else 0
    recipient = int(tokens[-1]) if is_number(tokens[-1]) else 0
    cscore = donor+recipient
    return cscore

def process_founder(founder, dirs):
    """Calculate distance between founder & inferred sequences under directories 'dirs'

    This routine calculates distance measures between 'founder' and 
    sequences inferred by several algorithms, including Prank and
    Beast.

    For each inferred founder file,
    identify the root node and calculate the score between the root
    and the founder.  Higher scores are better.
    """

    # preallocate the dataframe for efficiency
    df = pd.DataFrame(columns=('root', 'control',
                               'p_score', 'p_identity', 'p_gaps', 'p_len',
                               'b_score', 'b_identity', 'b_gaps', 'b_len'))

    df.to_csv(sys.stdout, index=False, header=True)
    for dir in dirs:
        (bscore, bidentity, bgaps, blength) = calculate_needle_score(founder, beast_iter(dir))
        (pscore, pidentity, pgaps, plength)  = calculate_needle_score(founder, prank_iter(dir))
        cscore = calculate_control_score(founder, dir)
        df.loc[0] = (dir, cscore, pscore, pidentity, pgaps, int(plength), bscore, bidentity, bgaps, int(blength))
        df.to_csv(sys.stdout, index=False, header=False)
        sys.stdout.flush()

    return(df)

def main():
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s:%(name)s: %(message)s')
    a = parse_args()

    if a.verbose:
        log.setLevel(logging.INFO)
    if a.debug:
        log.setLevel(logging.DEBUG)
    logging.info("Info")
    logging.debug("Debug")

    tbl = pd.DataFrame()
    matches = []
    for root, dirnames, filenames in os.walk('runs'):
        if 'founder.fa' in filenames:
            logging.debug(root)
            logging.info('found "founder.fa" in {}'.format(root))

            # read the single founder sequence that initiated this lineage.
            with open(os.path.join(root, "founder.fa")) as fh:
                founder = SeqIO.read(fh, "fasta")
            assert(founder is not None)
            assert(len(founder) > 1)

            # Find all directories below this point holding a 'donor.fasta' file
            # These are points at which donor and recipient samples are combined to infer a founder.
            dirs = [r for r, d, f in os.walk(root) if 'donor.fasta' in f]
            assert(len(dirs) != 0)
            
            process_founder(founder, dirs)
            dirnames = []  # prune rest of tree.

        
if __name__ == '__main__':
    main()


