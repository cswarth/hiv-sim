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
import csv
import multiprocessing


# logging.basicConfig(stream=sys.stdout)
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
    parser.add_argument('-p', '--processes', default=15)
    return parser.parse_args()

def needle_records(fh):
    """
    separate output form needleman-wunsch alignment into individual records.

    We are only interested in the statistics asociated with each
    alignment, so discard all the details of specific base alignment.

    # Length: 2740
    # Identity:    2740/2740 (100.0%)
    # Similarity:  2740/2740 (100.0%)
    # Gaps:           0/2740 ( 0.0%)
    # Score: 13700.0

    Wish there were a way to make needle be quiet and just output the
    parts we need.  In fact it would be much more effeicient if it were
    just calculating the parts we are interested in, instead of wasting
    time and memory keep track of enough information constrct the actual
    alignment.  All I need is the score.
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


from collections import namedtuple
NWScore = namedtuple('NWScore', 'score identity gaps length')

def needle_score(seq1, seq2):
    """Calculate needlman-wunsch score for aligning two sequences.
    """
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq1', delete=False) as fh1, \
         ntf(prefix='seq2', delete=False) as fh2, \
         ntf(prefix='align_out', delete=True) as outfile:
        SeqIO.write(seq1, fh1, 'fasta')
        fh1.flush()
        SeqIO.write(seq2, fh2, 'fasta')
        fh2.flush()

        # invoke Needleman-Wunsch global alignment of two sequences from Emboss toolkit
        # http://emboss.sourceforge.net/apps/release/6.3/emboss/apps/needle.html
        # expect to find this in /home/matsengrp/local/bin/needle
        # NB the -gapopen and -gapextend parameters are NOT optional!
        # needle claims default values for these, but will not run if they are not specified.
        # same with -outfile.
        cmd = ['needle', '-gapopen', '10.0', '-gapextend', '0.5',  
               '-outfile',  outfile.name,
               fh1.name, fh2.name]
        logging.info(' '.join(cmd))
        subprocess.check_call(cmd, stderr=outfile)

        score_pattern = re.compile(r'# Score: (.*)')
        score_pattern = re.compile(r'# Length: (.*)')
        gaps_pattern = re.compile(r'# Gaps:\s+(\d+)/(\d+)')
        ident_pattern = re.compile(r'# Identity:\s+(\d+)/(\d+)')

        for record in needle_records(outfile):
            # Parse numeric statistics values from N-W text record.
            score = score_pattern.search(record)
            score = float(score.group(1)) if score else 0

            length = length_pattern.search(record)
            length = float(length.group(1)) if length else 0

            gaps = gaps_pattern.search(record)
            gaps = float(gaps.group(1))/float(gaps.group(2)) if gaps else 0
        
            identity = ident_pattern.search(record)
            identity = float(identity.group(1))/float(identity.group(2)) if identity else 0
            yield NWScore(score, identity, gaps, length)


def calculate_needle_score(founder, seq_iter):
    """
    Calculate various scores for comparing founder to sequences from 'seq_iter'

    :param founder: a single SeqIO.record object
    :param seq_iter: generator of SeqIO.Record objects
    :returns: tuple of scores
    :rtype: tuple of floats
    """
    # needle_score() is a generator function that will yield scores for all comparisons between 'founder' and
    # sequences in 'seq_iter'
    scores = list(needle_score(founder, seq_iter))

    # throw the scores into a dataframe and take the average of each column.
    df = pd.DataFrame.from_records(scores)
    return df.mean()

def parse_log(fp, burnin=0.9):
    """
    generator function for parsing ancestral founder sequences from BEAST trait file.

    The file expected to be a set of ancestral founder sequences from Beast.
    Each non-comment line consists of an iteration count followed by some number or quoted founder sequences.
    """

    # first count how many sequences are in the file so we can figure out how many to skip for burn-in.
    count = 0
    for line in csv.reader(fp, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        if not line[0].isdigit():
            continue
        count += 1

    skip = round(count * burnin)
    logging.debug("beast skipping {}/{} sequences for burning, keeping {}".format(skip,count, count-skip))
    
    # Reposition pointer at the beginning once again
    fp.seek(0, 0);
    for line in csv.reader(fp, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        if not line[0].isdigit():
            continue
        if skip > 0:
            skip -= 1
            continue
        yield line


def beast_iter(root, burnin):
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

    with open(beast, "rb") as fh:
        # NB parse_log() will skip a burn-in period
        for seq in parse_log(fh, burnin=burnin):
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


# NB this routine MUST be a top-level routine!
#
# Because this routine is caslled via multiprocessing.apply_async(), it must
# not be buried inside anothr routine or it won't be visible.  You might this
# this restiction would only affect Windows systems that lack true fork
# functionality, but apparently this restriction is true on all systems.
def _process_dir(founder, dir):
    """
    Calculate distance between founder & inferred sequences under directories 'dirs'

    Calculates distance measures between 'founder' and sequences
    inferred by several algorithms, including Prank and Beast.  For
    Prank files, only a single inferred founder is extracted.  For BEAST
    files we calculate scores for a number of high-posterior-probabilty
    founders and average the results.

    Higher scores are better.
"""
    (bscore, bidentity, bgaps, blength) = calculate_needle_score(founder, beast_iter(dir, 0.4))
    (pscore, pidentity, pgaps, plength)  = calculate_needle_score(founder, prank_iter(dir))
    cscore = calculate_control_score(founder, dir)
    sys.stdout.flush()
    # return a tuple 
    return (dir, cscore, pscore, pidentity, pgaps, int(plength), bscore, bidentity, bgaps, int(blength))

def process_founder(founder, dirs):
    """Calculate distance between founder & inferred sequences under directories 'dirs'

    This routine calculates distance measures between 'founder' and 
    sequences inferred by several algorithms, including Prank and
    Beast.

    For each inferred founder file,
    identify the root node and calculate the score between the root
    and the founder.  Higher scores are better.
    """
    pool = multiprocessing.Pool(15)
    results = [pool.apply_async(_process_dir, args = (founder, dir)) for dir in dirs]
    pool.close()
    pool.join()
    result_list = [r.get() for r in results]

    df = pd.DataFrame.from_records(result_list, columns=('root', 'control', 'p_score', 'p_identity', 'p_gaps', 'p_len', 'b_score', 'b_identity', 'b_gaps', 'b_len'))
    return(df)

def main():
    loglevel = "ERROR"
    a = parse_args()

    if a.verbose:
        loglevel = "INFO"
    if a.debug:
        loglevel = "DEBUG"

    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, stream=sys.stderr)

    logging.info("Info")
    logging.debug("Debug")

    tbl = None
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
            # These are directories at which donor and recipient samples are combined to infer a founder.
            dirs = [r for r, d, f in os.walk(root) if 'donor.fasta' in f]
            assert(len(dirs) != 0)
            # compare the founder sequence to all samples taken from downstream lineages 
            df = process_founder(founder, dirs)
            tbl = df if tbl is None else tbl.append(df, ignore_index=True)
            dirnames = []  # prune rest of tree.


    from datetime import datetime
    
    print("# Created by distances.py on {}".format(datetime.now().strftime("%Y-%m-%d %H:%M")))
    print("# 'score' refers to Needleman-Wunsch pairwise alignment score.")
    print("# 'score' refers to Needleman-Wunsch global alignment between the infered founder sequence and the actual founder.")
    print("# The Needleman-Wunsch algorithm is implemented in the Emboss toolkit.")
    print("# The \"EDNAFULL\" scoring matrix is used to score DNA comparisons with a gap open penalty of 10.0 and gap extension penalty of 0.5.")
    print("# p_score = score btwn prank inferred and actual founder")
    print("# p_gaps = # gap positions between prank inferred and actual founder")
    print("# p_identity = # identical sites between prank inferred and actual founder")
    print("# p_len = avg. length of prank inferred founder")
    print("# b_score = score btwn beast inferred and actual founder")
    print("# b_gaps = mean # gap positions between Beast high-posterior inferred and actual founder")
    print("# b_identity = avg. # sites between Beast high-posterior inferred and actual founder")
    print("# b_len = avg. length of beast inferred founder")
            
    tbl.to_csv(sys.stdout, index=False, header=True)
    
if __name__ == '__main__':
    main()


