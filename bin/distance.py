#!/usr/bin/env python
"""
Classify mutations for read groups in a BAM file,
using mutations identified by a GFF3 file
"""
import pkg_resources
# Require scipy version > 0.12.0 in order to use logsumexp with a second argument.
# http://docs.scipy.org/doc/scipy/reference/generated/logsumexp.html
pkg_resources.require("scipy>=0.12.0")

import contextlib
import argparse
import collections
import logging
from Bio import SeqIO
from itertools import izip   # for python3, just use built-in zip
import numpy as np

from scipy.misc import logsumexp
from datetime import datetime, date, timedelta
import getpass

import math
import sys

import pysam
import os.path

# parse a beast trait output file.
# this is expected to be a a state output file from beast holding inferred ancestral sequences.
def parse_log(fp):
    for line in fp:
        if line[0] == '#':
            continue
        if not line[0].isdigit():
            continue
        line = line.rstrip()
        tokens = line.split()
        yield tokens
        
# calculate distance between two sequences
# as the number of mismatched positions divided by the length.
def sitediff(s1, s2):
    return(np.sum([a!=b for a,b in zip(s1,s2)]))


def test():
    a = [0.5, 0.5, 0.5]
    print(a)
    a = [math.log(m) for m in a]
    print(a)
    b = [10.0, 10.0, 10.0]
    n = logsumexp(a, b=b)
    print(type(n))
    print("{}  {}".format(str(n), math.exp(n)))
    
def test2():
    founder = "ABC"
    seqs = [ "APC", "ADC"]
    b = np.asarray([ dist(founder, s) for s in seqs], dtype=np.float64)
    print(b)
    a = [0.5, 0.005]
    a = [np.log(n) for n in a]
    n = logsumexp(a, b=b)
    print(type(n))
    print("{}  {}".format(str(n), math.exp(n)))


# Convert a list of log probability values into a log probability distribution over those values
def log_normalize(x):
    max_val = np.max(x);
    x = x - max_val
    s = logsumexp(x)
    x = x - s
    return(x)

# convert version string to tuples that can be compared
def versiontuple(v):
    return tuple(map(int, (v.split("."))))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('founder', metavar="founder.fasta", type=argparse.FileType('r'))
    p.add_argument('seqs', metavar='ancestral.log', type=argparse.FileType('r'))
    p.add_argument('post', metavar='posterior.log', type=argparse.FileType('r'))
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
                   default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s: %(message)s')

    # Read in the founder sequence.  We only read the first sequence in the fasta file and assume that is the
    # founder.  Additional sequences are ignored!
    founder = None
    with contextlib.closing(a.founder), a.founder as fp:
        for rec in SeqIO.parse(fp, "fasta"):
            if founder is None:
                founder = str(rec.seq)
            else:
                print("Warning: Additional founder sequences ignored in {}".format(a.founder.name))
            break

    with a.seqs as seqfp:
        changes = np.asarray([ sitediff(founder, seq[1][1:-1]) for seq in parse_log(seqfp)], dtype=np.float64)

    with a.post as postfp:
        post = np.asarray([ post[1] for post in parse_log(postfp)], dtype=np.float64)

    
    post = log_normalize(post)
    dist = changes/len(seq[1])
    wdist = logsumexp(post, b=dist)

    '''
    min = minumim number of mismatched between sequence and founder
    max = maximum number of mismatched between sequence and founder
    mean = mean number of mismatched between sequence and founder
    wdist = log expected distance between sequence and founder - weighted by normalized posterior probability
    exp.wdist = exp(wdist)
    '''

    print("# date: {}".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    print("# user: {}".format(getpass.getuser()))
    print("# command: {}".format(" ".join(sys.argv)))
    print("# workingdir: {}".format(os.getcwd()))
    print("min\tmax\tmean\twdist\texp.wdist")
    print("{}\t{}\t{:.3}\t{:.5}\t{:.5}".format(np.int(np.min(changes)), np.int(np.max(changes)), np.mean(changes), wdist, np.exp(wdist)))

if __name__ == '__main__':
    main()
       
