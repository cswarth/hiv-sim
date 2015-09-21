#!/usr/bin/env python
'''
Extracts and prints  ancestral sequence at the root of the phylogenetc tree inferred by PRANK.

usage:
	$ ./prankroot.py prank.best.anc.dnd prank.best.anc.fas 
'''
from __future__ import print_function
import sys
import argparse
import logging
from Bio import SeqIO, Alphabet
from Bio import Phylo	# for reading Prank trees

# logging.basicConfig(stream=sys.stdout)
log = logging.getLogger(__name__)
a = None # reserved for arguments


def prank_root_sequences(treefile, fastafile):
    """
    iterate over Prank founder sequence(s) inferred from CODON model
    See https://github.com/cswarth/hiv-sim/issues/2

    This routine knows how to retrieve and iterate over sequences
    produced by Prank. It will yield a single sequence corresponding to
    the root of the guide tree.

    :param root: string name of directory where sequence files can be found.
    :returns: iterator over seqeuences
    :rtype:
    """

    # identify the sequence associated with the root node.
    tree = Phylo.read(treefile, 'newick')
    rootid = str(tree.root)
    record_dict = SeqIO.index(fastafile, "fasta")
    return record_dict[rootid] 

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
    parser.add_argument('-p', '--processes', default=2)
    parser.add_argument('treefile')
    parser.add_argument('fastafile')
    return parser.parse_args()

def main():
    loglevel = "ERROR"
    global a
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

    print(prank_root_sequences(a.treefile, a.fastafile).seq)


    
if __name__ == '__main__':
    main()


    
