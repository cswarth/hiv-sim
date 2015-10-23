#!/usr/bin/env python
'''
Split the Prank trees at the root 

     usage: ./treesplit.py output/prank.best.anc.dnd samples/RV217_PDC_gp120_1M_aln.fa
     writing left.fasta
     RV217_PDC|1M|08WG|NFLG|2012/02/12
     writing right.fasta
     RV217_PDC|1M|10WG|NFLG|2012/02/12 RV217_PDC|1M|11WG|NFLG|2012/02/12 RV217_PDC|1M|09WG|NFLG|2012/02/12 RV217_PDC|1M|07WG|NFLG|2012/02/12 RV217_PDC|1M|06WG|NFLG|2012/02/12 RV217_PDC|1M|05WG|NFLG|2012/02/12 RV217_PDC|1M|04WG|NFLG|2012/02/12 RV217_PDC|1M|01WG|NFLG|2012/02/12 RV217_PDC|1M|02WG|NFLG|2012/02/12 RV217_PDC|1M|03WG|NFLG|2012/02/12

'''
from __future__ import print_function
import sys
import os.path
import argparse
import logging
from Bio import SeqIO, Alphabet
from Bio import Phylo	# for reading Newick trees like those produced by Prank

# logging.basicConfig(stream=sys.stdout)
log = logging.getLogger(__name__)
a = None # reserved for arguments


def split_sequences(treefile, fastafile, outdir):
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
    for f,child in zip(['left','right'], tree.root):
        fname = os.path.join(outdir, f+'.fasta')
        with open(fname, "w") as fh:
            sequences = [record_dict[k.name] for k in child.get_terminals()]
            SeqIO.write(sequences, fh, "fasta")


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
    parser.add_argument('-o', '--outdir', default='.')
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

    split_sequences(a.treefile, a.fastafile, a.outdir)


    
if __name__ == '__main__':
    main()


    
