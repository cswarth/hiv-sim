#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Calculate mean Ka/Ks ratio between founder sequences and descendants.

Usage:
    cd ~/src/matsen/hiv-sim/kaks
    . venv/bin/activate
    ../bin/paml.py noselection.fa mcc.tree

'''

from __future__ import print_function

from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML import baseml, yn00

from itertools import izip, chain, permutations, starmap
import operator
import numpy as np
import sys
import argparse
import os.path
from Bio import SeqIO, Phylo
import dendropy
import json

def codeml_action(args):
    workdir = './scratch'
    tmpaln = os.path.join(workdir, "align.fa")
    tmptree = os.path.join(workdir, "tree.nex")
    
    # preprocess the alignment file to remove sequences with stop codons
    with open(args.alignment, "rU") as fhin, open(tmpaln, "w") as fhout:
        for record in SeqIO.parse(fhin, "fasta"):
            if not '*' in record.seq.translate():
                record.id = record.id.replace('/','_')
                # codeml seems unhappy when sequences have descriptions
                # maybe it's only certain characters in descriptions?
                record.description = ""
                print(record.id)
                SeqIO.write(record, fhout, "fasta")


    # convert the newick tree to unrooted nexus format.
    tree = Phylo.read(args.tree, 'nexus')
    for t in tree.get_terminals():
        t.name = t.name.replace('/','_')

    Phylo.write([tree], tmptree, 'newick')

    cml = codeml.Codeml()
    cml.working_dir = "scratch"
    cml.alignment = tmpaln
    cml.tree = tmptree
    cml.out_file = "results.out"
    
    cml.set_options(noisy = 9) # 0,1,2,3,9: how much rubbish on the screen
    cml.set_options(verbose = (args.verbosity > 1)) # 1: detailed output, 0: concise output
    cml.set_options(runmode = -2)   # 0:user tree;  1:semi-automatic;  2:automatic
                      # 3:StepwiseAddition; (4,5):PerturbationNNI 

    cml.set_options(seqtype = 1) # 1:codons; 2:AAs; 3:codons-->AAs
    cml.set_options(CodonFreq = 2) # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    cml.set_options(clock = 0) # 0: no clock, unrooted tree, 1: clock, rooted tree
    cml.set_options(aaDist = 0) # 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
    cml.set_options(model = 0)

    cml.set_options(NSsites = [0])
    # 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
    # 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
    cml.set_options(icode = 0) # 0:standard genetic code; 1:mammalian mt; 2-10:see below
    cml.set_options(Mgene = 0) # 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    cml.set_options(fix_kappa = 0) # 1: kappa fixed, 0: kappa to be estimated
    cml.set_options(kappa = 1) # initial or fixed kappa
    cml.set_options(fix_omega = 0) # 1: omega or omega_1 fixed, 0: estimate
    cml.set_options(omega = 1) # initial or fixed omega, for codons or codon-based AAs
    cml.set_options(ncatG = 10) # # of categories in the dG or AdG models of rates

    cml.set_options(getSE = 0) # 0: don't want them, 1: want S.E.s of estimates
    cml.set_options(RateAncestor = 0) # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
    cml.set_options(Small_Diff = .45e-6)
    cml.set_options(cleandata = 1) # remove sites with ambiguity data (1:yes, 0:no)?
    cml.set_options(fix_blength = 0) # 0: ignore, -1: random, 1: initial, 2: fixed


    cml.print_options()
    cml.ctl_file = "control2.ctl"
    print(cml.ctl_file)
    cml.write_ctl_file()

    results = cml.run(verbose=(args.verbosity > 1))

    # print out the parsed results dict 
    if args.verbosity > 0:
        print(json.dumps(results, sort_keys=True, indent=4))

    return results

# Retrieve values form multidimensional dict
# Given a dict like d = {a: { b: { c: 3 } } }
# chained_get(d, 'a', 'b', 'c') will return 3
# http://stackoverflow.com/a/16003666/1135316
def chained_get(dct, *keys):
    SENTRY = object()
    def getter(level, key):
        return None if level is SENTRY else level.get(key, SENTRY)
    return reduce(getter, keys, dct)

'''
Retrieve values from hierarchical dict representing a pairwise comparison.
Given a hierarchical dict like this,

    "patient1_1000_1|1M|XXX|XXX|201": {
        "patient1_1000_27|1M|XXX|XXX|20": {
            "LPB93": {
                "dN": 0.125, 
                "dS": 0.045, 
                "w": 2.7905
            }, ...
        }
    }

this routine can retrieve a named value for each of the pairwise comparisons.
get_pairwise(d, 'LBP93', 'w') would return a list [2.7905 ...]
containing every 'w' value below 'LPB93' two levels down in the dict.
'''    
def get_pairwise(results, *names):
    return [chained_get(v2,*names) for p1,v in results.items() for p2,v2 in v.items()]



def yn00_action(args):
    workdir = './scratch'
    tmpaln = os.path.join(workdir, "align.fa")
    
    # Remove sequences with stop codons
    # yn00 will fail if there are stop codons in any of the sequences.
    n = 0
    with open(args.alignment, "rU") as fhin, open(tmpaln, "w") as fhout:
        for record in SeqIO.parse(fhin, "fasta"):
            if not '*' in record.seq.translate():
                n += 1
                SeqIO.write(record, fhout, "fasta")
    if n >= 2:
        yn = yn00.Yn00()
        yn.working_dir = "scratch"
        yn.alignment = tmpaln
        yn.out_file = "ynresults.out"
        yn.set_options(verbose = (args.verbosity > 1)) # 1: detailed output, 0: concise output

        results = yn.run(verbose=(args.verbosity > 1))

        # results = yn00.read('ynresults.out')
    else:
        results = {}

    # print out the parsed results dict 
    if args.verbosity > 0:
        print(json.dumps(results, sort_keys=True, indent=4))
        
    w = [w for w in get_pairwise(results, 'NG86', 'omega') if w > 0]
    if len(w) < 2:
        w = -1
    else:
        w = np.mean(w)
        
    print(w)

    sys.exit(0)

def parse_arguments(argv):
    """
    Build the command-line argument parser.
    """

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-v', '--verbose', dest='verbosity',
            action='count', default=0,
            help="Be more verbose. Specify -vv or -vvv for even more")

    # Subparsers
    # specification of subparsers cribbed from `seqmagick`
    # https://github.com/fhcrc/seqmagick/blob/master/seqmagick/scripts/cli.py
    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_help = subparsers.add_parser('help',
            help='Detailed help for actions using help <action>')

    parser_help.add_argument('action')

    # Add actions
    actions = {}
    subparser = subparsers.add_parser('codeml', help='PAML codeml',
            description='PAML codeml')
    subparser.add_argument('alignment', help='sequences in fasta format', type=existing_file)
    subparser.add_argument('tree', help='tree file in newick format', type=existing_file)
    actions['codeml'] = codeml_action
    
    subparser = subparsers.add_parser('yn00', help='PAML yn00',
            description='PAML yn00')
    subparser.add_argument('alignment', help='sequences in fasta format', type=existing_file)
    actions['yn00'] = yn00_action

    arguments = parser.parse_args(argv)
    arguments.argv = argv
    action = arguments.subparser_name

    if action == 'help':
        return parse_arguments([str(arguments.action), '-h'])

    return actions[action], arguments


def main(argv=sys.argv[1:]):
    '''
    Calculate mean ka/ks ration (aka dn/ds) between known founder and descendant sequences.
    '''
    action, arguments = parse_arguments(argv)

    return action(arguments)


    
if __name__ == "__main__":
   main(sys.argv[1:])
   



