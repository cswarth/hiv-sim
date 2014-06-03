#!/usr/bin/python
# make a new BEAST  config file by inserting FASTA sequences into a generic template.
from __future__ import print_function

from lxml import etree
import re
from Bio import SeqIO

import sys, getopt


def main(argv):
    inputfile = ''
    outputfile = ''
    prefix = None
    try:
        opts, args = getopt.getopt(argv,"hp:",["prefix="])
    except getopt.GetoptError:
        print('mksanta.py [-p <prefix>] <templatefile> <fastafile>', file=sys.stderr)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('mksanta.py [-p <prefix>] <templatefile> <fastafile>', file=sys.stderr)
            sys.exit()
        elif opt in ("-p", "--prefix"):
            prefix = arg

    templatefile = args[0]
    datafile = args[1]

    # Parse a generic template and insert sequences from a FASTA file into the middle, separated by the appropriate XML tags.

    tree = etree.parse(templatefile)

    # process the fasta sequences to get their length.
    # they should all have the same length
    with open(datafile, "rU") as handle:
        ll = [len(record.seq) for record in  SeqIO.parse(handle, "fasta")]

    if (len(set(ll)) > 1):
        print("error - lengths are not concistent")
        exit(1)
        
    length = etree.Element("length")
    length.text = str(ll[0])

    # slurp in all the fasta format data
    with open(datafile, "rU") as handle:
        data = handle.read()

    sequences = etree.Element("sequences")
    sequences.text = data

    # Look for the genome tag
    # fill it in with sequence definitions
    genome = tree.find(".//genome")
    if (genome is None):
        print("Cannot find genome")
        print(etree.tostring(tree, pretty_print=True))
        exit(1)

    genome.clear()
    genome.append(length)
    genome.append(sequences)

    # if a prefix is supplied on the command line,
    # replace all occurances of the word 'patient' in the text of any
    # tag with the supplied prefix.
    if (prefix is not None):
        fixme = tree.xpath("//*[starts-with(text(), 'patient')]") 
        for e in fixme:
            e.text = e.text.replace('patient', prefix)

    print(etree.tostring(tree, pretty_print=True))


if __name__ == "__main__":
   main(sys.argv[1:])
   


