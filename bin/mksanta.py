#!/usr/bin/python
# make a new SANTA config file by inserting FASTA sequences into a generic template.
from __future__ import print_function

from lxml import etree
import re
from Bio import SeqIO
from datetime import datetime, date, timedelta

import sys, getopt


# indent xml text for pretty-printing
# assumes text within elements is not significant.
def indent(elem, level=0):
    i = "\n" + level*"\t"
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "\t"
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def main(argv):
    inputfile = ''
    outputfile = ''
    data = ''
    datestr = None
    prefix = None
    usage = 'mksanta.py [-p <prefix>] [-d sampledate] <templatefile> <source> <generation> <count>'
    try:
        opts, args = getopt.getopt(argv,"hp:d:",["prefix=", "date="])
    except getopt.GetoptError:
        print(usage, file=sys.stderr)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage, file=sys.stderr)
            sys.exit()
        elif opt in ("-p", "--prefix"):
            prefix = arg
        elif opt in ("-d", "--date"):
            # if the user supplies a datestr on the commandline then we will not
            # use the starting date found in the source file
            datestr = arg




    argTemplate = args[0]
    argSourcefile = args[1]
    argGeneration = args[2]
    argCount = args[3]

    count = int(argCount)
    # Parse a generic template and insert sequences from a FASTA file
    # into the middle, separated by the appropriate XML tags.

    tree = etree.parse(argTemplate)

        
    name_regex = re.compile("(?P<patient>^[^_]*)_(?P<generation>.*)_")

    # find the sequences from the chosen <generation> and
    # randomly choose <count> of them to pass on to the next generation
    with open(argSourcefile, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta") :
            # import pdb; pdb.set_trace()

            fields = record.id.split('|')
            match = name_regex.search(fields[0])
            patient = match.group('patient')
            generation = match.group('generation')

            if generation == argGeneration:
                if count > 0:
                    # Randomly choose a sequence.  For now
                    # assume they are already randomy sorted and choose from the
                    # order in which they are presented.
                    #
                    # need to assert the all the lenghts are the same....
                    # if not, we will have to align all the sequences to one another...
                    ll = len(record.seq)
                    data += ">{}\n".format(record.id)
                    data += "{}\n".format(record.seq)

                    sampledGeneration = int(generation)
                    cumulativeGeneration  = int(fields[3])
                    count -= 1

            if count <= 0:
                # break out when we have collected all the sequences we need
                break

    sequences = etree.Element("sequences")
    sequences.text = data

    length = etree.Element("length")
    length.text = str(ll)


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

    if (prefix is None):
        prefix = "patient"

    # preprend the to the fasta output file name and the tag on each sequence.
    fixme = tree.xpath("//*[starts-with(text(), 'patient')]") 
    for e in fixme:
        e.text = e.text.replace('patient', prefix)

    # compose a label to indicate where the starting population from the sample
    label = tree.find("//samplingSchedule/sampler/alignment/label")
    cumulativeGeneration += sampledGeneration
    label.text = "{}_%g_%s|{}|{}|{}".format(prefix, argSourcefile, sampledGeneration, cumulativeGeneration)

    # pretty-print the tree
    indent(tree.getroot())
      
    print(etree.tostring(tree, pretty_print=True))


if __name__ == "__main__":
   main(sys.argv[1:])
   


