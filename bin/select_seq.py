#!/usr/bin/env python

"""                                                                                 
This script to selects a random sample of sequences without replacement from a particular generation.
The sequence ids are presumed to have the format <string>_<generation>_<string>
The selected sequennce(s) will be selected from those with an id matching the specified generation.
"""
# inspired by https://www.biostars.org/p/1709/#21500
                                                               
from Bio import SeqIO                                                               
import random                                                                          
import argparse                                                                          
import re                                                                          
import sys                                                                          
                                                                                    
p = argparse.ArgumentParser()
p.add_argument('fasta', metavar="founder.fasta", type=argparse.FileType('r'), help="fasta format file from which sequences are to be selected")
p.add_argument('generation', nargs='?', type=int, default=100, help='genertion to choose')
p.add_argument('count', nargs='?', type=int, default=1, help='number of sample sequences to extract')
a = p.parse_args()

pat = re.compile("[^_]+_{}_".format(a.generation))

seqiter = SeqIO.parse(a.fasta, 'fasta')                                    
SeqIO.write(random.sample([seq for seq in seqiter if re.match(pat, seq.id)], a.count), sys.stdout, "fasta")
