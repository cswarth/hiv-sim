#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    Convert a newick formatted tree into a nexus formatted tree.

    usage:

         newick2nexus.py output.newick input.nex 
'''
from __future__ import print_function
import dendropy

import sys
import argparse
import os.path

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('inputfile', nargs='?', help='input newick file', type=argparse.FileType('r'),
                     default=sys.stdin)
parser.add_argument('outputfile', nargs='*', help='output nexus file', type=argparse.FileType('w'),
                     default=sys.stdout)

a = parser.parse_args()

t = dendropy.Tree.get_from_stream(a.inputfile, "newick")

# Iterating Over Edges
# http://pythonhosted.org/DendroPy/tutorial/treenav.html#iterating-over-edges
t_len = t.length()
for edge in t.postorder_edge_iter():
    if edge.length is not None:
        edge.length = float(edge.length)/t_len 

t.write_to_stream(a.outputfile, "nexus", suppress_rooting=True)

