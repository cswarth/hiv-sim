#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    Convert a nexus formatted tree into a newick formatted tree.
	The output is specific for PRANK whic is very sensitive to the format of the tree.

    usage:

         nexus2newick.py input.nex output.newick
'''
from __future__ import print_function
import dendropy

import sys
import argparse
import os.path

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('inputfile', nargs='?', help='input nexus file', type=argparse.FileType('r'),
                     default=sys.stdin)
parser.add_argument('outputfile', nargs='*', help='output newick file', type=argparse.FileType('w'),
                     default=sys.stdout)

a = parser.parse_args()

t = dendropy.Tree.get_from_stream(a.inputfile, "nexus", preserve_underscores=True)

# Iterating Over Edges
# http://pythonhosted.org/DendroPy/tutorial/treenav.html#iterating-over-edges
t_len = t.length()
for edge in t.postorder_edge_iter():
    if edge.length is not None:
        edge.length = float(edge.length)/t_len 

t.write_to_stream(a.outputfile, "newick", suppress_rooting=True)

