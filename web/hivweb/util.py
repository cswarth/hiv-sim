
import os.path
import sys
import itertools

from Bio import SeqIO




def highlight(seqs):
    """Convert columns of identical residues into '.' characters
    		ACGGCTA			    		.C.G...
    		ACGCCTA		becomes >>>		.C.C...
    		AGGGCTA			    		.G.G...

    """
    def dots(col):
        if len(set(col)) == 1:
            return "."*len(col)
        else:
            return col
    print(type(zip(*map(dots,zip(*seqs)))))

    return map(''.join, zip(*map(dots,zip(*seqs))))
        

def is_number(s):
    """return true if argument looks like a python number, false otherwise.

    Used when sorting list of mixed numeric and string values.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


# Fast min/max function
# Does only a single pass so works with iterators as well as lists
# http://code.activestate.com/recipes/577916-fast-minmax-function/
def minmax(data):
    """Computes the minimum and maximum values in one-pass using only 1.5*len(data) comparisons"""
    it = iter(data)
    try:
        lo = hi = next(it)
    except StopIteration:
        raise ValueError('minmax() arg is an empty sequence')
    for x, y in itertools.izip_longest(it, it, fillvalue=lo):
        if x > y:
            x, y = y, x
        if x < lo:
            lo = x
        if y > hi:
            hi = y
    return lo, hi


# stolen from http://code.activestate.com/recipes/410692/
# This class provides the functionality we want. You only need to look at
# this if you want to know how this works. It only needs to be defined
# once, no need to muck around with its internals.
class switch(object):
    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration
    
    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args: # changed for v1.5, see below
            self.fall = True
            return True
        else:
            return False

