
import os.path
import sys

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
        
 
