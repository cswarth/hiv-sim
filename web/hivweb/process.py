from cStringIO import StringIO
import gzip
import os.path
import sys
import re
import subprocess
import tempfile

from Bio import SeqIO

# Calculate Needleman-Wunsch global alignment score between two sequences.
#
# I don't think this routie is actually used right now
# as the alignment score is available as part of the precomputed distance table
# passed on the commandline.
#
# Author: cswarth
def needle_score(seq1, seq2, verbose=False, keep=False):
    """
    get needlman-wunsch score for aligning two sequences
    """
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq1', delete = not keep) as fh1, \
         ntf(prefix='seq2', delete = not keep) as fh2, \
         ntf(prefix='align_out') as outfile:
        SeqIO.write(seq1, fh1, 'fasta')
        fh1.flush()
        SeqIO.write(seq2, fh2, 'fasta')
        fh2.flush()

        # invoke Needleman-Wunsch global alignment of two sequences from Emboss toolkit
        # http://emboss.sourceforge.net/apps/release/6.3/emboss/apps/needle.html
        # expect to find this in /home/matsengrp/local/bin/needle
        

        cmd = ['needle', '-outfile',  outfile.name,
               fh1.name, fh2.name]
        if verbose:
            print(' '.join(cmd))
        subprocess.check_call(cmd, stderr=subprocess.STDOUT)
        result = outfile.read()
        pattern = re.compile(r'# Score: (.*)')
        score = pattern.search(result)
        if score is not None:
            return float(score.group(1))
        return 0

# Author: cswarth
def multi_align(*seqs):
    """return multiple sequence alignment of 

    Spawns `muscle` multiple sequence alignment process 

    :returns: aligned sequences
    :rtype: list of SeqIO records

    """
    results = None
    ntf = tempfile.NamedTemporaryFile
    with ntf(prefix='seq', delete=True) as fh, \
      	 open(os.devnull) as dn, \
         ntf(prefix='align_out') as outfile:
         for seq in seqs:
             SeqIO.write(seq, fh, 'fasta')
	 fh.flush()
	 
         cmd = ['muscle', '-quiet', '-out',  outfile.name, '-in', fh.name]
	 print(' '.join(cmd))
         subprocess.check_call(cmd, stderr=dn)
         results = list(SeqIO.parse(outfile, 'fasta'))

    return results

# Author: cmccoy
def which(executable):
    paths = (os.path.join(p, executable)
             for p in os.environ['PATH'].split(':'))
    values = (p for p in paths if os.path.isfile(p) and os.access(p, os.X_OK))
    try:
        return next(values)
    except StopIteration:
        raise OSError("{0} not found in PATH".format(executable))

# Author: cmccoy
def guess_figtree_jar_path():
    """
    Try to guess the path to the figtree jar.

    Should be in $(which figtree)/../lib/figtree.jar

    raises OSError if path cannot be found.
    """
    figtree_path = which('figtree')

    jar_path = os.path.abspath(os.path.join(os.path.dirname(figtree_path),
                                            '..', 'lib', 'figtree.jar'))
    if os.path.exists(jar_path):
        return jar_path
    raise OSError("FigTree jar could not be found.")


figtree_params = """
Begin figtree;
set tipLabels.fontSize=18;
set trees.order=true;
set trees.orderType="decreasing";
End;
"""

# Author: cmccoy
def tree_svg(tree_path, width=1200, height=800, compress=True,
        figtree_path=None):
    """Create svg graphics file for phylogenetic tree file.

    Resulting graphic file is suitable for display in a web page.

    :param tree_path: 
    :param width: width of grpahic in pixels
    :param height: height of graphic in pixels
    :param compress: True to gzip compress the graphic
    				(not sure this works now)
    :param figtree_path: path to figtree executable 
           (defaults to searching relative to 'figtree' executable)
    :returns: SVG graphic string (may be compressed)
    :rtype: string

    """

    
    """
    Generate an SVG from a BEAST nexus file
    """
    name = os.path.basename(os.path.splitext(tree_path)[0]) + '.svg'
    if figtree_path is None:
        figtree_path = guess_figtree_jar_path()
    print(tree_path)
    with open(tree_path, 'r') as infile, tempfile.NamedTemporaryFile() as treefile:
        tree = infile.read()
        treefile.write(tree)
        treefile.write(figtree_params)
        treefile.flush()
	 
        with tempfile.NamedTemporaryFile() as tf,open(os.devnull) as dn:
            # Run headless so graphics operations work without a display.
            # Java memory options taken from the figtree script
            cmd = ['java', '-client', '-Djava.awt.headless=true', '-Xms64m', '-Xmx512m',
               '-jar', figtree_path,
               '-graphic', 'SVG', '-width', width, '-height', height,
               treefile.name, tf.name]
            cmd = map(str, cmd)
            subprocess.check_call(cmd, stdout=dn)
            if compress:
                s = StringIO()
                with gzip.GzipFile(name, 'wb', fileobj=s) as fp:
                    fp.writelines(tf)
                    return s.getvalue()
            else:
                return tf.read()

