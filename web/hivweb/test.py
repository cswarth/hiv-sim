import matplotlib

matplotlib.use('GTKAgg')

from Bio import Phylo
import pylab

tree = Phylo.read('apaf.xml', 'phyloxml')

Phylo.draw_graphviz(tree)
pylab.show()
