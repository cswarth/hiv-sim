#!/usr/bin/python
# make a new BEAST  config file by inserting FASTA sequences into a generic template.
#
# csw 6/3/2014
# The way to use this script:
#
#      mkbeast.py -p foo template.xml sequences.fasta  >beast_in.xml
#
# This will take the fasta sequences and insert them into template.xml
# to produce an XML file that is suitable to pass to BEAST.  The
# prefix "foo" is used to configure the names of the BEAST output
# files.
#
# The fasta sequences should have been generated by a simulation or
# from patients.  Each sequence is expected to have a label like
# "patient1_200_3" which is laid out as
# "<prefix>_<generation>_<index>".  The <generation> is used as a tip
# date when configuring BEAST.
#
# Once the BEAST config file is generated, you would run,
#
#      beast beast_in.xml
#
# This will produce various output files, among which is foo.trees.
# That file gets fed to the 'annotatetrees' program and the output of that gets
# visualized with 'figtree'.
#
from __future__ import print_function

from lxml import etree
import re
from Bio import SeqIO
from collections import defaultdict

import sys, getopt

# from http://en.wikipedia.org/wiki/Autovivification#Python
def autovivify(levels=1, final=dict):
    '''Returns a nested defaultdict with a set number of levels and defined final structure.
    '''
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))
 
'''usage:
   words = autovivify(5, int)
   words["sam"][2012][5][25]["hello"] += 1
   words["sue"][2012][5][24]["today"] += 1
'''

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
    prefix = None
    try:
        opts, args = getopt.getopt(argv,"h")
    except getopt.GetoptError:
        print('mkbeast.py <templatefile> <fastafile>', file=sys.stderr)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('usage: mkbeast.py [-p <prefix>] <templatefile> <fastafile>', file=sys.stderr)
            sys.exit()
        

    templatefile = args[0]
    datafiles = args[1:]

    # datemap stores string dates assigned to {patients, generation} combinations.
    # patients stores xml element corresponding to a patient taxa.
    #
    datemap = defaultdict()
    patients = defaultdict()
    dates = []		# list of date elements
    mstat = []		# list of monophylyStatistic elements
    tmrca = []		# list of tmrcaStatistic elements

    # Parse a generic template and insert sequences from a FASTA file into the middle, separated by the appropriate XML tags.

    tree = etree.parse(templatefile)
    root = tree.getroot()
    
    # Eliminate any existing taxa elements
    # we will replace all these.
    for taxa in tree.xpath("/beast/taxa"):
        taxa.getparent().remove(taxa)

    # Other tags expect there to be two global taxa.
    # we will fill these in with toxons below
    taxa_taxa = etree.Element("taxa", id='taxa')
    taxa_root = etree.Element("taxa", id='root')

    # Eliminate any existing alignment elements
    for alignment in tree.xpath("/beast/alignment"):
        alignment.getparent().remove(alignment)
    alignment = etree.Element("alignment", id='alignment', dataType="nucleotide")

    # eliminate any existing date elements
    # we will replace all these.
    for date in tree.xpath("/beast/date"):
        date.getparent().remove(date)   
    
    # define a regex to extract the generation number from the fasta id string
    # we use this to provide tip dates to BEAST.
    name_regex = re.compile("(?P<patient>^[^_]*)_(?P<generation>.*)_")

    # for each fasta sequence in the data file, create a taxon node and a sequence node.
    for datafile in datafiles:
        with open(datafile, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta") :
                # extract the patient id and generation from the fasta name.
                match = name_regex.search(record.id)
                patient = match.group('patient')
                generation = match.group('generation')
                
                # dates in the beast config correspond to generations in the
                # santa-sim output.
                # we need a mapping from generation to date....
                #
                dateId = "date_"+patient+"_"+generation
#                import pdb; pdb.set_trace()

                if not dateId in datemap:
                    date = datemap[dateId] = generation
                    date = etree.Element("date", id=dateId, value=date, units="days")
                    date.tail = "\n"
                    dates.append(date)

                # create a taxon tag, embed the date within the taxon,
                # append to global taxa element called 'taxa'
                # <taxon id="patient3_1000_3">
                #   <date units="days" direction="forwards" value="1000"/>
                # </taxon>

                taxon = etree.Element("taxon", id=record.id)
                taxon.append(etree.Element("date", idref=dateId))
                taxa_taxa.append(taxon)
                taxa_root.append(etree.Element("taxon", idref=record.id))

                # under <alignment>, create <sequence> with refid to appropriate <taxon>
                # <sequence>
                #   <taxon idref="patient1200"/>
                #   TTTTTTGCAACAGGAGATATAATAGGAAATA
                # </sequence>


                xml = ( '<sequence>'
                        '<taxon idref="{id}"/>'
                        '{sequence}'
                        '</sequence>' ).format(id=record.id, sequence=str(record.seq))
                # sequence = etree.Element("sequence")
                # taxon = etree.Element("taxon", idref=record.id)
                # taxon.tail = str(record.seq) + "\n"
                # sequence.append(taxon)
                # sequence.tail="\n"
                alignment.append(etree.XML(xml))


                # Append a taxon reference to the patient taxa
                # <taxa id="patientFG">
                #   <taxon idref="F02cl15"/>
                # </taxa>
                if patient not in patients:
                    patients[patient] = newPatient(tree, patient)

                p = patients[patient]
                t = etree.Element("taxon", idref=record.id)
                t.tail = "\n"
                p.append(t)


    # place the tmrca and mstat tags after the treeModel tag
    # http://stackoverflow.com/a/7475897/1135316
    treemodel = tree.find(".//treeModel")
    parent = treemodel.getparent()
    for m in reversed(tmrca):
        parent.insert(parent.index(treemodel)+1, m)
    for m in reversed(mstat):
        parent.insert(parent.index(treemodel)+1, m)

    # place the date and patient taxa definitions at the start of the tree.
    pos = 0
    for d in dates:
        root.insert(pos, d)
        pos += 1  # wish python had autoincrement
    root.insert(pos, taxa_taxa)
    pos += 1
    root.insert(pos, taxa_root)
    pos += 1
    for p in patients.values():
        root.insert(pos, p)
        pos += 1
    root.insert(pos, alignment)
    pos += 1
    
    # pretty-print the tree
    indent(root)
    
    print(etree.tostring(tree, pretty_print=True))



# When encountering a sequence from a patient that has not been seen
# before,  define a monopylyetic clade for samples from this patient
# and define some likelihood constraints.
#
def newPatient(tree, patient):
    root = tree.getroot()

    # if we haven't seen this patient before,
    # define a clade for all sequences sampled from this patient
    p = etree.Element("taxa", id=patient)
    p.text = "\n"
    root.append(p)

    # Monitore the monophyly of the specified clade
    # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l17
    #
    # <monophylyStatistic id="monophyly(patientA)">
    # 	<mrca>
    # 		<taxa idref="patientA"/>
    # 	</mrca>
    # 	<treeModel idref="treeModel"/>
    # </monophylyStatistic>
    
    xml = ( '<monophylyStatistic id="monophyly({patient})">'
            '<mrca><taxa idref="{patient}"/></mrca>'
            '<treeModel idref="treeModel"/>'
            '</monophylyStatistic>' ).format(patient=patient)
    treemodel = tree.find(".//treeModel")
    parent = treemodel.getparent()
    parent.insert(parent.index(treemodel)+1, etree.XML(xml))

    # Enforce the monophyly of the specified clade
    #
    # <booleanLikelihood id="booleanLikelihood1">
    # 	<monophylyStatistic idref="monophyly(patientF)"/>
    # </booleanLikelihood>
    prior = tree.find(".//mcmc/posterior/prior")
    xml = ( '<booleanLikelihood id="likelihood({patient})">\n'
            '<monophylyStatistic idref="monophyly({patient})"/>\n'
            '</booleanLikelihood>\n' ).format(patient=patient)
    #prior.append(etree.Comment(etree.tostring(etree.XML(xml))))
    prior.append(etree.XML(xml))


    # Monitor the age of the specified clade
    # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l23
    #
    # <tmrcaStatistic id="tmrca(patientEK)">
    # 	<mrca>
    # 		<taxa idref="patientEK"/>
    # 	</mrca>
    # 	<treeModel idref="treeModel"/>
    # </tmrcaStatistic>
    # http://stackoverflow.com/a/7475897/1135316
    xml = ( '<tmrcaStatistic id="tmrca({patient})">\n'
            '<mrca><taxa idref="{patient}"/></mrca>\n'
            '<treeModel idref="treeModel"/>\n'
            '</tmrcaStatistic>\n' ).format(patient=patient)
    treemodel = tree.find(".//treeModel")
    parent = treemodel.getparent()
    parent.insert(parent.index(treemodel)+1, etree.XML(xml))

    return(p)

if __name__ == "__main__":
   main(sys.argv[1:])
   


