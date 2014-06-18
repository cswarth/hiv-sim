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
from datetime import date
from datetime import timedelta

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

    # Eliminate any existing alignment elements
    for alignment in tree.xpath("/beast/alignment"):
        alignment.getparent().remove(alignment)

    # eliminate any existing date elements
    for date in tree.xpath("/beast/date"):
        date.getparent().remove(date)   
    
    # eliminate existing coalescent taxa tags
    for cs in tree.xpath("/beast/coalescentSimulator[@id]/coalescentSimulator"):
        cs.getparent().remove(cs)   
    
    # insert some tags that we will populate later
    root.insert(0, etree.Element("taxa", id='taxa'))
    root.insert(1, etree.Element("taxa", id='root'))
    root.insert(2, etree.Element("alignment", id='alignment', dataType="nucleotide"))

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

                # Append a taxon reference to the patient taxa
                # <taxa id="patientFG">
                #   <taxon idref="F02cl15"/>
                # </taxa>
                if patient not in patients:
                    patients[patient] = newPatient(tree, patient)

                p = patients[patient]
                p.addSequence(generation, str(record.seq))


    # place the tmrca and mstat tags after the treeModel tag
    # http://stackoverflow.com/a/7475897/1135316
    treemodel = tree.find(".//treeModel")
    parent = treemodel.getparent()
    for m in reversed(tmrca):
        parent.insert(parent.index(treemodel)+1, m)
    for m in reversed(mstat):
        parent.insert(parent.index(treemodel)+1, m)

    # pretty-print the tree
    indent(root)
    
    print(etree.tostring(tree, pretty_print=True))


# When encountering a sequence from a patient that has not been seen
# before,  define a monopylyetic clade for samples from this patient
# and define some likelihood constraints.
#
def newPatient(tree, patient):
    p = Patient(tree, patient)
    return(p)
    
# define a patient class.  I would have preferred to make this a
# subclass of etree._Element, but I couldn't get that to work.
class Patient():

    def __init__(self, tree, patientid):
        self._id = patientid
        self._tree = tree
        self._sindex = 1
        self._sampledate = date(2014,6, 13)

        self.monophyly()
        self.tmrca()
        patientTaxa = etree.Element("taxa", id=patientid)

        # insert patient taxa right before <alignment>
        a = tree.find("./alignment[@id='alignment']")
        self._tree.getroot().insert(a.getparent().index(a), patientTaxa)


    def monophyly(self):
        # Monitor the monophyly of the specified clade
        # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l17
        # http://stackoverflow.com/a/7475897/1135316
        xml = ( '<monophylyStatistic id="monophyly({patient})">'
                '<mrca><taxa idref="{patient}"/></mrca>'
                '<treeModel idref="treeModel"/>'
                '</monophylyStatistic>' ).format(patient=self._id)
        treemodel = self._tree.find(".//treeModel")
        parent = treemodel.getparent()
        parent.insert(parent.index(treemodel)+1, etree.XML(xml))

        # Enforce the monophyly of the specified patient taxa
        prior = self._tree.find(".//mcmc/posterior/prior")
        xml = ( '<booleanLikelihood id="likelihood({patient})">\n'
                '<monophylyStatistic idref="monophyly({patient})"/>\n'
                '</booleanLikelihood>\n' ).format(patient=self._id)
        #prior.append(etree.Comment(etree.tostring(etree.XML(xml))))
        prior.append(etree.XML(xml))

        # Constrain the starting tree to make each patient monophyletic.
        xml = ( '<coalescentSimulator>'
                '  <taxa idref="{id}"/>'
                '  <constantSize idref="constant"/>'
                '</coalescentSimulator>').format(id=self._id)
        cs = self._tree.find("./coalescentSimulator[@id]")
        cs.append(etree.XML(xml))

    def tmrca(self):
        # Monitor the age of the specified clade
        # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l23
        # http://stackoverflow.com/a/7475897/1135316
        xml = ( '<tmrcaStatistic id="tmrca({patient})">\n'
                '<mrca><taxa idref="{patient}"/></mrca>\n'
                '<treeModel idref="treeModel"/>\n'
                '</tmrcaStatistic>\n' ).format(patient=self._id)
        treemodel = self._tree.find(".//treeModel")
        parent = treemodel.getparent()
        parent.insert(parent.index(treemodel)+1, etree.XML(xml))

    def addSequence(self, generation, sequence):

        sequenceId = self._id+"_"+generation+"_"+str(self._sindex)
        dateId = "date_"+self._id+"_"+generation
        self._sindex += 1

        # add a sample date if we don't already have one.
        sampleDate = "./date[@id='{}']".format(dateId)
        if self._tree.find(sampleDate) is None:
            # dates in the beast config correspond to generations in the
            # santa-sim output.
            # we need a mapping from generation to date....
            #
            date = self._sampledate + timedelta(int(generation))
            date = etree.Element("date", id=dateId, value=date.strftime('%d/%m/%Y'), units="days")
            self._tree.getroot().insert(0, date)


        # create a taxon tag, embed the date within the taxon,
        # append to global taxa element called 'taxa'
        xml = ( '<taxon id="{}">'
                '  <date idref="{}"/>'
                '</taxon>' ).format(sequenceId, dateId)
        self._tree.find("./taxa[@id='taxa']").append(etree.XML(xml))
        self._tree.find("./taxa[@id='root']").append(etree.Element("taxon", idref=sequenceId))

        # append a taxon to the per-patient taxa
        patientTaxa = "./taxa[@id='{}']".format(self._id)
        self._tree.find(patientTaxa).append(etree.Element("taxon", idref=sequenceId))

        # Create a <sequence> with refid to appropriate <taxon>
        #
        xml = ( '<sequence>'
                '  <taxon idref="{id}"/>'
                '  {sequence}'
                '</sequence>' ).format(id=sequenceId, sequence=str(sequence))
        self._tree.find("./alignment").append(etree.XML(xml))




    
if __name__ == "__main__":
   main(sys.argv[1:])
   


