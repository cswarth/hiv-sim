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
    startDate = datetime.strptime("6/13/2014", "%m/%d/%Y")
    prefix = None
    try:
        opts, args = getopt.getopt(argv,"hd:")
    except getopt.GetoptError:
        print('mkbeast.py <templatefile> <fastafile>', file=sys.stderr)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('usage: mkbeast.py [-p <prefix>] <templatefile> <fastafile> [<fastafile> ...]', file=sys.stderr)
            sys.exit()
        elif opt == '-d':
            startDate = date.strptime(arg, "%m/%d/%Y")

    templatefile = args[0]
    datafiles = args[1:]

    # patients dict stores instance of class Patient keyed by patient
    # id string.  we create a new one when we run into a patient id
    # string that we haven't seem before.  patient id string is
    # derived from the first part of the name on an alignment,
    # e.g. 'patient1_100_1' yields patient id 'patient1'
    #
    patients = defaultdict()

    # Parse a generic template and insert sequences from a FASTA file into the middle,
    # separated by the appropriate XML tags.

    tree = etree.parse(templatefile)
    
    # Eliminate any existing taxa elements
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
    
    # insert some tags right at the top of the XML tree.   The tags will be populated later.
    root = tree.getroot()
    root.insert(2, etree.Element("alignment", id='alignment', dataType="nucleotide"))

    # define a regex to extract the generation number from the fasta id string
    # we use this to provide tip dates to BEAST.
    name_regex = re.compile("(?P<patient>^[^_]*)_(?P<generation>.*)_")

    # for each fasta sequence in the data file, create a taxon node and a sequence node.
    for datafile in datafiles:
        with open(datafile, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta") :
                # split the tag into fields separated by whitespace
                fields = record.id.split('|')
                # extract the patient id and generation from the fasta name.
                match = name_regex.search(fields[0])
                patient = match.group('patient')
                generation = match.group('generation')

                # Create a patient instance if we haven't seen this patient id before
                if patient not in patients:
                    patients[patient] = Patient(tree, patient)
                p = patients[patient]

                # >patient3_100_2|patient2.fa|400|600
                #  <sequence identifier>|<parent file>|<parent generation>|<cumulative generation>
                # the date of this sequence is calculated as,
                # 	 cumulative generation of parent
                #  + sampled generation from current patient.
                #

                sampleDate = startDate
                sampleDate += timedelta(int(fields[3]))
                sampleDate += timedelta(int(generation))
                p.addSequence(generation, str(record.seq), sampleDate)


    # pretty-print the tree
    indent(root)
    
    print(etree.tostring(tree, pretty_print=True))




class BeastTaxa(object):		# new style class inherits from object
    							# which means subclasses can use super()
                                # https://docs.python.org/2/library/functions.html#super
    '''
        Define a class to hold data and methods associated with taxa definitions in
        BEAST config files.

        We usually define two or three kinds of taxa in BEAST config files:

        	global taxa like 'taxa' and 'root',
            per-patient taxa like 'patient1', 
            per-patient,per-generation taxa like 'patient1_200'

        This class is a general BEAST taxa class.  Given a taxa ID it knows
        how to define a taxa, and how to add a taxon reference to the taxa.
        A coupld of helper functions to defining and enforcing various constraints
        on taxa are also packages in this class. 
        This class does NOT know how to define a single taxon - that is done elsewhere.
    '''

    def __init__(self, tree, id, position=None):
        self._tree = tree
        self._id = id

        # insert taxa after all existing taxa definitions, or at start of tree if no taxa exist.
        if position is None:
            position = 0
            e = tree.find("./taxa[last()]")
            if e is not None:
                position = e.getparent().index(e)
                position += 1
        taxa = etree.Element("taxa", id=self._id)
        self._tree.getroot().insert(position, taxa)

    def addDate(self, dateId, sampleDate, position=None):
        # dates in the beast config correspond to generations in the
        # santa-sim output.

        if position is None:
            position = 0
            e = self._tree.find("./date[last()]")
            if e is not None:
                position = e.getparent().index(e)
                position += 1
        date = etree.Element("date", id=dateId, value=sampleDate.strftime('%d/%m/%Y'), units="days")
        self._tree.getroot().insert(position, date)


    def addTaxon(self, taxonId, dateId=None):
        # create a taxon tag, embed the date within the taxon,
        if dateId is None:
            xml = ( '<taxon id="{}"/>' ).format(taxonId)
        else:
            xml = ( '<taxon id="{}">'
                    '  <date idref="{}"/>'
                    '</taxon>' ).format(taxonId, dateId)
        xpath = "./taxa[@id='{}']".format(self._id)
        self._tree.find(xpath).append(etree.XML(xml))

    def addTaxonRef(self, taxonId):
        '''
        Append a taxon reference to the per-patient-generation taxa
        '''
        xpath = "./taxa[@id='{}']".format(self._id)
        self._tree.find(xpath).append(etree.Element("taxon", idref=taxonId))

    def monophyly(self):
        # Monitor the monophyly of the specified taxa
        # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l17
        # http://stackoverflow.com/a/7475897/1135316
        xml = ( '<monophylyStatistic id="monophyly({id})">'
                '<mrca><taxa idref="{id}"/></mrca>'
                '<treeModel idref="treeModel"/>'
                '</monophylyStatistic>' ).format(id=self._id)
        treemodel = self._tree.find(".//treeModel")
        parent = treemodel.getparent()
        parent.insert(parent.index(treemodel)+1, etree.XML(xml))

        # Enforce the monophyly of the specified patient taxa
        prior = self._tree.find(".//mcmc/posterior/prior")
        xml = ( '<booleanLikelihood id="likelihood({id})">\n'
                '<monophylyStatistic idref="monophyly({id})"/>\n'
                '</booleanLikelihood>\n' ).format(id=self._id)
        #prior.append(etree.Comment(etree.tostring(etree.XML(xml))))
        prior.append(etree.XML(xml))

        # Constrain the starting tree to make this taxa monophyletic.
        xml = ( '<coalescentSimulator>'
                '  <taxa idref="{id}"/>'
                '  <constantSize idref="constant"/>'
                '</coalescentSimulator>').format(id=self._id)
        cs = self._tree.find("./coalescentSimulator[@id]")
        cs.append(etree.XML(xml))

    def tmrca(self):
        # Monitor the age of the specified clade
        # http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l23
        # http://stackoverflow.com/a/7475897/1135316
        xml = ( '<tmrcaStatistic id="tmrca({patient})">\n'
                '<mrca><taxa idref="{patient}"/></mrca>\n'
                '<treeModel idref="treeModel"/>\n'
                '</tmrcaStatistic>\n' ).format(patient=self._id)
        treemodel = self._tree.find(".//treeModel")
        parent = treemodel.getparent()
        parent.insert(parent.index(treemodel)+1, etree.XML(xml))


    def ancestralTrait(self):
        '''
        Log the ancestral sequences for this taxa.
        '''
        # Log the inferred ancestral sequence of this patient.
        # see http://bodegaphylo.wikispot.org/3._Editing_XML_Input_File#l17
        # http://stackoverflow.com/a/7475897/1135316
        xml = ( '<ancestralTrait name="{patient}"  traitName="states">'
                '<treeModel idref="treeModel"/>'
				'<ancestralTreeLikelihood idref="treeLikelihood"/>'
                '</ancestralTrait>' ).format(patient=self._id)
        node = self._tree.find(".//ancestralTrait[@traitName='states']")
        node.getparent().append(etree.XML(xml))
        


class Generation(BeastTaxa):
    _classInit = False
    _taxaTaxa = None
    _rootTaxa = None

    def __init__(self, tree, generationId, sampleDate):
        super(Generation, self).__init__(tree, generationId)
        self._sindex = 1		# unique index per taxon 

        if not Generation._classInit:
            Generation._taxaTaxa = BeastTaxa(self._tree, 'taxa', 0)
            Generation._rootTaxa = BeastTaxa(self._tree, 'root', 1)
            Generation._classInit = True

        # each generation per patient has a unique date.
        dateId = "date("+self._id+")"
        self.addDate(dateId, sampleDate)

        self.tmrca()
        self.monophyly()   # declatre each generation to be monophyletic

        
    def addTaxon(self):
        '''
        the tag has format
        >patient2_100_1 patient1 200 6/13/1995
        '''
        dateId = "date("+self._id+")"
        taxonId = self._id+"_"+str(self._sindex)
        self._sindex += 1

        Generation._taxaTaxa.addTaxon(taxonId, dateId)
        Generation._rootTaxa.addTaxonRef(taxonId)

        super(Generation,self).addTaxonRef(taxonId)
        return(taxonId)


# define a patient class.  I would have preferred to make this a
# subclass of etree._Element, but I couldn't get that to work.
class Patient(BeastTaxa):
    def __init__(self, tree, patientid):
        super(Patient, self).__init__(tree, patientid)

        self._generations = defaultdict()
        
        self.ancestralTrait()


    def addSequence(self, generation, sequence, sampleDate):
        '''
        the tag has format
        >patient2_100_1 patient1 200 6/13/1995
        '''
        generationId = self._id+"_"+generation
        # Create a generation instance if we haven't seen this generation id before
        if generationId not in self._generations:
            self._generations[generationId] = Generation(self._tree, generationId, sampleDate)
        g = self._generations[generationId]

        taxonId = g.addTaxon()
        self.addTaxonRef(taxonId)
        

        # Create a <sequence> with refid to appropriate <taxon>
        #
        xml = ( '<sequence>'
                '  <taxon idref="{id}"/>'
                '  {sequence}'
                '</sequence>' ).format(id=taxonId, sequence=str(sequence))
        self._tree.find("./alignment").append(etree.XML(xml))




    
if __name__ == "__main__":
   main(sys.argv[1:])
   


