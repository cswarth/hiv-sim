# Build a multiple alignment file using sumulated evolved sequences.
#
# SANTA creates the sequences according to the parameters in the santa_test.xml file.
# The output is massages into fasta format that is then fed to clustalw for multiple alignment.
#
# The next trick is the get BEAST to infer a phylogenic tree from these sequences. 
# 

TOOLS=~/src/matsen/tools
PREFIX=hiv-pol

all: patient2_best.tree

combined_best.tree: combined.trees
	treeannotator $^ >$@

combined.trees: combined_beast.xml
	beast -overwrite $^

combined_beast.xml : patient1.fa patient2.fa
	./mkbeast.py -p combined patient_beast_template.xml $^ >$@

patient2_santa.xml : patient1_samples.fa
	./mksanta.py -p patient2 patient_santa_template.xml $^ >$@

patient2.fa: patient2_santa.xml
	java -jar ${TOOLS}/santa-sim/dist/santa.jar $^

patient2_beast.xml : patient2.fa
	./mkbeast.py -p patient2 patient_beast_template.xml $^ >$@

patient2.trees : patient2_beast.xml
	beast -overwrite $^

patient2_best.tree : patient2.trees
	treeannotator $^ >$@


patient1_samples.fa: patient1_santa.xml
	# create a series of mutated sequences.
	java -jar ${TOOLS}/santa-sim/dist/santa.jar $^


${PREFIX}.aln : ${PREFIX}.fa
	# mutiple alignments!
	muscle <${PREFIX}.fa >${PREFIX}.aligned.fa



clean:   # clean intermediate files
	rm -f patient2.*
	rm -f patient2_*
	rm -f *.ops
	rm -f *.log
	rm -f *.trees
	rm -f *.nex
	rm -f *.aln


clobber:  clean
	rm -f ${PREFIX}.aln ${PREFIX}.dnd
