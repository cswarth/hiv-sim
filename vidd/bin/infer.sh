#!/usr/bin/env bash
set -e # exit on error
function usage() {
    echo "infer.sh <sequences.fasta>";
}

function echotty() {
    # echo to stderr if bound to a tty, otherwise stay silent.
    if [ -t 2 ] ; then
	echo "$@" 1>&2
    fi
}

PATH=../bin:${PATH}
if [[ $# > 0 ]]; then
   sample="$1";
else
    usage
    exit 1
fi

# see http://unix.stackexchange.com/a/84980/100709
# for a mktmpdir that works on both linux and OS X
mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

# inserts sample sequences into a beast configurarion file.
# parses date fromt he sequence ids and converts to tip date for BEAST
echotty "Make BEAST config file..."
mkbeast_rv217.py --template ../templates/relaxed-clock.xml ${sample} > ${mytmpdir}/beast_in.xml

echotty "Running BEAST..."
beast -working -overwrite -beagle ${mytmpdir}/beast_in.xml >${mytmpdir}/beastcmd.log  2>&1 

echotty "Exracting guide tree..." 
treeannotator ${mytmpdir}/beastout.trees > ${mytmpdir}/mcc.tree 2>${mytmpdir}/treeannotator.log

tr "/" "-" < ${mytmpdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${mytmpdir}/guidetree.tree

echotty "Inferring ancestral states with PRANK..." 
prank -d=${sample} -t=${mytmpdir}/guidetree.tree -o=${mytmpdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${mytmpdir}/prankcmd.log 2>&1 

# extract sequences at root of guide tree and output in FASTA format.
echo ">founder"
./prankroot.py ${mytmpdir}/prank.best.anc.dnd ${mytmpdir}/prank.best.anc.fas 

