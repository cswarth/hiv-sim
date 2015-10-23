#!/usr/bin/env bash
set -e # exit on error


export PATH=~matsengrp/local/bin/:~/src/matsen/hiv-sim/bin:$PATH
export TEMPLATES=~/src/matsen/hiv-sim/templates



function usage() {
    echo "infer.sh <sequences.fasta>";
}

function echotty() {
    # echo to stderr if bound to a tty, otherwise stay silent.
    if [  "$verbose" = true ] && [ -t 2 ] ; then
	echo "$@" 1>&2
    fi
}

function echocmd() {
    if [ "$verbose" = true ]; then
	echo "$@" 1>&2
    fi
    eval $@
}

# get full directory name of the script
# http://stackoverflow.com/a/246128/1135316
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )"  >/dev/null && pwd )"

PATH=${DIR}:${DIR}/../../bin:${PATH}
files=()

while [[ $# > 0 ]]
do
key="$1"
echocmd=false
case $key in
    -v|--verbose)
	verbose=true
    ;;
    -*)
	echo 'Unknown argument - "$key"'
	usage
	exit 1
	;;
    *)
	files+=($key)
    ;;
esac
shift # past argument
done


if [[ ${#files[@]} == 0 ]]; then
    usage
    exit 1
fi

# see http://unix.stackexchange.com/a/84980/100709
# for a mktmpdir that works on both linux and OS X
outdir=`mktemp -d 2>/dev/null || mktemp -d -t 'tmpdir'`
outdir="output"
echotty "outdir= ${outdir}"
for sample in "${files[@]}"
do
    label=$(basename $sample)
    label=${label%.*}

    echotty "Creating alignment with PRANK..." 
    echocmd "prank -d=${sample} -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"

    echotty "Testing for multiple founders..."
    moi=$(keele ${sample})
    if [[ "${moi}" == *"multiple infection"* ]]
    then
	# multiple founders
	# split just below the root and sequences separately.
	echocmd "treesplit.py -o  ${outdir} ${outdir}/prank.best.dnd ${sample}"
	echocmd "prank -d=${outdir}/left.fasta -o=${outdir}/left -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"
	echocmd "prank -d=${outdir}/right.fasta -o=${outdir}/right -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"

	founder1=$(prankroot.py ${outdir}/left.best.anc.dnd ${outdir}/left.best.anc.fas)
	founder2=$(prankroot.py ${outdir}/right.best.anc.dnd ${outdir}/right.best.anc.fas)
	# if [[ "${founder1}" == "${founder2}" ]]
	# then
	#     echo "Founders are identical"
	# else
	#     echo "Founders are different"
	# fi
	echo ">" ${label}_founder_1
	echo ${founder1}
	echo ">" ${label}_founder_2
	echo ${founder2}

	echocmd "prankroot.py ${outdir}/right.best.anc.dnd ${outdir}/right.best.anc.fas"
	
    else
	# single founder - process all sequences together
	# extract sequences at root of guide tree and output in FASTA format.
	echo ">" ${label}_founder
	echocmd "prankroot.py ${outdir}/prank.best.anc.dnd ${outdir}/prank.best.anc.fas"
    fi
    
	
    # inserts sample sequences into a beast configurarion file.
    # parses date fromt he sequence ids and converts to tip date for BEAST
    echotty "Make BEAST config file..."
    echocmd "mkbeast_rv217.py --template ${TEMPLATES}/beast_strict.template ${outdir}/prank.best.fas > ${outdir}/beast_in.xml"

    # sigh...if using the -working option of beast, the configuration file must be specified with an absoute path.
    echotty "Running BEAST..."
    echocmd "beast -working -overwrite -beagle ${PWD}/${outdir}/beast_in.xml >${outdir}/beastcmd.log  2>&1"

    echotty 'Extracting estimated time of infection'
    echocmd 'posterior_toi.py ${outdir}/beastout.log $sample'

    # echotty "Exracting guide tree..." 
    # echocmd "treeannotator ${outdir}/beastout.trees > ${outdir}/mcc.tree 2>${outdir}/treeannotator.log"

    # tr "/" "-" < ${outdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${outdir}/guidetree.tree

    # echotty "Inferring ancestral states with PRANK..." 
    # echocmd "prank -d=${sample} -t=${outdir}/guidetree.tree -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${outdir}/prankcmd.log 2>&1"

done

