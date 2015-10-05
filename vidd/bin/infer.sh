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

function echocmd() {
    if [ "$verbose" = true ]; then
	echo "$@" 1>&2
    fi
    eval $@
}

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
mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

for sample in "${files[@]}"
do
    
    # inserts sample sequences into a beast configurarion file.
    # parses date fromt he sequence ids and converts to tip date for BEAST
    echotty "Make BEAST config file..."
    echocmd "mkbeast_rv217.py --template ../templates/relaxed-clock.xml ${sample} > ${mytmpdir}/beast_in.xml"

    echotty "Running BEAST..."
    echocmd "beast -working -overwrite -beagle ${mytmpdir}/beast_in.xml >${mytmpdir}/beastcmd.log  2>&1"

    echotty "Exracting guide tree..." 
    echocmd "treeannotator ${mytmpdir}/beastout.trees > ${mytmpdir}/mcc.tree 2>${mytmpdir}/treeannotator.log"

    tr "/" "-" < ${mytmpdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${mytmpdir}/guidetree.tree

    echotty "Inferring ancestral states with PRANK..." 
    echocmd "prank -d=${sample} -t=${mytmpdir}/guidetree.tree -o=${mytmpdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${mytmpdir}/prankcmd.log 2>&1"

    # extract sequences at root of guide tree and output in FASTA format.
    echo ">" $(basename $sample)
    echocmd "prankroot.py ${mytmpdir}/prank.best.anc.dnd ${mytmpdir}/prank.best.anc.fas"
done

