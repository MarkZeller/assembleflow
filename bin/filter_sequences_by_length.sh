#!/usr/bin/env bash

usage() { echo "Usage: filter_sequences_by_length.sh [-i input fasta file] [-s minimum sequence length] [-o output fasta file]" 1>&2; exit 1; }

while getopts ":i:s:o:" arg; do
    case "$arg" in
        i)
            i=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${i}" ] || [ -z "${s}" ] || [ -z "${o}" ]; then
    usage
fi

awk -v n=${s} '/^>/{ if(l>n) print b; b=$0;l=0;next }{l+=length;b=b ORS $0}END{if(l>n) print b }' ${i} > ${o}
