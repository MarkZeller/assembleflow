#!/usr/bin/env python3 

import sys
import argparse

def find_circular_scaffolds(filein,fileout,kmer,min_len):
    found = 0
    label=''
    seq=''
    for line in filein:
        if line[0] == '>':
            if len(seq) > min_len:
                if seq[0:kmer] == seq[(len(seq)-kmer):len(seq)]:
                    found+=1
                    fileout.write(label+'\n'+seq+'\n')
            label=line.strip()
            seq=''
        else:
            seq+=line.strip()
    return found

if __name__ == '__main__':

    usage = '%prog {options} sequences'

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str, required=True)
    parser.add_argument("--output", "-o", type=str, required=True)
    parser.add_argument("--kmer", "-k", type=str, required=True)
    parser.add_argument("--min_len", "-l", type=str, required=True)
    args = parser.parse_args()

    with open(args.input) as filein, open(args.output, 'a') as fileout:
        num_circular = find_circular_scaffolds(filein,fileout,int(args.kmer),int(args.min_len))
        print("Found " + str(num_circular) + " circular contigs.")

