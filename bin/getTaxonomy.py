#!/usr/bin/env python3

from itertools import count
import csv
import argparse
from ete3 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)
args = parser.parse_args()

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, 'NA') for rank in desired_ranks}

if __name__ == '__main__':
    taxids = []
    desired_ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = list()
    with open(args.input) as diamondhits:
        next(diamondhits) 
        for hit in diamondhits:
            queryID = hit.split(",")[0].strip()
            taxid = hit.split(",")[1].strip()
            if taxid:
                taxid = int(float(taxid))
            else:
                taxid = 1
            similarity = hit.split(",")[2].strip()
            magnitude = hit.split(",")[3].strip()
            results.append(list())
            results[-1].append(str(queryID))
            results[-1].append(str(taxid))
            results[-1].append(str(similarity))
            results[-1].append(str(magnitude))
            ranks = get_desired_ranks(taxid, desired_ranks)
            for key, rank in ranks.items():
                if rank != 'NA':
                    results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
                else:
                    results[-1].append(rank)

    with open (args.output, 'w') as out_file:
        csv = []
        header = ['queryID','taxid','similarity','magnitude']
        header.extend(desired_ranks)
        csv.append(header)
        print(','.join(header), file=out_file)
    
    
        for result in results:
            csv.append(result)
            print(','.join(result), file=out_file)
