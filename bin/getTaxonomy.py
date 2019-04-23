#!/usr/bin/env python3

from itertools import count
import csv
import argparse
from ete3 import NCBITaxa

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
args = parser.parse_args()

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

if __name__ == '__main__':
    taxids = []
    desired_ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = list()
    with open(args.input) as diamondhits:
        for hit in diamondhits:
            queryID = hit.split("\t")[0].strip()
            taxid = hit.split("\t")[1].strip()
            similarity = hit.split("\t")[2].strip()
            magnitude = hit.split("\t")[3].strip()
            results.append(list())
            results[-1].append(str(queryID))
            results[-1].append(str(taxid))
            results[-1].append(str(similarity))
            results[-1].append(str(magnitude))
            try:
                ranks = get_desired_ranks(taxid, desired_ranks)
                for key, rank in ranks.items():
                    if rank != '<not present>':
                        results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
                    else:
                        results[-1].append(rank)
          
            except:
                ranks = ['NA','NA','NA','NA','NA','NA','NA','NA']
                for rank in ranks:
                    results[-1].append(rank)
    
    for result in results:
        for index, obj in zip(count(len(result) - 1, -1), reversed(result)):
            if index > 2:
                if obj != "NA":
                    if result[index - 1] == '<not present>':
                        result[index - 1] = result[index]

    #generate the header
    header = ['queryID','taxid','similarity','magnitude']
    header.extend(desired_ranks)
    print('\t'.join(header))

    #print the results
    for result in results:
        print('\t'.join(result))

