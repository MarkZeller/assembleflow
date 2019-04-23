#!/usr/bin/env python3

import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
args = parser.parse_args()

def summarize_taxonomy():
    file = open(args.input)
    ranks = ['superkingdom','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    df = pd.read_csv(file, sep="\t", header=0)
    df['magnitude'] = df['magnitude'].apply(pd.to_numeric)
    
    #summarize viruses only
    virus = df[df.superkingdom == "Viruses"]
    for rank in ranks:
        level = virus.groupby(rank)['magnitude'].agg((len, sum))
        level.to_csv('virus_'+rank+".tab", sep='\t')

    #summarize all superkingdoms
    for rank in ranks:
        level = df.groupby(rank)['magnitude'].agg((len, sum))
        level.to_csv(rank+".tab", sep='\t')


if __name__ == "__main__":
    summarize_taxonomy()

