#!/usr/bin/env python3

import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--magnitudes", "-m", type=str, required=True)
parser.add_argument("--lca", "-l", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)

args = parser.parse_args()


def merge_data_magnitudes():
    file1 = open(args.magnitudes)
    file2 = open(args.lca)
    output =  args.output
    magnitudes = pd.read_csv(file1, sep="\t", header=None)
    magnitudes.columns = ["#queryID", "magnitude"]
    lca = pd.read_csv(file2, sep="\t", header=0)
    combined_df=pd.merge(lca,magnitudes, how="inner", on="#queryID")
    combined_df.to_csv(output, sep=',',index=False)

if __name__ == "__main__":
    merge_data_magnitudes()
