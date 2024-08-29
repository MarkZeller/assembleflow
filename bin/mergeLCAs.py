#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import glob

def process_file(file_path):
    sample_name=os.path.basename(file_path).replace("_lca_summary.csv", "")
    df=pd.read_csv(file_path)
    tax_cols=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']
    df=df[tax_cols+["Count"]]
    df=df.rename(columns={"Count": sample_name})
    return df

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output", type=str, required=True)
    args=parser.parse_args()

    input_files=glob.glob("*_lca_summary.csv")
    if not input_files:
        print("No lca_summary files found in the current directory")
        return

    merged_df=pd.DataFrame()
    for file_path in input_files:
        df=process_file(file_path)

        if merged_df.empty:
            merged_df=df
        else:
            merged_df=pd.merge(merged_df, df, on=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'], how='outer')

    tax_cols=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']
    count_cols=[col for col in merged_df.columns if col not in tax_cols]
    merged_df["Total Count"]=merged_df[count_cols].sum(axis=1)
    merged_df[count_cols]=merged_df[count_cols].fillna(0).astype(int)
    merged_df["Total Count"]=merged_df["Total Count"].astype(int)

    for col in tax_cols:
        merged_df[col] = merged_df[col].fillna("NA")

    merged_df=merged_df.sort_values(by=tax_cols)

    final_cols=tax_cols+count_cols+["Total Count"]
    merged_df=merged_df[final_cols]

    merged_df.to_csv(args.output, index=False)

if __name__=="__main__":
    main()
