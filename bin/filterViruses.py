import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args=parser.parse_args()

df=pd.read_csv(args.input)
virus_df=df[df['Superkingdom']=='Viruses']
virus_df=virus_df.fillna("NA")

virus_df.to_csv(args.output, index=False)