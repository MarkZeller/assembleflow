import pandas as pd
import argparse


parser=argparse.ArgumentParser()
parser.add_argument('sample_ids', nargs='+')
parser.add_argument('files', nargs='+')
parser.add_argument('--output', type=str, required=True)
args=parser.parse_args()

merged_df = pd.DataFrame()

# Loop over each file and merge them into a single DataFrame
for i, file_path in enumerate(args.files):
    sample_name = args.sample_ids[i]
    
    # Load the lca_summary CSV
    df = pd.read_csv(file_path)
    
    # Rename the 'Count' column to reflect the sample
    df = df.rename(columns={'Count': sample_name})
    
    # Merge with the existing DataFrame, matching on taxonomy columns
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on=[
            'Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 
            'Family', 'Genus', 'Species'], how='outer')


# Add a 'Total Count' column by summing across all samples
merged_df['Total'] = merged_df[sample_columns].sum(axis=1)

sample_columns = [col for col in merged_df.columns if col != 'Root' and col != 'Cellular Organisms' and col != 'Superkingdom' and col != 'Phylum' and col != 'Class' and col != 'Superfamily' and col != 'Family' and col != 'Genus' and col != 'Species' and col != 'Total Count']
merged_df[sample_columns]=merged_df[sample_columns].fillna(0).astype(int)
merged_df['Total'] = merged_df['Total'].astype(int)

for col in ['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']:
    merged_df[col] = merged_df[col].fillna("NA")

merged_df = merged_df.sort_values(by=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])

merged_df.to_csv(args.output, index=False)
