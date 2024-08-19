#!/usr/bin/env python3

import pandas as pd
from ete3 import NCBITaxa
import argparse

ncbi=NCBITaxa()

#function which gets lowest common ancestor for a list of target accession IDs
def get_lca(tax_ids):
	try:
		tax_ids=[int(float(tax_id)) for tax_id in tax_ids if pd.notna(tax_id)]

		lca=ncbi.get_topology(tax_ids, intermediate_nodes=False)
		
		if lca:
			return lca.taxid
		else:
			return None 
	except Exception as e:
		print(f"Error finding LCA: {e}")
		return None


def translate_taxid(taxid):
	try:
		
		taxid=str(int(float(taxid)))
		lineage=ncbi.get_lineage(taxid)
		lineage_names=ncbi.get_taxid_translator(lineage)
		translated_names=[lineage_names.get(taxid, "Unknown") for taxid in lineage]
		ranks=ncbi.get_rank(lineage)
		rank_names=[ranks.get(taxid, "Unknown") for taxid in lineage]

		root=cellular_organisms=superkingdom=phylum=class_=family=superfamily=genus=species="NA"

		rank_dict=dict(zip(lineage, rank_names))
		name_dict=dict(zip(lineage, translated_names))

		for tax in lineage:
			rank=rank_dict.get(tax, "NA")
			if tax==1:
				root=name_dict.get(tax, "root")
			elif tax==131567:
				cellular_organisms=name_dict.get(tax, "cellular organisms")
			elif rank=="superkingdom":
				superkingdom=name_dict.get(tax, "NA")
			elif rank=="phylum":
				phylum=name_dict.get(tax, "NA")
			elif rank=="class":
				class_=name_dict.get(tax, "NA")
			elif rank=="superfamily":
				superfamily=name_dict.get(tax, "NA")
			elif rank=="family":
				family=name_dict.get(tax, "NA")
			elif rank=="genus":
				genus=name_dict.get(tax, "NA")
			elif rank=="species":
				species=name_dict.get(tax, "NA")

		return root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species
	except Exception as e:
		print(f"Error translating taxid {taxid}: {e}")
		return "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

def main(diamond_out, magnitudes, output_table, output_lca_summary, output_hits_summary, output_virus_summary):
	df=pd.read_csv(diamond_out, sep='\t', header=None)
	magnitudes_df=pd.read_csv(magnitudes, sep='\t', header=None, names=['Query', 'Coverage'])

	query_hits={}

	for _, row in df.iterrows():
		query=row[0]
		hits=row[1]
		hit_taxid=row[2]
		percent_id=row[3]
		length=row[4]
		mismatch=row[5]
		evalue=row[7]
		bitscore=row[8]

		if query not in query_hits:
			query_hits[query]={"hits": [], "tax_ids": [], "percent_ids": [], "lengths": [], "mismatches": [], "evalues": [], "bitscores": []}
		if pd.notna(hit_taxid):
			try:
				tax_id=str(int(float(hit_taxid)))
				query_hits[query]["hits"].append(hits)
				query_hits[query]["tax_ids"].append(tax_id)
				query_hits[query]["percent_ids"].append(percent_id)
				query_hits[query]["lengths"].append(length)
				query_hits[query]["mismatches"].append(mismatch)
				query_hits[query]["evalues"].append(evalue)
				query_hits[query]["bitscores"].append(bitscore)
			except ValueError:
				print(f"Invalid taxonomy ID: {hit_taxid}")


	output_summary=[]
	lca_summary=[]
	virus_summary=[]

	for query, data in query_hits.items():
		hits=data["hits"]
		tax_ids=data["tax_ids"]
		percent_ids=data["percent_ids"]
		lengths=data["lengths"]
		mismatches=data["mismatches"]
		evalues=data["evalues"]
		bitscores=data["bitscores"]

		lca = get_lca(tax_ids)

		if lca:
			try:

				lca_name=ncbi.get_taxid_translator([lca]).get(lca)
				lca_name_id=f"{lca_name} ({lca})"
				root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species=translate_taxid(lca)
				taxonomy_info=', '.join([root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species])

				lca_info={
					'Root': root,
					'Cellular Organisms': cellular_organisms,
					'Superkingdom': superkingdom,
					'Phylum': phylum,
					'Class': class_,
					'Superfamily': superfamily,
					'Family': family,
					'Genus': genus,
					'Species': species,
					'Count': 1
				}
				lca_summary.append(lca_info)

				if superkingdom=="Viruses":
					virus_summary.append(lca_info)

				avg_percents=round(sum(percent_ids)/len(percent_ids), 2) if percent_ids else 0
				avg_lengths=round(sum(lengths)/len(lengths), 2) if lengths else 0
				avg_mismatches=round(sum(mismatches)/len(mismatches), 2) if mismatches else 0
				avg_evalues=round(sum(evalues)/len(evalues), 2) if evalues else 0
				avg_bitscores=round(sum(bitscores)/len(bitscores), 2) if bitscores else 0


				output_summary.append({
					'Query': query,
					'Hits': ', '.join(hits),
					'Hit TaxIDs': ', '.join(tax_ids),
					'LCA': lca_name_id,
					'LCA Tax Info': taxonomy_info,
					'Av % Identity': avg_percents,
					'Av Align Lengths': avg_lengths,
					'Av Mismatches': avg_mismatches,
					'Av e-value': avg_evalues,
					'Av Bit Score': avg_bitscores
				})
				
				
			except Exception as e:
				print(f"Error processing LCA {lca} for query {query}: {e}")
		else:
			print(f"No LCA found for query {query}")
		
	output_df=pd.DataFrame(output_summary)

	output_df['Query']=output_df['Query'].astype(str)
	magnitudes_df['Query']=magnitudes_df['Query'].astype(str)
	output_df['Query']=output_df['Query'].str.strip().str.upper()
	magnitudes_df['Query']=magnitudes_df['Query'].str.strip().str.upper()

	merged_output_df=pd.merge(output_df, magnitudes_df, on='Query', how='inner')
	merged_output_df.to_csv(output_table, sep=',', index=False)

	lca_summary_df=pd.DataFrame(lca_summary)
	lca_summary_df=lca_summary_df.groupby(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']).sum().reset_index()
	lca_summary_df=lca_summary_df.sort_values(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])
	lca_summary_df.to_csv(output_lca_summary, sep=',', index=False)

	if virus_summary:
		virus_summary_df=pd.DataFrame(virus_summary)
		virus_summary_df=virus_summary_df.groupby(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']).sum().reset_index()
		virus_summary_df=virus_summary_df.sort_values(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])
		virus_summary_df.to_csv(output_virus_summary, sep=',', index=False)
	else:
		print("No virus LCAs found")

if __name__ == "__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', '--input', type=str, required=True)
	parser.add_argument('-m', '--magnitudes', type=str, required=True)
	parser.add_argument('-o', '--output', type=str, required=True)
	parser.add_argument('-l', '--lca', type=str, required=True)
	parser.add_argument('-v', '--virus', type=str, required=True)

	args=parser.parse_args()

	main(args.input, args.magnitudes, args.output, args.lca, args.hits, args.virus)

