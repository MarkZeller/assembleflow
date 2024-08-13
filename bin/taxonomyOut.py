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
		rank_levels=ncbi.get_rank(lineage)

		root=cellular_organisms=superkingdom=phylum=class_=family=superfamily=genus=species="NA"

		for tax in lineage:
			rank=rank_levels.get(tax, "NA")
			if tax==1:
				root=lineage_names.get(tax, "root")
			elif tax==131567:
				cellular_organisms=lineage_names.get(tax, "cellular organisms")
			elif rank=="superkingdom":
				superkingdom=lineage_names.get(tax, "NA")
			elif rank=="phylum":
				phylum=lineage_names.get(tax, "NA")
			elif rank=="class":
				class_=lineage_names.get(tax, "NA")
			elif rank=="superfamily":
				superfamily=lineage_names.get(tax, "NA")
			elif rank=="family":
				family=lineage_names.get(tax, "NA")
			elif rank=="genus":
				genus=lineage_names.get(tax, "NA")
			elif rank=="species":
				species=lineage_names.get(tax, "NA")

		return root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species
	except Exception as e:
		print(f"Error translating taxid {taxid}: {e}")
		return "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"

def parse_tax_info(taxonomy_info):
	ranks=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']
	taxa=taxonomy_info.split(', ')
	return dict(zip(ranks, taxa))

def process_hits(hit, taxonomy_info):
	ranks=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']
	taxa=taxonomy_info.split(', ')
	return dict(zip(ranks, taxa + ['NA']*(len(ranks) - len(taxa))))

def process_virus_hits(hit, taxonomy_info):
	ranks=['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']
	taxa=taxonomy_info.split(', ')
	hit_dict=dict(zip(ranks, taxa+['NA']*(len(ranks)-len(taxa))))

	if hit_dict['Superkingdom']=="Viruses":
		hit_dict['Count']=1
		return hit_dict
	else:
		return None


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


	output_data=[]

	for query, data in query_hits.items():
		hits=data["hits"]
		tax_ids=data["tax_ids"]
		percent_ids=data["percent_ids"]
		lengths=data["lengths"]
		mismatches=data["mismatches"]
		evalues=data["evalues"]
		bitscores=data["bitscores"]

		lca = get_lca(tax_ids)

		if lca is not None:
			try:

				lca_name=ncbi.get_taxid_translator([lca]).get(lca)
				lca_name_id=f"{lca_name} ({lca})"
				root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species=translate_taxid(lca)
				taxonomy_info=', '.join([root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species])

		
				avg_percents=round(sum(percent_ids)/len(percent_ids), 2) if percent_ids else 0
				avg_lengths=round(sum(lengths)/len(lengths), 2) if lengths else 0
				avg_mismatches=round(sum(mismatches)/len(mismatches), 2) if mismatches else 0
				avg_evalues=round(sum(evalues)/len(evalues), 2) if evalues else 0
				avg_bitscores=round(sum(bitscores)/len(bitscores), 2) if bitscores else 0


				output=[
					query,
					', '.join(hits),
					', '.join(tax_ids),
					lca_name_id,
					taxonomy_info,
					avg_percents,
					avg_lengths,
					avg_mismatches,
					avg_evalues,
					avg_bitscores
				]
				output_data.append(output)
			except Exception as e:
				print(f"Error processing LCA {lca} for query {query}: {e}")
		else:
			print(f"No LCA found for query {query}")

		for tax_id in data["tax_ids"]:
			root, cellular_organisms, superkingdom, phylum, class_, superfamily, family, genus, species=translate_taxid(tax_id)
		

	output_df=pd.DataFrame(output_data, columns=["Query", "Hits", "Hit TaxIDs", "LCA", "LCA Taxonomy Info (Root, Cellular Organisms, Superkingdom, Phylum, Class, Superfamily, Family, Genus, Species)", "Av % Identity", "Av Align Length", "Av Mismatches", 
		"Av e-value", "Av Bit Score"])

	output_df['Query']=output_df['Query'].astype(str)
	magnitudes_df['Query']=magnitudes_df['Query'].astype(str)
	output_df['Query']=output_df['Query'].str.strip().str.upper()
	magnitudes_df['Query']=magnitudes_df['Query'].str.strip().str.upper()

	merged_output_df=pd.merge(output_df, magnitudes_df, on='Query', how='inner')
	merged_output_df.to_csv(output_table, index=False)


	lca_summary=[]

	for _, row in output_df.iterrows():
		taxonomy_dict=parse_tax_info(row['LCA Taxonomy Info (Root, Cellular Organisms, Superkingdom, Phylum, Class, Superfamily, Family, Genus, Species)'])
		taxonomy_dict['Count']=1
		lca_summary.append(taxonomy_dict)

	lca_summary_df=pd.DataFrame(lca_summary)
	lca_summary_df=lca_summary_df.groupby(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']).sum().reset_index()
	lca_summary_df=lca_summary_df.sort_values(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])
	lca_summary_df.to_csv(output_lca_summary, sep='\t', index=False)


	hits_summary=[]
	for _, row in output_df.iterrows():
		hits=row['Hits'].split(', ')
		taxonomy_info=row['LCA Taxonomy Info (Root, Cellular Organisms, Superkingdom, Phylum, Class, Superfamily, Family, Genus, Species)']

		for hit in hits:
			hit_dict=process_hits(hit, taxonomy_info)
			hit_dict['Count']=1
			hits_summary.append(hit_dict)

	hits_summary_df=pd.DataFrame(hits_summary)
	hits_summary_df=hits_summary_df.groupby(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']).sum().reset_index()
	hits_summary_df=hits_summary_df.sort_values(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])
	hits_summary_df.to_csv(output_hits_summary, sep='\t', index=False)

	virus_hits_summary=[]
	for _, row in output_df.iterrows():
		hits=row['Hits'].split(', ')
		taxonomy_info=row['LCA Taxonomy Info (Root, Cellular Organisms, Superkingdom, Phylum, Class, Superfamily, Family, Genus, Species)']

		for hit in hits:
			virus_hit_dict=process_virus_hits(hit, taxonomy_info)
			if virus_hit_dict:
				virus_hits_summary.append(virus_hit_dict)

	virus_hits_summary_df=pd.DataFrame(virus_hits_summary)

	if not virus_hits_summary_df.empty:
		virus_hits_summary_df=virus_hits_summary_df.groupby(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species']).sum().reset_index()
		virus_hits_summary_df=virus_hits_summary_df.sort_values(['Root', 'Cellular Organisms', 'Superkingdom', 'Phylum', 'Class', 'Superfamily', 'Family', 'Genus', 'Species'])
		virus_hits_summary_df.to_csv(output_virus_summary, sep='\t', index=False)
	else:
		print("No virus hits found")


if __name__ == "__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('-i', '--input', type=str, required=True)
	parser.add_argument('-m', '--magnitudes', type=str, required=True)
	parser.add_argument('-o', '--output', type=str, required=True)
	parser.add_argument('-l', '--lca', type=str, required=True)
	parser.add_argument('-s', '--hits', type=str, required=True)
	parser.add_argument('-v', '--virus', type=str, required=True)

	args=parser.parse_args()

	main(args.input, args.magnitudes, args.output, args.lca, args.hits, args.virus)

