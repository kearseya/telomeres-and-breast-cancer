#!/bin/python3

import os
import csv
import pandas
from pandas import DataFrame
from statistics import mean
import numpy as np
from itertools import chain


cancer_data = pandas.read_csv("data/COSMIC_breast.tsv", sep="\t")
cancer = set(cancer_data['Gene name'])

sample_dirs = os.listdir("snpeff-out")
#print(sample_dirs)


def get_both():
	global data_dict
	global both_dict
	global both_counts
	global both_percent

	data_dict = {}
	both_dict = {}
	both_counts = {}
	both_percent = {}

	for sample in sample_dirs:
		file_open = open("snpeff-out/"+sample+"/gene_list_extract_10", "r").readlines()
		data = set()
		for i in range(len(file_open)):
		    data.add(file_open[i].strip())

		data_dict[sample] = data

		both = [g for g in data if g in cancer]
		both_dict[sample] = both

		both_counts[sample] = len(both)
		#print(sample)
		#print("Data:	", len(data))
		#print("Cancer:	", len(cancer))
		try:
			#print("In both:	", len(both))
			#print("Percent:	", (len(both)/len(data)*100))
			both_percent[sample] = (len(both)/len(data)*100)
		except:
			#print("In both:	", 0)
			#print("Percent:	", 0)
			both_percent[sample] = 0

		#print("")

	#print(both_counts)
	#print(both_percent)

	#print("mean:")
	#print(mean(both_percent[k] for k in both_percent))

get_both()



def make_both_value_table():
	cancer_total_dict = {}
	cancer_total_dict["cancer"] = {}
	cancer_total_dict["total"] = {}
	for sample in sample_dirs:
		cancer_total_dict["cancer"][sample] = len(both_dict[sample])
		cancer_total_dict["total"][sample] = len(data_dict[sample])

	cancer_total_table = pandas.DataFrame(cancer_total_dict["cancer"].values(), cancer_total_dict["cancer"].keys())
	cancer_total_table["total"] = cancer_total_dict["total"].values()
	cancer_total_table.index.name = "sample"
	cancer_total_table.columns = ["cancer", "total"]
	cancer_total_table.to_csv("cancer_total_10.csv")
	print(cancer_total_table)

make_both_value_table()




def sample_count_all(write_csv=False):
	global all_genes
	all_genes = set(chain.from_iterable(list(data_dict.values())))

	gene_dict = {}
	gene_count_dict = {}

	for g in all_genes:
		gene_dict[g] = set()
		for s in data_dict.keys():
			if g in data_dict[s]:
				gene_dict[g].add(s)

	for g in gene_dict:
		gene_count_dict[g] = len(gene_dict[g])

	gene_count_dict = {k: v for k, v in sorted(gene_count_dict.items(), key=lambda item: item[1])}
	#print(gene_count_dict)
	#print(len(all_genes))

	gene_table = pandas.DataFrame(gene_count_dict.values(), gene_count_dict.keys())
	gene_table.index.name = "gene"
	gene_table.columns = ["count"]

	for s in sorted(sample_dirs):
		gene_table[s] = 0
		for g in all_genes:
			if g in data_dict[s]:
				gene_table.loc[g][s] = 1

	gene_table = gene_table.drop(['Gene_Name', 'Percent_of_transcripts_affected\'">'])

	gene_table["COSMIC_breast"] = 0
	for g in all_genes:
		if g in cancer:
			gene_table.loc[g]["COSMIC_breast"] = 1

	print(gene_table)

	if write_csv == True:
		gene_table.to_csv("sample_gene_table_10.csv")

sample_count_all(write_csv=True)
