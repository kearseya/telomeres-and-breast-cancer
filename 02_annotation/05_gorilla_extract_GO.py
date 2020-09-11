import sys
import os
import pandas as pd
from bs4 import BeautifulSoup
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

input_dir = "/path/to/project/GORILLA-out"
os.chdir(input_dir)
files = list(file for file in os.listdir(input_dir)
         if os.path.isfile(os.path.join(input_dir, file)) and file.endswith('.html'))

def dict_of_tables(files):
    data = {}
    for f in files:
        test_open = open(f, "r")
        test = BeautifulSoup(test_open, "lxml")

        try:
            table = test.find('table', {"id": "table1"})

            go_table = pd.DataFrame(columns=range(0,6), index = range(0, len(table.findChildren(['th', 'tr']))))#[os.path.splitext(test_file)[0]]) # I know the size

            row_marker = 0
            for row in table.find_all('tr'):
                column_marker = 0
                columns = row.find_all('td')
                for index, column in enumerate(columns):
                    if index != 5:
                        go_table.iat[row_marker,column_marker] = column.get_text().replace('\n', ' ').strip()
                        column_marker += 1
                    elif index == 5 and row_marker == 0:
                        go_table.iat[row_marker,column_marker] = column.get_text().replace('\n', ' ').strip()
                        column_marker += 1
                    elif index == 5 and row_marker != 0:
                        gene_dict = {}
                        gene_list = [g for g in column.get_text().splitlines() if g not in ["", " ", "[+] Show genes"]]
                        for gene in gene_list:
                            sep = list(gene.split("-", 1))
                            gene_dict[sep[0].strip()] = sep[-1].strip()
                        go_table.iat[row_marker,column_marker] = gene_dict
                        column_marker += 1
                row_marker += 1

            new_header = go_table.iloc[0]
            go_table = go_table[1:]
            go_table.columns = new_header

            perc_area = str(test.select('h1')[0])
            perc = re.search(r"Note that only (\d+\.\d+)", perc_area)

            number_area = str(test.select('p')[5])
            used = re.search(r"recognized (\d+)", number_area)
            total = re.search(r"genes out of (\d+)", number_area)

            go_table.insert(6, "perc_genes_used", float(perc.group(1).strip()))
            go_table.insert(7, "used", int(used.group(1).strip()))
            go_table.insert(8, "total", int(total.group(1).strip()))
            go_table.insert(9, "sample", os.path.splitext(f)[0])

            go_table.columns = go_table.columns.str.replace(" ", "_")

            data[os.path.splitext(f)[0]] = go_table
        except:
            continue
    return data

data = dict_of_tables(files)

#print(data)

#table of all, unique GO ids and list of samples in column

def combine_table(data):
    df = pd.concat(data.values())
    small = pd.DataFrame(columns=["GO_term", "sample"])
    #df.columns = df.columns.str.replace(" ", "_")
    all_go_terms = df.GO_term.unique()
    for g in all_go_terms:
        samples_in = pd.DataFrame({"GO_term": g, "sample": set(df.loc[df["GO_term"] == g]["sample"])})
        small = pd.concat([small, samples_in], ignore_index = True)
    small = small.drop_duplicates(subset="GO_term").reset_index()
    small["frequency"] = small["sample"].str.len()
    small = small.sort_values("frequency", ascending = False)
    #print(small.loc[small["frequency"] >1 ])
    #print(len(all_go_terms))
    #print(small[small["frequency"] == small["frequency"].max()]["GO_term"].iloc[0])
    #print(small["frequency"].max())
    #print(df)
    print(small)
    return small

small = combine_table(data)



def read_tel_len():
    sample_pairs = pd.read_csv("/path/to/project/sample_pairs.csv")
    return sample_pairs

sample_pairs = read_tel_len()



def plot_go_bar(min_num_samples=2, write_csv=False):
    #fig = plt.figure()
    sample_in_multi = sorted(set(s for n in list(small.loc[small["frequency"] >= min_num_samples ]["sample"]) for s in n))
    #print(sample_in_multi)
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, len(sample_in_multi)))
    melted = pd.DataFrame(columns=sample_in_multi, index=small[small["frequency"] >= min_num_samples]["GO_term"])
    for g in melted.index.values.tolist():
        for s in sample_in_multi:
            if s in small.loc[small["GO_term"] == g]["sample"].iloc[0]:
                melted.loc[g][s] = 1
            else:
                melted.loc[g][s] = 0

    if write_csv == True:
        melted.to_csv("sample_GO_terms_10.csv")

    melted.plot.bar(stacked=True, color=colors)
    plt.show()

plot_go_bar(2, write_csv=False)
