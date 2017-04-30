#!/usr/bin/env python
import pandas as pd
import numpy as np

df = pd.read_csv("data/mct-v3/otutable-subset.txt", skiprows=1, sep="\t")

taxa = df['taxonomy']
df = df.drop('taxonomy', 1)

df_taxa = pd.DataFrame([[name[3:] for name in row.split('; ')] for row in taxa])\

df_taxa.index = df['#OTU ID']
df_taxa.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

df_taxa.to_csv("results/taxatable-subset.txt", sep="\t")

df = pd.read_csv("data/mct-v3/foodtaxa/foodtable_for_taxonomy-subset_L3.txt", skiprows=1, sep="\t")

taxa = df['#OTU ID']

df_taxa = pd.DataFrame([[name[3:] for name in row.split(';')] for row in taxa])\

df_taxa.index = taxa
df_taxa.columns = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']

df_taxa.to_csv("results/taxatable-foodtable_for_taxonomy-subset_L3.txt", sep="\t")
