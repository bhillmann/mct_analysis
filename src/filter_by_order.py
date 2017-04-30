#!/usr/bin/env python
import pandas as pd
import numpy as np
from collections import Counter

order_set = set()
with open('results/amg2.order.txt') as inf:
    for line in inf:
        order_set.add(line.strip())

print(order_set)

df = pd.read_csv("data/mct-v3/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre.txt", skiprows=1, sep="\t")

subset = np.zeros(df.shape[0], dtype=np.bool)
hits = []
for i, taxa in enumerate(df['taxonomy']):
    split = taxa.split('; ')
    if len(split) > 3:
        order = split[3]
        if len(order) > 3:
            if order[3:] in order_set:
                subset[i] = True
                hits.append(order[3:])

df_sub = df.ix[subset ,:]

taxa = df_sub['taxonomy']
df_sub = df_sub.drop('taxonomy', 1)
df_sub['#OTU ID'] = range(df_sub.shape[0])

df_taxa = pd.DataFrame([[name[3:] for name in row.split('; ')] for row in taxa])

df_taxa.columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']

df_sub.to_csv("results/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre-amg2-order.txt", sep="\t", index=False)
df_taxa.to_csv("results/taxatable-subset-n50000-s10-norm-no-soylent-no-transition-pre-amg2-order.txt", sep="\t")
