#!/usr/bin/env python
import pandas as pd
import numpy as np
from collections import Counter

df = pd.read_csv("data/mct-v3/taxa/otutable-subset-n50000-s10-norm-no-soylent-no-transition-pre_L7.txt", skiprows=1, sep="\t", index_col=0)
df_tax = pd.read_csv("results/taxtable-amgut2-filt.txt", sep="\t")

tax_map = {}
for i, row in df_tax.iterrows():
    row = [v for v in row if len(v) > 3]
    tax_map['; '.join(row)] = i

tax = df.index
c = []
for row in tax:
    split = row.split(';')
    if split[-1][0] == 's':
        split[-1] = 's__' + split[-1][3:].replace('_', ' ')
        row = ';'.join(split)
    name = '; '.join([v for v in row.split(';') if not 'Other' in v])
    if name in tax_map:
        c.append(name)

print(len(c))
# too few map, only 7
