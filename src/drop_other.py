#!/usr/bin/env python
import os
import csv
import glob
import pandas as pd

if not os.path.exists("results\\taxa"):
    os.mkdir("results\\taxa")

for file in glob.glob("data\\mct-v3\\taxa\\*.txt"):
    df = pd.read_csv(file, sep='\t', index_col=0, skiprows=1)
    df = df[~df.index.str.contains("Other")]
    outf = file[:-4] + "-no-other.txt"
    outf = os.path.basename(outf)
    outf = os.path.join("results", "taxa", outf)
    df.to_csv(outf)
