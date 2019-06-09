#!/usr/bin/env python

import pandas as pd
import re

df = pd.read_csv("result.csv", sep=',', index_col = [1])

# make column for individual ids
pattern = "([A-Z]{2}[0-9]+)"
iids = []
for vcf in df.index.tolist():
    iid = re.search(pattern, vcf)
    iids.append(iid.group())

df['IID'] = iids

# reorder columns
cols = list(df.columns.values)
new_cols = []
main_cols = ['FID','IID','PAT','MAT','SEX']
for col in main_cols:
    cols.remove(col)
    new_cols.append(col)

for col in cols:
    new_cols.append(col)

df = df[new_cols]

df.to_csv('sample.phe', sep=' ', index=False)
