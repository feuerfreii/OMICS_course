#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 17:52:06 2026

@author: dasha
"""

import pandas as pd
import numpy as np

chunk_size = 100000
chunks = []
total_rows = 0

for chunk in pd.read_csv(
    'data/bismark/MoPh7_1_bismark_bt2_pe.bismark.cov.gz',
    sep='\t',
    header=None,
    names=['chrom', 'start', 'end', 'methylation_percentage', 
           'count_methylated', 'count_unmethylated'],
    chunksize=chunk_size
):
    chunk['coverage'] = chunk['count_methylated'] + chunk['count_unmethylated']
    chunk_filtered = chunk[chunk['coverage'] >= 5]
    if len(chunk_filtered) > 0:
        chunks.append(chunk_filtered)
    total_rows += len(chunk)
    print(f"   Обработано {total_rows:,} строк...")

df_coverage = pd.concat(chunks, ignore_index=True)
df_coverage.to_csv('data/bismark/MoPh7_coverage.tsv.gz', sep='\t', index=False, compression='gzip')
