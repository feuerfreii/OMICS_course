#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 20:50:13 2026

@author: dasha
"""

import pandas as pd
import numpy as np
from pathlib import Path

input_file = '/home/dasha/data/bismark/MoPh7_beta.tsv.gz'
output_file = '/home/dasha/data/bismark/MoPh7_cpg_methylation_values.tsv.gz'

if not Path(input_file).exists():
    exit()

chunk_size = 100000
first_chunk = True
total_processed = 0

sum_m = 0
sum_m_sq = 0
min_m = float('inf')
max_m = float('-inf')
m_values = []

for chunk in pd.read_csv(
    input_file,
    sep='\t',
    compression='gzip',
    chunksize=chunk_size
):

    chunk['m_value'] = np.log2(
        (chunk['count_methylated'] + 1) / 
        (chunk['count_unmethylated'] + 1)
    )
    
    # Собираем статистику
    n = len(chunk)
    total_processed += n
    
    sum_m += chunk['m_value'].sum()
    sum_m_sq += (chunk['m_value'] ** 2).sum()
    
    chunk_min = chunk['m_value'].min()
    chunk_max = chunk['m_value'].max()
    if chunk_min < min_m:
        min_m = chunk_min
    if chunk_max > max_m:
        max_m = chunk_max
    
    if len(m_values) < 10000:
        m_values.extend(chunk['m_value'].sample(min(100, len(chunk))).tolist())
    
    # Сохраняем
    chunk.to_csv(
        output_file,
        sep='\t',
        mode='a',
        header=first_chunk,
        index=False,
        compression='gzip'
    )
    first_chunk = False
    
    if total_processed % 1_000_000 == 0:
        print(f"   Обработано {total_processed:,} строк...")


mean_m = sum_m / total_processed
print(f"   Среднее: {mean_m:.4f}")

variance = (sum_m_sq / total_processed) - (mean_m ** 2)
std_m = np.sqrt(variance)
print(f"   Стандартное отклонение: {std_m:.4f}")

print(f"   Min: {min_m:.4f}")
print(f"   Max: {max_m:.4f}")

m_values = np.array(m_values)
print(f"   Медиана (по выборке): {np.median(m_values):.4f}")
print(f"   25-й перцентиль: {np.percentile(m_values, 25):.4f}")
print(f"   75-й перцентиль: {np.percentile(m_values, 75):.4f}")
