#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ШАГ 5: СОЗДАНИЕ BIGWIG ТРЕКОВ (ЧЕРЕЗ BEDGRAPH + BEDGRAPHTOBIGWIG)
"""

import pandas as pd
import numpy as np  # ← ДОБАВЛЕНО!
import subprocess
import os
from pathlib import Path

print("=" * 70)
print("ШАГ 5: СОЗДАНИЕ BIGWIG ТРЕКОВ (ЧЕРЕЗ BEDGRAPH)")
print("=" * 70)

input_file = '/home/dasha/data/bismark/MoPh7_cpg_methylation_values.tsv.gz'
chrom_sizes_file = '/home/dasha/data/bismark/chrom.sizes'
bedgraph_to_bigwig = '/home/dasha/Desktop/bedGraphToBigWig'

if not Path(bedgraph_to_bigwig).exists():
    print(f"❌ bedGraphToBigWig не найден по пути: {bedgraph_to_bigwig}")
    exit()
    
print("\n1. Загружаем chrom.sizes...")
chrom_lengths = {}
with open(chrom_sizes_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            chrom_lengths[parts[0]] = int(parts[1])
print(f"   Загружено {len(chrom_lengths)} хромосом")


Path("results/tracks").mkdir(parents=True, exist_ok=True)

# 3. Треки
tracks = [
    {'col': 'beta', 'name': 'beta_methylation'},
    {'col': 'm_value', 'name': 'm_value'},
    {'col': 'coverage', 'name': 'coverage'}
]

print("\n2. Создаем bedGraph и bigWig треки...")

for track in tracks:
    print(f"\n   Создаем {track['name']} трек...")
    
    bedgraph_file = f'results/tracks/MoPh7_{track["name"]}.bedGraph'
    sorted_bedgraph = f'results/tracks/MoPh7_{track["name"]}.sorted.bedGraph'
    bigwig_file = f'results/tracks/MoPh7_{track["name"]}.bw'
    
    for f in [bedgraph_file, sorted_bedgraph, bigwig_file]:
        if Path(f).exists():
            Path(f).unlink()

    with open(bedgraph_file, 'w') as out:
        chunk_size = 500000
        total_processed = 0
        
        for chunk in pd.read_csv(
            input_file,
            sep='\t',
            compression='gzip',
            usecols=['chrom', 'start', 'end', track['col']],
            chunksize=chunk_size
        ):

            chunk['end'] = chunk['start'] + 1
            
            # Фильтруем по chrom.sizes
            for chrom in chunk['chrom'].unique():
                if chrom not in chrom_lengths:
                    continue
                
                chrom_data = chunk[chunk['chrom'] == chrom].copy()
                max_end = chrom_lengths[chrom]
                
         
                chrom_data = chrom_data[chrom_data['end'] <= max_end]
                chrom_data = chrom_data[chrom_data['start'] < chrom_data['end']]
                chrom_data = chrom_data[chrom_data['start'] >= 0]
                chrom_data = chrom_data[~chrom_data[track['col']].isna()]
                
                if len(chrom_data) == 0:
                    continue

                for _, row in chrom_data.iterrows():
 
                    val = row[track['col']]
                    if pd.isna(val) or np.isinf(val):
                        val = 0.0
                    out.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{float(val)}\n")
            
            total_processed += len(chunk)
            if total_processed % 5_000_000 == 0:
                print(f"      Обработано {total_processed:,} строк...")
    
    print(f"       Сохранен bedGraph: {bedgraph_file}")
    
    print(f"      Сортируем...")
    sort_cmd = f"sort -k1,1 -k2,2n {bedgraph_file} > {sorted_bedgraph}"
    subprocess.run(sort_cmd, shell=True, check=True)
    print(f"      Отсортирован: {sorted_bedgraph}")

    print(f"      Конвертируем в bigWig...")
    convert_cmd = f"{bedgraph_to_bigwig} {sorted_bedgraph} {chrom_sizes_file} {bigwig_file}"
    result = subprocess.run(convert_cmd, shell=True)
    
    if result.returncode == 0:
        file_size = Path(bigwig_file).stat().st_size / (1024 * 1024)
        print(f"      Создан bigWig: {bigwig_file} ({file_size:.1f} MB)")
    else:
        print(f"      Ошибка конвертации")
        print(f"      Последние 5 строк bedGraph:")
        os.system(f"tail -5 {bedgraph_file}")
    
    if Path(sorted_bedgraph).exists():
        os.remove(sorted_bedgraph)
    if Path(bedgraph_file).exists():
        os.remove(bedgraph_file)

for track in tracks:
    f = f'results/tracks/MoPh7_{track["name"]}.bw'
    if Path(f).exists() and Path(f).stat().st_size > 0:
        size = Path(f).stat().st_size / (1024 * 1024)
        print(f"   {f} ({size:.1f} MB)")
    else:
        print(f"  проблема!")