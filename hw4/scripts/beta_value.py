import pandas as pd
import numpy as np
from pathlib import Path


input_file = '/home/dasha/data/bismark/MoPh7_coverage.tsv.gz'
output_file = '/home/dasha/data/bismark/MoPh7_beta.tsv.gz'

# Проверяем входной файл
if not Path(input_file).exists():
    print(f"❌ Файл {input_file} не найден!")
    exit()

print(f"\n1. Обрабатываем файл по частям...")

chunk_size = 100000
first_chunk = True
total_processed = 0

sum_beta = 0
sum_beta_sq = 0 
min_beta = float('inf')
max_beta = float('-inf')
beta_values = [] 

print("   Обработка чанков...")

for chunk in pd.read_csv(
    input_file,
    sep='\t',
    compression='gzip',
    chunksize=chunk_size
):

    chunk['beta'] = chunk['count_methylated'] / chunk['coverage']

    n = len(chunk)
    total_processed += n

    sum_beta += chunk['beta'].sum()
    sum_beta_sq += (chunk['beta'] ** 2).sum()

    chunk_min = chunk['beta'].min()
    chunk_max = chunk['beta'].max()
    if chunk_min < min_beta:
        min_beta = chunk_min
    if chunk_max > max_beta:
        max_beta = chunk_max
    
    if len(beta_values) < 10000:
        beta_values.extend(chunk['beta'].sample(min(100, len(chunk))).tolist())
    
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

mean_beta = sum_beta / total_processed
print(f"   Среднее: {mean_beta:.4f}")
variance = (sum_beta_sq / total_processed) - (mean_beta ** 2)
std_beta = np.sqrt(variance)
print(f"   Стандартное отклонение: {std_beta:.4f}")
print(f"   Min: {min_beta:.4f}")
print(f"   Max: {max_beta:.4f}")

beta_values = np.array(beta_values)
print(f"   Медиана (по выборке): {np.median(beta_values):.4f}")
print(f"   25-й перцентиль: {np.percentile(beta_values, 25):.4f}")
print(f"   75-й перцентиль: {np.percentile(beta_values, 75):.4f}")
print(f"\nСохранено: {output_file}")
