#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

input_file = '/home/dasha/data/bismark/MoPh7_cpg_methylation_values.tsv.gz'

if not Path(input_file).exists():
    exit()

Path("results/figures").mkdir(parents=True, exist_ok=True)

chunk_size = 200000
total_rows = 0

coverage_bins = np.arange(0, 201, 1)
coverage_hist = np.zeros(len(coverage_bins) - 1)

beta_bins = np.arange(0, 1.01, 0.01)
beta_hist = np.zeros(len(beta_bins) - 1)

m_bins = np.arange(-10, 10.1, 0.1)
m_hist = np.zeros(len(m_bins) - 1)

sum_beta = 0
sum_beta_sq = 0
sum_m = 0
sum_m_sq = 0
sum_coverage = 0
sum_coverage_sq = 0

min_beta = float('inf')
max_beta = float('-inf')
min_m = float('inf')
max_m = float('-inf')
min_cov = float('inf')
max_cov = float('-inf')    


for chunk in pd.read_csv(
    input_file,
    sep='\t',
    compression='gzip',
    chunksize=chunk_size
):
    n = len(chunk)
    total_rows += n
    
    cov_vals = chunk['coverage'].values
    beta_vals = chunk['beta'].values
    m_vals = chunk['m_value'].values

    hist_cov, _ = np.histogram(cov_vals, bins=coverage_bins)
    coverage_hist += hist_cov
    
    hist_beta, _ = np.histogram(beta_vals, bins=beta_bins)
    beta_hist += hist_beta
    
    hist_m, _ = np.histogram(m_vals, bins=m_bins)
    m_hist += hist_m

    sum_beta += beta_vals.sum()
    sum_beta_sq += (beta_vals ** 2).sum()
    sum_m += m_vals.sum()
    sum_m_sq += (m_vals ** 2).sum()
    sum_coverage += cov_vals.sum()
    sum_coverage_sq += (cov_vals ** 2).sum()

    if cov_vals.min() < min_cov:
        min_cov = cov_vals.min()
    if cov_vals.max() > max_cov:
        max_cov = cov_vals.max()

    if beta_vals.min() < min_beta:
        min_beta = beta_vals.min()
    if beta_vals.max() > max_beta:
        max_beta = beta_vals.max()

    if m_vals.min() < min_m:
        min_m = m_vals.min()
    if m_vals.max() > max_m:
        max_m = m_vals.max()
    
    if total_rows % 2_000_000 == 0:
        print(f"   Обработано {total_rows:,} строк...")

mean_beta = sum_beta / total_rows
std_beta = np.sqrt((sum_beta_sq / total_rows) - (mean_beta ** 2))
mean_m = sum_m / total_rows
std_m = np.sqrt((sum_m_sq / total_rows) - (mean_m ** 2))
mean_cov = sum_coverage / total_rows
std_cov = np.sqrt((sum_coverage_sq / total_rows) - (mean_cov ** 2))

print(f"   Coverage:")
print(f"     - Среднее: {mean_cov:.1f}")
print(f"     - Std: {std_cov:.1f}")
print(f"     - Min: {min_cov:.0f}")    
print(f"     - Max: {max_cov:.0f}")     
print(f"   Beta-value:")
print(f"     - Среднее: {mean_beta:.4f}")
print(f"     - Std: {std_beta:.4f}")
print(f"     - Min: {min_beta:.4f}")
print(f"     - Max: {max_beta:.4f}")
print(f"   M-value:")
print(f"     - Среднее: {mean_m:.4f}")
print(f"     - Std: {std_m:.4f}")
print(f"     - Min: {min_m:.4f}")
print(f"     - Max: {max_m:.4f}")

print("\n3. Строим QC-графики...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('WGBS Quality Control - MoPh7 (All Data)', fontsize=18, fontweight='bold')

ax = axes[0, 0]
x_cov = (coverage_bins[:-1] + coverage_bins[1:]) / 2
mask = x_cov <= 100
ax.bar(x_cov[mask], coverage_hist[mask], width=0.8, color='steelblue', edgecolor='black', alpha=0.7)
ax.set_xlabel('Coverage', fontsize=12)
ax.set_ylabel('Number of CpG sites', fontsize=12)
ax.set_title('Coverage distribution (≤ 100x)', fontsize=14)
ax.axvline(mean_cov, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_cov:.1f}')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
x_beta = (beta_bins[:-1] + beta_bins[1:]) / 2
ax.bar(x_beta, beta_hist, width=0.008, color='forestgreen', edgecolor='black', alpha=0.7)
ax.set_xlabel('Beta-value', fontsize=12)
ax.set_ylabel('Number of CpG sites', fontsize=12)
ax.set_title('Beta-value distribution', fontsize=14)
ax.axvline(mean_beta, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_beta:.3f}')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.bar(x_cov[mask], coverage_hist[mask], width=0.8, color='steelblue', edgecolor='black', alpha=0.7)
ax.set_xlabel('Coverage', fontsize=12)
ax.set_ylabel('Number of CpG sites (log scale)', fontsize=12)
ax.set_title('Coverage distribution (≤ 100x, log scale)', fontsize=14)
ax.set_yscale('log')
ax.axvline(mean_cov, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_cov:.1f}')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
x_m = (m_bins[:-1] + m_bins[1:]) / 2
ax.bar(x_m, m_hist, width=0.08, color='darkorange', edgecolor='black', alpha=0.7)
ax.set_xlabel('M-value', fontsize=12)
ax.set_ylabel('Number of CpG sites', fontsize=12)
ax.set_title('M-value distribution', fontsize=14)
ax.axvline(mean_m, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_m:.3f}')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/methylation_qc_plots.png', dpi=150, bbox_inches='tight')

fig2, axes = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle('Cumulative Distribution', fontsize=16, fontweight='bold')

ax = axes[0]
cum_beta = np.cumsum(beta_hist) / np.sum(beta_hist)
ax.plot(x_beta, cum_beta, color='forestgreen', linewidth=2)
ax.set_xlabel('Beta-value', fontsize=12)
ax.set_ylabel('Cumulative fraction', fontsize=12)
ax.set_title('Cumulative distribution of beta-values', fontsize=14)
ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Median')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
cum_cov = np.cumsum(coverage_hist[mask]) / np.sum(coverage_hist)
ax.plot(x_cov[mask], cum_cov, color='steelblue', linewidth=2)
ax.set_xlabel('Coverage', fontsize=12)
ax.set_ylabel('Cumulative fraction', fontsize=12)
ax.set_title('Cumulative distribution of coverage (≤ 100x)', fontsize=14)
ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='Median')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/methylation_cumulative.png', dpi=150, bbox_inches='tight')

plt.show()

print("  results/figures/methylation_qc_plots.png")
print("  results/figures/methylation_cumulative.png")
print(f" CpG сайтов: {total_rows:,}")
print(f" Среднее метилирование: {mean_beta:.2%}")
print(f" Среднее покрытие: {mean_cov:.1f}x")
print(f" Диапазон coverage: {min_cov:.0f} - {max_cov:.0f}")
