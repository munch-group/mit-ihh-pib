#!/usr/bin/env python3
"""
Compare two IHS result files and generate comprehensive statistics and visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import sys

# File paths
file1 = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_with_dummy/ALL.chrX.ihs_reformatted.tsv"
file2 = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs/ALL.chrX.ihs.tsv"

print("=" * 80)
print("IHS Results Comparison")
print("=" * 80)
print(f"\nFile 1 (fixed_par1): {file1}")
print(f"File 2 (original):   {file2}")
print()

# Load data
print("Loading data...")
df1 = pd.read_csv(file1, sep="\t")
df2 = pd.read_csv(file2, sep="\t")

print(f"File 1: {len(df1)} variants")
print(f"File 2: {len(df2)} variants")
print()

# Basic overlap analysis
print("=" * 80)
print("Overlap Analysis")
print("=" * 80)
locations1 = set(df1['Location'])
locations2 = set(df2['Location'])

common = locations1 & locations2
only_file1 = locations1 - locations2
only_file2 = locations2 - locations1

print(f"Variants in both files:  {len(common)}")
print(f"Only in file 1:          {len(only_file1)}")
print(f"Only in file 2:          {len(only_file2)}")
print()

# Merge on common locations
df_merged = pd.merge(df1, df2, on='Location', suffixes=('_fixed', '_orig'))
print(f"Merged dataset: {len(df_merged)} variants")
print()

# Compare statistics for each metric
metrics = ['iHH_0', 'iHH_1', 'iHS', 'Std iHS']

print("=" * 80)
print("Statistical Summary")
print("=" * 80)

for metric in metrics:
    print(f"\n{metric}:")
    print("-" * 40)

    col_fixed = f"{metric}_fixed"
    col_orig = f"{metric}_orig"

    # Remove any infinite or NaN values for comparison
    valid_mask = np.isfinite(df_merged[col_fixed]) & np.isfinite(df_merged[col_orig])
    valid_data = df_merged[valid_mask]

    if len(valid_data) == 0:
        print("No valid data for comparison")
        continue

    # Calculate differences
    diff = valid_data[col_fixed] - valid_data[col_orig]
    abs_diff = np.abs(diff)

    # Pearson correlation
    corr, p_value = stats.pearsonr(valid_data[col_fixed], valid_data[col_orig])

    # Summary statistics
    print(f"  Valid pairs:           {len(valid_data)}")
    print(f"  Correlation (r):       {corr:.6f} (p={p_value:.2e})")
    print(f"  Mean difference:       {diff.mean():.6f}")
    print(f"  Std difference:        {diff.std():.6f}")
    print(f"  Median difference:     {diff.median():.6f}")
    print(f"  Mean abs difference:   {abs_diff.mean():.6f}")
    print(f"  Max abs difference:    {abs_diff.max():.6f}")
    print(f"  Q95 abs difference:    {abs_diff.quantile(0.95):.6f}")

    # Count significant differences (e.g., > 1% relative difference)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff = np.abs(diff / valid_data[col_orig]) * 100
        rel_diff = rel_diff[np.isfinite(rel_diff)]

    if len(rel_diff) > 0:
        sig_diff = (rel_diff > 1).sum()
        print(f"  Variants with >1% diff: {sig_diff} ({sig_diff/len(rel_diff)*100:.2f}%)")
        print(f"  Mean relative diff:     {rel_diff.mean():.4f}%")
        print(f"  Median relative diff:   {rel_diff.median():.4f}%")

# Identify variants with largest differences in standardized iHS
print()
print("=" * 80)
print("Variants with Largest Differences in Standardized iHS")
print("=" * 80)

valid_mask = np.isfinite(df_merged['Std iHS_fixed']) & np.isfinite(df_merged['Std iHS_orig'])
df_valid = df_merged[valid_mask].copy()
df_valid['diff_std_ihs'] = np.abs(df_valid['Std iHS_fixed'] - df_valid['Std iHS_orig'])
df_top_diff = df_valid.nlargest(10, 'diff_std_ihs')

print("\nTop 10 variants by absolute difference:")
for idx, row in df_top_diff.iterrows():
    print(f"\n{row['Location']}")
    print(f"  Fixed:    {row['Std iHS_fixed']:8.4f}")
    print(f"  Original: {row['Std iHS_orig']:8.4f}")
    print(f"  Diff:     {row['diff_std_ihs']:8.4f}")

# Create visualizations
print()
print("=" * 80)
print("Generating Visualizations")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('IHS Results Comparison: Fixed PAR1 vs Original', fontsize=16, fontweight='bold')

for idx, metric in enumerate(metrics):
    ax = axes[idx // 2, idx % 2]

    col_fixed = f"{metric}_fixed"
    col_orig = f"{metric}_orig"

    # Filter valid data
    valid_mask = np.isfinite(df_merged[col_fixed]) & np.isfinite(df_merged[col_orig])
    valid_data = df_merged[valid_mask]

    if len(valid_data) == 0:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
        ax.set_title(metric)
        continue

    # Scatter plot with density
    x = valid_data[col_orig]
    y = valid_data[col_fixed]

    # Calculate correlation
    corr, _ = stats.pearsonr(x, y)

    # Plot
    ax.hexbin(x, y, gridsize=50, cmap='viridis', mincnt=1, alpha=0.7)

    # Add diagonal line
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='y=x')

    ax.set_xlabel(f'{metric} (Original)', fontsize=11)
    ax.set_ylabel(f'{metric} (Fixed PAR1)', fontsize=11)
    ax.set_title(f'{metric}\n(r={corr:.4f}, n={len(valid_data)})', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = '/home/vanbruggenmit//mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_with_dummy/comparison_to_OctIHS/ihs_comparison_scatter.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Saved: {output_file}")

# Create difference distribution plots
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
fig.suptitle('Distribution of Differences: Fixed PAR1 - Original', fontsize=16, fontweight='bold')

for idx, metric in enumerate(metrics):
    ax = axes[idx // 2, idx % 2]

    col_fixed = f"{metric}_fixed"
    col_orig = f"{metric}_orig"

    valid_mask = np.isfinite(df_merged[col_fixed]) & np.isfinite(df_merged[col_orig])
    valid_data = df_merged[valid_mask]

    if len(valid_data) == 0:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
        ax.set_title(metric)
        continue

    diff = valid_data[col_fixed] - valid_data[col_orig]

    # Histogram with KDE
    ax.hist(diff, bins=50, alpha=0.6, color='steelblue', edgecolor='black', density=True)

    # Add KDE
    from scipy.stats import gaussian_kde
    try:
        kde = gaussian_kde(diff)
        x_range = np.linspace(diff.min(), diff.max(), 200)
        ax.plot(x_range, kde(x_range), 'r-', linewidth=2, label='KDE')
    except:
        pass

    # Add vertical line at zero
    ax.axvline(0, color='black', linestyle='--', linewidth=2, label='Zero')

    # Statistics
    mean_diff = diff.mean()
    median_diff = diff.median()
    ax.axvline(mean_diff, color='green', linestyle='-', linewidth=2, label=f'Mean={mean_diff:.4f}')
    ax.axvline(median_diff, color='orange', linestyle='-', linewidth=2, label=f'Median={median_diff:.4f}')

    ax.set_xlabel(f'Difference in {metric}', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'{metric}\n(mean={mean_diff:.4f}, std={diff.std():.4f})', fontsize=12)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_with_dummy/comparison_to_OctIHS/ihs_comparison_differences.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Saved: {output_file}")

# Save detailed comparison to CSV
output_csv = '/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_with_dummy/comparison_to_OctIHS/ihs_comparison_detailed_highsig.csv'
df_output = df_valid[['Location', 'iHH_0_fixed', 'iHH_0_orig', 'iHH_1_fixed', 'iHH_1_orig',
                       'iHS_fixed', 'iHS_orig', 'Std iHS_fixed', 'Std iHS_orig']].copy()

for metric in metrics:
    df_output[f'diff_{metric}'] = df_output[f'{metric}_fixed'] - df_output[f'{metric}_orig']
    df_output[f'abs_diff_{metric}'] = np.abs(df_output[f'diff_{metric}'])

df_output = df_output.sort_values('abs_diff_Std iHS', ascending=False)
df_output.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")

print()
print("=" * 80)
print("Comparison Complete!")
print("=" * 80)
