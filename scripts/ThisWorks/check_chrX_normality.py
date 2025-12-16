#!/usr/bin/env python3
"""
Script: Check normality of Std iHS distribution on X chromosome
Purpose: Visualize and test if |Std iHS| follows a normal distribution
         across all 26 sub-populations

Author: MIT IHH PIB Project
Date: 2025-12-16
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")

populations = {
    # African populations
    'YRI': base_dir / "results/ihs_YRI/YRI.chrX.ihs.tsv",
    'LWK': base_dir / "results/ihs_LWK/LWK.chrX.ihs.tsv",
    'GWD': base_dir / "results/ihs_GWD/GWD.chrX.ihs.tsv",
    'MSL': base_dir / "results/ihs_MSL/MSL.chrX.ihs.tsv",
    'ESN': base_dir / "results/ihs_ESN/ESN.chrX.ihs.tsv",
    'ASW': base_dir / "results/ihs_ASW/ASW.chrX.ihs.tsv",
    'ACB': base_dir / "results/ihs_ACB/ACB.chrX.ihs.tsv",
    # European populations
    'CEU': base_dir / "results/ihs_CEU/CEU.chrX.ihs.tsv",
    'TSI': base_dir / "results/ihs_TSI/TSI.chrX.ihs.tsv",
    'FIN': base_dir / "results/ihs_FIN/FIN.chrX.ihs.tsv",
    'GBR': base_dir / "results/ihs_GBR/GBR.chrX.ihs.tsv",
    'IBS': base_dir / "results/ihs_IBS/IBS.chrX.ihs.tsv",
    # East Asian populations
    'CHB': base_dir / "results/ihs_CHB/CHB.chrX.ihs.tsv",
    'JPT': base_dir / "results/ihs_JPT/JPT.chrX.ihs.tsv",
    'CHS': base_dir / "results/ihs_CHS/CHS.chrX.ihs.tsv",
    'CDX': base_dir / "results/ihs_CDX/CDX.chrX.ihs.tsv",
    'KHV': base_dir / "results/ihs_KHV/KHV.chrX.ihs.tsv",
    # South Asian populations
    'GIH': base_dir / "results/ihs_GIH/GIH.chrX.ihs.tsv",
    'PJL': base_dir / "results/ihs_PJL/PJL.chrX.ihs.tsv",
    'BEB': base_dir / "results/ihs_BEB/BEB.chrX.ihs.tsv",
    'STU': base_dir / "results/ihs_STU/STU.chrX.ihs.tsv",
    'ITU': base_dir / "results/ihs_ITU/ITU.chrX.ihs.tsv",
    # American populations
    'MXL': base_dir / "results/ihs_MXL/MXL.chrX.ihs.tsv",
    'PUR': base_dir / "results/ihs_PUR/PUR.chrX.ihs.tsv",
    'CLM': base_dir / "results/ihs_CLM/CLM.chrX.ihs.tsv",
    'PEL': base_dir / "results/ihs_PEL/PEL.chrX.ihs.tsv"
}

output_dir = base_dir / "results/chrX_normality_check"
output_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("X Chromosome iHS Normality Check")
print("=" * 70)
print()

# Load data for all populations
data = {}
for pop, file_path in populations.items():
    print(f"Loading {pop} data from chrX...", end=" ")
    df = pd.read_csv(file_path, sep='\t')
    data[pop] = df
    print(f"{len(df):,} variants")

print()

# Set plot style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

# Define colors for populations by super-population
pop_groups = {
    'AFR': ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB'],
    'EUR': ['CEU', 'TSI', 'FIN', 'GBR', 'IBS'],
    'EAS': ['CHB', 'JPT', 'CHS', 'CDX', 'KHV'],
    'SAS': ['GIH', 'PJL', 'BEB', 'STU', 'ITU'],
    'AMR': ['MXL', 'PUR', 'CLM', 'PEL']
}

# Generate colors using colormap
colors = {}
group_colormaps = {
    'AFR': plt.cm.Greens,
    'EUR': plt.cm.Oranges,
    'EAS': plt.cm.Blues,
    'SAS': plt.cm.Purples,
    'AMR': plt.cm.Reds
}

for group, pops in pop_groups.items():
    cmap = group_colormaps[group]
    for i, pop in enumerate(pops):
        colors[pop] = cmap(0.4 + 0.5 * (i / max(len(pops) - 1, 1)))

print("=" * 70)
print("Normality Statistics")
print("=" * 70)
print()

normality_results = []

for pop, df in data.items():
    std_ihs = df['Std iHS'].values
    abs_std_ihs = np.abs(std_ihs)

    # Calculate statistics
    mean_val = np.mean(std_ihs)
    median_val = np.median(std_ihs)
    std_val = np.std(std_ihs)
    skew_val = stats.skew(std_ihs)
    kurt_val = stats.kurtosis(std_ihs)

    # Normality tests
    # Shapiro-Wilk (use subset if too large)
    if len(std_ihs) > 5000:
        shapiro_stat, shapiro_p = stats.shapiro(np.random.choice(std_ihs, 5000, replace=False))
    else:
        shapiro_stat, shapiro_p = stats.shapiro(std_ihs)

    # Kolmogorov-Smirnov test against standard normal
    ks_stat, ks_p = stats.kstest(std_ihs, 'norm', args=(0, 1))

    # Anderson-Darling test
    anderson_result = stats.anderson(std_ihs, dist='norm')

    print(f"{pop} Population:")
    print(f"  N variants: {len(std_ihs):,}")
    print(f"  Mean: {mean_val:.4f} (expected: 0)")
    print(f"  Median: {median_val:.4f}")
    print(f"  Std Dev: {std_val:.4f} (expected: 1)")
    print(f"  Skewness: {skew_val:.4f} (expected: 0)")
    print(f"  Kurtosis: {kurt_val:.4f} (expected: 0)")
    print(f"  Shapiro-Wilk test: W={shapiro_stat:.4f}, p={shapiro_p:.2e}")
    print(f"  Kolmogorov-Smirnov test: D={ks_stat:.4f}, p={ks_p:.2e}")
    print(f"  Anderson-Darling test: statistic={anderson_result.statistic:.4f}")
    print()

    normality_results.append({
        'population': pop,
        'n_variants': len(std_ihs),
        'mean': mean_val,
        'median': median_val,
        'std': std_val,
        'skewness': skew_val,
        'kurtosis': kurt_val,
        'shapiro_w': shapiro_stat,
        'shapiro_p': shapiro_p,
        'ks_d': ks_stat,
        'ks_p': ks_p,
        'anderson_stat': anderson_result.statistic
    })

# Create detailed plots for each super-population group
for group_name, group_pops in pop_groups.items():
    n_pops = len(group_pops)
    fig = plt.figure(figsize=(5 * n_pops, 12))

    for idx, pop in enumerate(group_pops, 1):
        if pop not in data:
            continue

        df = data[pop]
        std_ihs = df['Std iHS'].values
        abs_std_ihs = np.abs(std_ihs)

        # Plot 1: Histogram with normal overlay (row 1)
        ax1 = plt.subplot(3, n_pops, idx)
        ax1.hist(std_ihs, bins=100, density=True, alpha=0.6, color=colors[pop], label=f'{pop} data')

        # Overlay theoretical normal distribution
        x_range = np.linspace(std_ihs.min(), std_ihs.max(), 100)
        ax1.plot(x_range, stats.norm.pdf(x_range, 0, 1), 'k--', linewidth=2, label='N(0,1)')

        ax1.set_xlabel('Std iHS')
        ax1.set_ylabel('Density')
        ax1.set_title(f'{pop}: Std iHS Distribution')
        ax1.legend(fontsize=8)
        ax1.axvline(0, color='gray', linestyle=':', alpha=0.5)

        # Plot 2: Q-Q plot (row 2)
        ax2 = plt.subplot(3, n_pops, idx + n_pops)
        stats.probplot(std_ihs, dist="norm", plot=ax2)
        ax2.set_title(f'{pop}: Q-Q Plot')
        ax2.get_lines()[0].set_color(colors[pop])
        ax2.get_lines()[0].set_markersize(2)
        ax2.get_lines()[0].set_alpha(0.5)

        # Plot 3: Histogram of |Std iHS| (row 3)
        ax3 = plt.subplot(3, n_pops, idx + 2 * n_pops)
        ax3.hist(abs_std_ihs, bins=100, density=True, alpha=0.6, color=colors[pop], label=f'{pop} data')

        # Overlay theoretical half-normal distribution
        x_range_pos = np.linspace(0, abs_std_ihs.max(), 100)
        ax3.plot(x_range_pos, 2 * stats.norm.pdf(x_range_pos, 0, 1), 'k--', linewidth=2, label='Half-normal')

        ax3.set_xlabel('|Std iHS|')
        ax3.set_ylabel('Density')
        ax3.set_title(f'{pop}: |Std iHS| Distribution')
        ax3.legend(fontsize=8)

        # Add threshold lines
        for thresh in [2.0, 2.5, 3.0]:
            ax3.axvline(thresh, color='red', linestyle='--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    plt.savefig(output_dir / f"chrX_normality_{group_name}.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / f'chrX_normality_{group_name}.png'}")

print()

# Save normality test results
normality_df = pd.DataFrame(normality_results)
normality_df.to_csv(output_dir / "chrX_normality_statistics.tsv", sep='\t', index=False)
print(f"Saved: {output_dir / 'chrX_normality_statistics.tsv'}")
print()

# Create comparison plots by super-population group
for group_name, group_pops in pop_groups.items():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot 1: All populations in group overlaid - Std iHS
    for pop in group_pops:
        if pop not in data:
            continue
        df = data[pop]
        axes[0].hist(df['Std iHS'], bins=100, density=True, alpha=0.4,
                     color=colors[pop], label=pop)
    x_range = np.linspace(-6, 6, 100)
    axes[0].plot(x_range, stats.norm.pdf(x_range, 0, 1), 'k--', linewidth=2, label='N(0,1)')
    axes[0].set_xlabel('Std iHS')
    axes[0].set_ylabel('Density')
    axes[0].set_title(f'{group_name}: Std iHS Distribution')
    axes[0].legend(fontsize=8)
    axes[0].axvline(0, color='gray', linestyle=':', alpha=0.5)

    # Plot 2: All populations in group overlaid - |Std iHS|
    for pop in group_pops:
        if pop not in data:
            continue
        df = data[pop]
        axes[1].hist(np.abs(df['Std iHS']), bins=100, density=True, alpha=0.4,
                     color=colors[pop], label=pop)
    x_range_pos = np.linspace(0, 6, 100)
    axes[1].plot(x_range_pos, 2 * stats.norm.pdf(x_range_pos, 0, 1), 'k--',
                 linewidth=2, label='Half-normal')
    axes[1].set_xlabel('|Std iHS|')
    axes[1].set_ylabel('Density')
    axes[1].set_title(f'{group_name}: |Std iHS| Distribution')
    axes[1].legend(fontsize=8)
    for thresh in [2.0, 2.5, 3.0]:
        axes[1].axvline(thresh, color='red', linestyle='--', alpha=0.3, linewidth=1)

    # Plot 3: Empirical CDFs
    for pop in group_pops:
        if pop not in data:
            continue
        df = data[pop]
        sorted_vals = np.sort(np.abs(df['Std iHS']))
        y_vals = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        axes[2].plot(sorted_vals, y_vals, label=pop, color=colors[pop], linewidth=2)

    # Add theoretical CDF for half-normal
    x_theory = np.linspace(0, 6, 1000)
    y_theory = 2 * stats.norm.cdf(x_theory, 0, 1) - 1
    axes[2].plot(x_theory, y_theory, 'k--', linewidth=2, label='Half-normal')

    axes[2].set_xlabel('|Std iHS|')
    axes[2].set_ylabel('Cumulative Probability')
    axes[2].set_title(f'{group_name}: Empirical CDF Comparison')
    axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3)
    for thresh in [2.0, 2.5, 3.0]:
        axes[2].axvline(thresh, color='red', linestyle='--', alpha=0.3, linewidth=1)

    plt.tight_layout()
    plt.savefig(output_dir / f"chrX_comparison_{group_name}.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / f'chrX_comparison_{group_name}.png'}")

print()

# Create an overall comparison plot with all populations
fig_all, axes_all = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: All populations - Std iHS (use different styles for readability)
for group_name, group_pops in pop_groups.items():
    for pop in group_pops:
        if pop not in data:
            continue
        df = data[pop]
        axes_all[0].hist(df['Std iHS'], bins=100, density=True, alpha=0.2,
                         color=colors[pop], label=pop)
x_range = np.linspace(-6, 6, 100)
axes_all[0].plot(x_range, stats.norm.pdf(x_range, 0, 1), 'k--', linewidth=3, label='N(0,1)')
axes_all[0].set_xlabel('Std iHS')
axes_all[0].set_ylabel('Density')
axes_all[0].set_title('Std iHS Distribution - All 26 Populations')
axes_all[0].axvline(0, color='gray', linestyle=':', alpha=0.5)

# Plot 2: All populations - |Std iHS|
for group_name, group_pops in pop_groups.items():
    for pop in group_pops:
        if pop not in data:
            continue
        df = data[pop]
        axes_all[1].hist(np.abs(df['Std iHS']), bins=100, density=True, alpha=0.2,
                         color=colors[pop])
x_range_pos = np.linspace(0, 6, 100)
axes_all[1].plot(x_range_pos, 2 * stats.norm.pdf(x_range_pos, 0, 1), 'k--',
                 linewidth=3, label='Half-normal')
axes_all[1].set_xlabel('|Std iHS|')
axes_all[1].set_ylabel('Density')
axes_all[1].set_title('|Std iHS| Distribution - All 26 Populations')
axes_all[1].legend(fontsize=8)
for thresh in [2.0, 2.5, 3.0]:
    axes_all[1].axvline(thresh, color='red', linestyle='--', alpha=0.3, linewidth=1)

# Plot 3: Empirical CDFs by super-population group (aggregate)
for group_name, group_pops in pop_groups.items():
    all_vals = []
    for pop in group_pops:
        if pop in data:
            all_vals.extend(np.abs(data[pop]['Std iHS']).tolist())
    if all_vals:
        sorted_vals = np.sort(all_vals)
        y_vals = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        # Use the first population's color for the group
        group_color = colors[group_pops[0]] if group_pops[0] in colors else 'gray'
        axes_all[2].plot(sorted_vals, y_vals, label=group_name, linewidth=2)

x_theory = np.linspace(0, 6, 1000)
y_theory = 2 * stats.norm.cdf(x_theory, 0, 1) - 1
axes_all[2].plot(x_theory, y_theory, 'k--', linewidth=3, label='Half-normal')
axes_all[2].set_xlabel('|Std iHS|')
axes_all[2].set_ylabel('Cumulative Probability')
axes_all[2].set_title('Empirical CDF by Super-Population')
axes_all[2].legend(fontsize=10)
axes_all[2].grid(True, alpha=0.3)
for thresh in [2.0, 2.5, 3.0]:
    axes_all[2].axvline(thresh, color='red', linestyle='--', alpha=0.3, linewidth=1)

plt.tight_layout()
plt.savefig(output_dir / "chrX_comparison_all.png", dpi=300, bbox_inches='tight')
print(f"Saved: {output_dir / 'chrX_comparison_all.png'}")
print()

print("=" * 70)
print("Interpretation")
print("=" * 70)
print()
print("Under the null hypothesis of neutrality, Std iHS should be ~ N(0, 1)")
print()
print("Deviations from normality may indicate:")
print("  - Widespread selection on the X chromosome")
print("  - Population-specific demographic history")
print("  - Technical artifacts in iHS calculation")
print()
print("For the normality tests:")
print("  - p < 0.05 suggests significant deviation from normality")
print("  - Q-Q plots show how well data matches theoretical normal distribution")
print("  - Skewness and kurtosis should be close to 0 for normal distribution")
print()
print("If distributions are non-normal, consider using empirical percentiles")
print("rather than fixed thresholds (|Std iHS| >= 2.5) for defining candidates")
print()
