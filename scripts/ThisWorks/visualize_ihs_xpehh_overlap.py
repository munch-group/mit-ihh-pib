#!/usr/bin/env python3
"""
Create visualizations showing where iHS and XP-EHH signals overlap and where they don't
"""

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path
import numpy as np

# Setup paths
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
comparison_dir = base_dir / "plots/ihs_xpehh_comparison"
output_dir = comparison_dir

print("=" * 80)
print("Visualizing iHS vs XP-EHH Signal Overlap")
print("=" * 80)
print()

# Load data
print("Loading data...")
high_conf = pd.read_csv(comparison_dir / "high_confidence_regions_EUR_vs_EAS.tsv", sep='\t')
ihs_only = pd.read_csv(comparison_dir / "ihs_only_regions.tsv", sep='\t')
xpehh_only = pd.read_csv(comparison_dir / "xpehh_only_regions_EUR_vs_EAS.tsv", sep='\t')

print(f"  High confidence (both): {len(high_conf)} regions")
print(f"  iHS only: {len(ihs_only)} regions")
print(f"  XP-EHH only: {len(xpehh_only)} regions")
print()

# Get unique iHS and XP-EHH regions from high confidence
# Extract unique regions (deduplicate by region boundaries)
ihs_regions_dict = {}
xpehh_regions_dict = {}

for _, row in high_conf.iterrows():
    # Parse ihs_region and xpehh_region columns
    ihs_range = row['ihs_region'].split('-')
    xpehh_range = row['xpehh_region'].split('-')

    ihs_start, ihs_end = int(ihs_range[0]), int(ihs_range[1])
    xpehh_start, xpehh_end = int(xpehh_range[0]), int(xpehh_range[1])

    # Use region boundaries as key to deduplicate
    ihs_key = (ihs_start, ihs_end)
    xpehh_key = (xpehh_start, xpehh_end)

    # Keep the max score for each unique region
    if ihs_key not in ihs_regions_dict or abs(row['ihs_max']) > abs(ihs_regions_dict[ihs_key][2]):
        ihs_regions_dict[ihs_key] = (ihs_start, ihs_end, row['ihs_max'])

    if xpehh_key not in xpehh_regions_dict or abs(row['xpehh_max']) > abs(xpehh_regions_dict[xpehh_key][2]):
        xpehh_regions_dict[xpehh_key] = (xpehh_start, xpehh_end, row['xpehh_max'])

ihs_regions_high_conf = list(ihs_regions_dict.values())
xpehh_regions_high_conf = list(xpehh_regions_dict.values())

print(f"Unique iHS regions in high confidence: {len(ihs_regions_high_conf)}")
print(f"Unique XP-EHH regions in high confidence: {len(xpehh_regions_high_conf)}")
print(f"Note: {len(high_conf)} overlap instances from these unique regions")
print()

# ============================================================================
# Visualization 1: Genome-wide overview showing all signal types
# ============================================================================
print("Creating Visualization 1: Genome-wide overview...")

fig, ax = plt.subplots(figsize=(16, 8))

# Plot each type of region
y_offset = 0

# High confidence regions (overlapping) - plot both iHS and XP-EHH
for start, end, score in ihs_regions_high_conf:
    ax.barh(y_offset, end - start, left=start, height=0.8,
            color='darkgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    y_offset += 1

# Add gap
y_offset += 2

# XP-EHH from high confidence (to show overlap)
for start, end, score in xpehh_regions_high_conf:
    ax.barh(y_offset, end - start, left=start, height=0.8,
            color='darkblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    y_offset += 1

# Add gap
y_offset += 2

# iHS only regions
for _, row in ihs_only.iterrows():
    ax.barh(y_offset, row['end'] - row['start'], left=row['start'], height=0.8,
            color='orange', alpha=0.7, edgecolor='black', linewidth=0.5)
    y_offset += 1

# Add gap
y_offset += 2

# XP-EHH only regions
for _, row in xpehh_only.iterrows():
    ax.barh(y_offset, row['end'] - row['start'], left=row['start'], height=0.8,
            color='red', alpha=0.7, edgecolor='black', linewidth=0.5)
    y_offset += 1

ax.set_xlabel('Position on Chromosome X (bp)', fontsize=12, fontweight='bold')
ax.set_ylabel('Selection Signal Regions', fontsize=12, fontweight='bold')
ax.set_title('iHS vs XP-EHH Selection Signal Comparison\nChromosome X',
             fontsize=14, fontweight='bold')

# Format x-axis
ax.ticklabel_format(style='plain', axis='x')
ax.set_xlim(0, 160000000)  # X chromosome length

# Create legend
legend_elements = [
    mpatches.Patch(facecolor='darkgreen', edgecolor='black', alpha=0.7,
                   label=f'iHS (High Confidence, n={len(ihs_regions_high_conf)})'),
    mpatches.Patch(facecolor='darkblue', edgecolor='black', alpha=0.7,
                   label=f'XP-EHH (High Confidence, n={len(xpehh_regions_high_conf)})'),
    mpatches.Patch(facecolor='orange', edgecolor='black', alpha=0.7,
                   label=f'iHS Only (n={len(ihs_only)})'),
    mpatches.Patch(facecolor='red', edgecolor='black', alpha=0.7,
                   label=f'XP-EHH Only (n={len(xpehh_only)})')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

# Remove y-ticks (not meaningful for stacked regions)
ax.set_yticks([])

plt.tight_layout()
outfile = output_dir / "signal_overlap_genome_wide.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# Visualization 2: Summary statistics (corrected for unique regions)
# ============================================================================
print("Creating Visualization 2: Summary statistics...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# LEFT PANEL: Unique regions
n_ihs_total = len(ihs_regions_high_conf) + len(ihs_only)
n_xpehh_total = len(xpehh_regions_high_conf) + len(xpehh_only)
n_ihs_confirmed = len(ihs_regions_high_conf)
n_xpehh_confirmed = len(xpehh_regions_high_conf)

categories1 = ['Total iHS\nRegions', 'iHS Regions\nConfirmed', 'iHS Only\nRegions']
counts1 = [n_ihs_total, n_ihs_confirmed, len(ihs_only)]
colors1 = ['lightcoral', 'darkgreen', 'orange']

bars1 = ax1.bar(categories1, counts1, color=colors1, alpha=0.7, edgecolor='black', linewidth=2)

for bar in bars1:
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height)}',
            ha='center', va='bottom', fontsize=14, fontweight='bold')

ax1.set_ylabel('Number of Regions', fontsize=12, fontweight='bold')
ax1.set_title('iHS Regions: Validation by XP-EHH', fontsize=13, fontweight='bold')
ax1.set_ylim(0, max(counts1) * 1.15)

# Add validation rate
if n_ihs_total > 0:
    validation_rate = (n_ihs_confirmed / n_ihs_total) * 100
    ax1.text(0.5, 0.90, f'Validation Rate: {validation_rate:.1f}%',
            transform=ax1.transAxes, ha='center', va='top',
            fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

# RIGHT PANEL: XP-EHH regions
categories2 = ['Total XP-EHH\nRegions', 'XP-EHH Regions\nConfirmed', 'XP-EHH Only\nRegions']
counts2 = [n_xpehh_total, n_xpehh_confirmed, len(xpehh_only)]
colors2 = ['lightblue', 'darkgreen', 'red']

bars2 = ax2.bar(categories2, counts2, color=colors2, alpha=0.7, edgecolor='black', linewidth=2)

for bar in bars2:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height)}',
            ha='center', va='bottom', fontsize=14, fontweight='bold')

ax2.set_ylabel('Number of Regions', fontsize=12, fontweight='bold')
ax2.set_title('XP-EHH Regions: Validation by iHS', fontsize=13, fontweight='bold')
ax2.set_ylim(0, max(counts2) * 1.15)

# Add validation rate
if n_xpehh_total > 0:
    xpehh_validation_rate = (n_xpehh_confirmed / n_xpehh_total) * 100
    ax2.text(0.5, 0.90, f'Validation Rate: {xpehh_validation_rate:.1f}%',
            transform=ax2.transAxes, ha='center', va='top',
            fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

plt.suptitle('Selection Signal Validation: iHS vs XP-EHH (EUR vs EAS)\n'
             f'Note: {len(high_conf)} total overlap instances between these unique regions',
             fontsize=14, fontweight='bold')
plt.tight_layout()
outfile = output_dir / "signal_overlap_summary.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# Visualization 3: Regional detail view of top signals
# ============================================================================
print("Creating Visualization 3: Top regions detail view...")

# Get top 10 high confidence regions by iHS score
top_regions = high_conf.nlargest(10, 'ihs_max')

fig, axes = plt.subplots(5, 2, figsize=(16, 20))
axes = axes.flatten()

for idx, (_, region) in enumerate(top_regions.iterrows()):
    if idx >= 10:
        break

    ax = axes[idx]

    # Parse region boundaries
    ihs_range = region['ihs_region'].split('-')
    xpehh_range = region['xpehh_region'].split('-')

    ihs_start, ihs_end = int(ihs_range[0]), int(ihs_range[1])
    xpehh_start, xpehh_end = int(xpehh_range[0]), int(xpehh_range[1])

    # Plot regions
    ax.barh(1, ihs_end - ihs_start, left=ihs_start, height=0.5,
            color='orange', alpha=0.7, edgecolor='black', linewidth=2,
            label=f'iHS (max={region["ihs_max"]:.2f})')

    ax.barh(0, xpehh_end - xpehh_start, left=xpehh_start, height=0.5,
            color='steelblue', alpha=0.7, edgecolor='black', linewidth=2,
            label=f'XP-EHH (max={region["xpehh_max"]:.2f})')

    # Highlight overlap region
    overlap_start = region['overlap_start']
    overlap_end = region['overlap_end']
    ax.axvspan(overlap_start, overlap_end, alpha=0.3, color='green', zorder=0)

    # Add peak markers
    ax.plot(region['ihs_peak_pos'], 1, 'r*', markersize=15, markeredgecolor='black',
            markeredgewidth=1, label='iHS peak', zorder=5)
    ax.plot(region['xpehh_peak_pos'], 0, 'r*', markersize=15, markeredgecolor='black',
            markeredgewidth=1, label='XP-EHH peak', zorder=5)

    # Format
    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['XP-EHH', 'iHS'])
    ax.set_xlabel('Position (bp)', fontsize=9)
    ax.set_title(f'Region {idx+1}: X:{overlap_start:,}-{overlap_end:,}\n'
                 f'Overlap: {region["overlap_kb"]:.1f} kb',
                 fontsize=10, fontweight='bold')
    ax.ticklabel_format(style='plain', axis='x')
    ax.legend(loc='upper right', fontsize=7)
    ax.grid(True, alpha=0.3, axis='x')

plt.suptitle('Top 10 High Confidence Regions: iHS and XP-EHH Overlap Detail',
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
outfile = output_dir / "top_regions_detail.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# Visualization 4: Signal strength comparison
# ============================================================================
print("Creating Visualization 4: Signal strength comparison...")

fig, ax = plt.subplots(figsize=(10, 10))

# Scatter plot of iHS vs XP-EHH scores for overlapping regions
ihs_scores = high_conf['ihs_max'].abs()
xpehh_scores = high_conf['xpehh_max'].abs()

scatter = ax.scatter(ihs_scores, xpehh_scores, c=high_conf['overlap_kb'],
                    s=100, alpha=0.6, cmap='viridis', edgecolors='black', linewidth=0.5)

# Add colorbar
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Overlap Size (kb)', fontsize=11, fontweight='bold')

# Add diagonal reference line
max_val = max(ihs_scores.max(), xpehh_scores.max())
ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=2, label='Equal strength')

# Add threshold lines
ax.axvline(2.0, color='orange', linestyle=':', linewidth=2, alpha=0.5, label='iHS threshold (2.0)')
ax.axhline(3.0, color='steelblue', linestyle=':', linewidth=2, alpha=0.5, label='XP-EHH threshold (3.0)')

ax.set_xlabel('|iHS Score|', fontsize=12, fontweight='bold')
ax.set_ylabel('|XP-EHH Score|', fontsize=12, fontweight='bold')
ax.set_title('Signal Strength Comparison: iHS vs XP-EHH\nHigh Confidence Regions',
             fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)

# Add correlation text
correlation = np.corrcoef(ihs_scores, xpehh_scores)[0, 1]
ax.text(0.95, 0.05, f'Correlation: {correlation:.3f}',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout()
outfile = output_dir / "signal_strength_comparison.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# Create summary statistics file
# ============================================================================
print("\nCreating summary statistics file...")

summary = {
    'total_ihs_regions': n_ihs_total,
    'unique_ihs_regions_confirmed': n_ihs_confirmed,
    'ihs_only_regions': len(ihs_only),
    'total_xpehh_regions': n_xpehh_total,
    'unique_xpehh_regions_confirmed': n_xpehh_confirmed,
    'xpehh_only_regions': len(xpehh_only),
    'total_overlap_instances': len(high_conf),
    'ihs_validation_rate_percent': (n_ihs_confirmed / n_ihs_total * 100) if n_ihs_total > 0 else 0,
    'xpehh_validation_rate_percent': (n_xpehh_confirmed / n_xpehh_total * 100) if n_xpehh_total > 0 else 0,
    'mean_overlap_size_kb': high_conf['overlap_kb'].mean(),
    'median_overlap_size_kb': high_conf['overlap_kb'].median(),
    'mean_ihs_score_high_conf': high_conf['ihs_max'].abs().mean(),
    'mean_xpehh_score_high_conf': high_conf['xpehh_max'].abs().mean(),
    'correlation_ihs_xpehh': correlation
}

summary_df = pd.DataFrame([summary]).T
summary_df.columns = ['value']
summary_df.to_csv(output_dir / "overlap_summary_statistics.tsv", sep='\t')
print(f"  Saved: overlap_summary_statistics.tsv")

print("\n" + "=" * 80)
print("Visualization Complete!")
print("=" * 80)
print(f"\nGenerated files:")
print(f"  1. signal_overlap_genome_wide.png - Genome-wide view of all signals")
print(f"  2. signal_overlap_summary.png - Bar chart summary")
print(f"  3. top_regions_detail.png - Detailed view of top 10 regions")
print(f"  4. signal_strength_comparison.png - iHS vs XP-EHH correlation")
print(f"  5. overlap_summary_statistics.tsv - Numerical summary")
print(f"\nKey Findings:")
print(f"  - {n_ihs_total} total unique iHS regions detected")
print(f"  - {n_xpehh_total} total unique XP-EHH regions detected")
print(f"  - {n_ihs_confirmed} iHS regions confirmed by XP-EHH overlap")
print(f"  - {n_xpehh_confirmed} XP-EHH regions confirmed by iHS overlap")
print(f"  - {len(high_conf)} total overlap instances (one region can have multiple overlaps)")
print(f"  - {(n_ihs_confirmed / n_ihs_total * 100):.1f}% of iHS regions validated by XP-EHH")
print(f"  - Correlation between iHS and XP-EHH scores: {correlation:.3f}")
