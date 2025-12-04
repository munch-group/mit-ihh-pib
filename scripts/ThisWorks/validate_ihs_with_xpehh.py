#!/usr/bin/env python3
"""
Validate October iHS findings using XP-EHH cross-population comparisons

Goal: Determine which iHS selection signals from October analysis are:
1. HIGH CONFIDENCE: Validated by XP-EHH (both methods agree)
2. INCOMPLETE SWEEP: iHS signal only (selection in progress, not near fixation)
3. FALSE POSITIVE: iHS signal without cross-population support

For EUR-specific selection:
- iHS shows signal in EUR samples (from October)
- XP-EHH shows positive signal for EUR vs EAS
- XP-EHH shows positive signal for EUR vs AFR
→ Strong evidence for EUR-specific positive selection

Author: MIT IHH PIB Project
Date: 2025-12-04
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from pathlib import Path
import numpy as np

# Setup
base_dir = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib")
comparison_dir = base_dir / "plots/ihs_xpehh_comparison"
output_dir = comparison_dir

print("=" * 80)
print("VALIDATING OCTOBER iHS FINDINGS WITH XP-EHH")
print("=" * 80)
print()
print("Background:")
print("  - iHS analysis (October): identified selection signals in EUR samples")
print("  - XP-EHH analysis: cross-population comparison to validate")
print()
print("Interpretation Framework:")
print("  HIGH CONFIDENCE = iHS + XP-EHH both show signal")
print("    → Strong evidence for selection")
print("  iHS ONLY = iHS shows signal, XP-EHH doesn't")
print("    → Possible incomplete sweep OR false positive")
print("  XP-EHH ONLY = XP-EHH shows signal, iHS doesn't")
print("    → Possible near-fixation (iHS loses power)")
print("=" * 80)
print()

# Load data
print("Loading comparison results...")
high_conf = pd.read_csv(comparison_dir / "high_confidence_regions_EUR_vs_EAS.tsv", sep='\t')
ihs_only = pd.read_csv(comparison_dir / "ihs_only_regions.tsv", sep='\t')
xpehh_only = pd.read_csv(comparison_dir / "xpehh_only_regions_EUR_vs_EAS.tsv", sep='\t')
eur_specific = pd.read_csv(comparison_dir / "EUR_specific_signals.tsv", sep='\t')

print(f"  High confidence overlaps: {len(high_conf)} instances")
print(f"  iHS only regions: {len(ihs_only)}")
print(f"  XP-EHH only regions: {len(xpehh_only)}")
print(f"  EUR-specific signals (positive in both EUR comparisons): {len(eur_specific)}")
print()

# Extract unique regions from high confidence overlaps
print("Extracting unique regions from overlap data...")
ihs_regions_dict = {}
xpehh_eur_eas_dict = {}

for _, row in high_conf.iterrows():
    ihs_range = row['ihs_region'].split('-')
    xpehh_range = row['xpehh_region'].split('-')

    ihs_start, ihs_end = int(ihs_range[0]), int(ihs_range[1])
    xpehh_start, xpehh_end = int(xpehh_range[0]), int(xpehh_range[1])

    ihs_key = (ihs_start, ihs_end)
    xpehh_key = (xpehh_start, xpehh_end)

    if ihs_key not in ihs_regions_dict or abs(row['ihs_max']) > abs(ihs_regions_dict[ihs_key]['score']):
        ihs_regions_dict[ihs_key] = {
            'start': ihs_start, 'end': ihs_end, 'score': row['ihs_max'],
            'peak': row['ihs_peak_pos']
        }

    if xpehh_key not in xpehh_eur_eas_dict or abs(row['xpehh_max']) > abs(xpehh_eur_eas_dict[xpehh_key]['score']):
        xpehh_eur_eas_dict[xpehh_key] = {
            'start': xpehh_start, 'end': xpehh_end, 'score': row['xpehh_max'],
            'peak': row['xpehh_peak_pos']
        }

n_ihs_validated = len(ihs_regions_dict)
n_ihs_only = len(ihs_only)
n_ihs_total = n_ihs_validated + n_ihs_only

print(f"\nUnique iHS regions from October analysis:")
print(f"  Total: {n_ihs_total}")
print(f"  Validated by XP-EHH: {n_ihs_validated} ({n_ihs_validated/n_ihs_total*100:.1f}%)")
print(f"  Not validated (iHS only): {n_ihs_only} ({n_ihs_only/n_ihs_total*100:.1f}%)")
print()

# ============================================================================
# VISUALIZATION 1: iHS Validation Summary
# ============================================================================
print("Creating Visualization 1: iHS Validation Summary...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# LEFT: Pie chart of validation status
sizes = [n_ihs_validated, n_ihs_only]
labels = [f'Validated by XP-EHH\n({n_ihs_validated} regions)',
          f'Not Validated\n({n_ihs_only} regions)']
colors = ['#2ecc71', '#e74c3c']
explode = (0.05, 0)

wedges, texts, autotexts = ax1.pie(sizes, explode=explode, labels=labels, colors=colors,
                                     autopct='%1.1f%%', startangle=90, textprops={'fontsize': 12})
for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
    autotext.set_fontsize(14)

ax1.set_title(f'Validation Status of {n_ihs_total} iHS Regions\n(October Analysis)',
              fontsize=14, fontweight='bold', pad=20)

# RIGHT: Bar chart with interpretation
categories = ['Total iHS\nRegions\n(Oct 2024)',
              'HIGH CONFIDENCE\n(iHS + XP-EHH)',
              'INCOMPLETE/\nFALSE POSITIVE\n(iHS only)']
counts = [n_ihs_total, n_ihs_validated, n_ihs_only]
colors_bar = ['#3498db', '#2ecc71', '#e74c3c']

bars = ax2.bar(categories, counts, color=colors_bar, alpha=0.8, edgecolor='black', linewidth=2)

for bar in bars:
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
             f'{int(height)}',
             ha='center', va='bottom', fontsize=16, fontweight='bold')

ax2.set_ylabel('Number of Regions', fontsize=13, fontweight='bold')
ax2.set_title('Breakdown by Evidence Strength', fontsize=14, fontweight='bold', pad=20)
ax2.set_ylim(0, max(counts) * 1.2)
ax2.grid(axis='y', alpha=0.3, linestyle='--')

# Add interpretation box
interpretation = (
    "INTERPRETATION:\n"
    f"• {n_ihs_validated}/{n_ihs_total} iHS signals VALIDATED by XP-EHH\n"
    "  → Strong evidence for positive selection\n"
    f"• {n_ihs_only}/{n_ihs_total} iHS signals NOT validated\n"
    "  → Incomplete sweeps or false positives"
)
ax2.text(0.98, 0.98, interpretation,
         transform=ax2.transAxes, ha='right', va='top',
         fontsize=10, family='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
outfile = output_dir / "iHS_validation_summary.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# VISUALIZATION 2: Genome-wide view with clear labels
# ============================================================================
print("Creating Visualization 2: Genome-wide validation map...")

fig, ax = plt.subplots(figsize=(18, 6))

y_pos = 0
region_data = []

# Plot validated iHS regions (HIGH CONFIDENCE)
for key, region in ihs_regions_dict.items():
    ax.barh(y_pos, region['end'] - region['start'], left=region['start'],
            height=0.6, color='#2ecc71', alpha=0.7, edgecolor='black', linewidth=0.5)
    region_data.append(('Validated', region['start'], region['end'], region['score']))
    y_pos += 1

y_pos += 1  # Gap

# Plot iHS-only regions (NOT VALIDATED)
for _, region in ihs_only.iterrows():
    ax.barh(y_pos, region['end'] - region['start'], left=region['start'],
            height=0.6, color='#e74c3c', alpha=0.7, edgecolor='black', linewidth=0.5)
    region_data.append(('Not Validated', region['start'], region['end'], region['max_ihs']))
    y_pos += 1

# Formatting
ax.set_xlabel('Position on Chromosome X (bp)', fontsize=13, fontweight='bold')
ax.set_ylabel('iHS Selection Regions', fontsize=13, fontweight='bold')
ax.set_title('Validation of October iHS Findings Using XP-EHH Cross-Population Comparison\n'
             'iHS: Selection signals from EUR samples | XP-EHH: EUR vs EAS comparison',
             fontsize=14, fontweight='bold', pad=20)

ax.set_xlim(0, 160000000)
ax.ticklabel_format(style='plain', axis='x')
ax.set_yticks([])

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#2ecc71', edgecolor='black', alpha=0.7,
                   label=f'HIGH CONFIDENCE: Validated by XP-EHH (n={n_ihs_validated})'),
    mpatches.Patch(facecolor='#e74c3c', edgecolor='black', alpha=0.7,
                   label=f'INCOMPLETE/FALSE POSITIVE: iHS only (n={n_ihs_only})')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)

ax.grid(axis='x', alpha=0.3, linestyle='--')

plt.tight_layout()
outfile = output_dir / "iHS_validation_genome_wide.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# VISUALIZATION 3: Detailed regional view of top validated regions
# ============================================================================
print("Creating Visualization 3: Top validated regions detail...")

# Get top 10 validated regions by iHS score
top_validated = high_conf.nlargest(10, 'ihs_max')

fig, axes = plt.subplots(5, 2, figsize=(18, 22))
axes = axes.flatten()

for idx, (_, region) in enumerate(top_validated.iterrows()):
    if idx >= 10:
        break

    ax = axes[idx]

    ihs_range = region['ihs_region'].split('-')
    xpehh_range = region['xpehh_region'].split('-')

    ihs_start, ihs_end = int(ihs_range[0]), int(ihs_range[1])
    xpehh_start, xpehh_end = int(xpehh_range[0]), int(xpehh_range[1])

    # Plot regions
    ax.barh(1.5, ihs_end - ihs_start, left=ihs_start, height=0.4,
            color='#3498db', alpha=0.8, edgecolor='black', linewidth=2,
            label=f'iHS region (|iHS|={abs(region["ihs_max"]):.2f})')

    ax.barh(0.5, xpehh_end - xpehh_start, left=xpehh_start, height=0.4,
            color='#e67e22', alpha=0.8, edgecolor='black', linewidth=2,
            label=f'XP-EHH EUR vs EAS (|XP-EHH|={abs(region["xpehh_max"]):.2f})')

    # Highlight overlap
    overlap_start = region['overlap_start']
    overlap_end = region['overlap_end']
    ax.axvspan(overlap_start, overlap_end, alpha=0.2, color='green', zorder=0,
               label=f'Overlap ({region["overlap_kb"]:.1f} kb)')

    # Peak markers
    ax.plot(region['ihs_peak_pos'], 1.5, 'r*', markersize=18,
            markeredgecolor='black', markeredgewidth=1.5, zorder=10)
    ax.plot(region['xpehh_peak_pos'], 0.5, 'r*', markersize=18,
            markeredgecolor='black', markeredgewidth=1.5, zorder=10)

    # Labels
    ax.text(region['ihs_peak_pos'], 1.5, f'  iHS peak',
            va='center', fontsize=8, fontweight='bold')
    ax.text(region['xpehh_peak_pos'], 0.5, f'  XP-EHH peak',
            va='center', fontsize=8, fontweight='bold')

    ax.set_ylim(0, 2)
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(['XP-EHH\n(EUR vs EAS)', 'iHS\n(EUR)'], fontsize=9)
    ax.set_xlabel('Position (bp)', fontsize=10)
    ax.set_title(f'Region {idx+1}: chrX:{overlap_start:,}-{overlap_end:,}\n'
                 f'HIGH CONFIDENCE - Both methods agree',
                 fontsize=11, fontweight='bold', color='darkgreen')
    ax.ticklabel_format(style='plain', axis='x')
    ax.legend(loc='upper right', fontsize=7)
    ax.grid(True, alpha=0.2, axis='x')

plt.suptitle('Top 10 Validated iHS Regions: Detailed View\n'
             'Green overlap = High confidence selection signal',
             fontsize=16, fontweight='bold', y=0.998)
plt.tight_layout()
outfile = output_dir / "top_validated_regions_detail.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# VISUALIZATION 4: Signal strength comparison
# ============================================================================
print("Creating Visualization 4: Signal strength correlation...")

fig, ax = plt.subplots(figsize=(11, 10))

ihs_scores = high_conf['ihs_max'].abs()
xpehh_scores = high_conf['xpehh_max'].abs()

scatter = ax.scatter(ihs_scores, xpehh_scores,
                    s=high_conf['overlap_kb']*0.5,
                    c=high_conf['overlap_kb'],
                    cmap='viridis', alpha=0.6,
                    edgecolors='black', linewidth=1)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Overlap Size (kb)', fontsize=12, fontweight='bold')

# Reference lines
max_val = max(ihs_scores.max(), xpehh_scores.max())
ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.4, linewidth=2,
        label='Equal strength line')
ax.axvline(2.0, color='blue', linestyle=':', linewidth=2, alpha=0.6,
          label='iHS threshold (|iHS| ≥ 2.0)')
ax.axhline(3.0, color='orange', linestyle=':', linewidth=2, alpha=0.6,
          label='XP-EHH threshold (|XP-EHH| ≥ 3.0)')

ax.set_xlabel('|iHS Score| (EUR samples, October analysis)',
              fontsize=13, fontweight='bold')
ax.set_ylabel('|XP-EHH Score| (EUR vs EAS)',
              fontsize=13, fontweight='bold')
ax.set_title('Correlation Between iHS and XP-EHH Signals\n'
             'High Confidence Regions Only',
             fontsize=14, fontweight='bold', pad=20)
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)

# Correlation
correlation = np.corrcoef(ihs_scores, xpehh_scores)[0, 1]
ax.text(0.95, 0.05,
        f'Pearson correlation: {correlation:.3f}\n'
        f'n = {len(high_conf)} overlap instances',
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

plt.tight_layout()
outfile = output_dir / "iHS_XPEHH_correlation.png"
plt.savefig(outfile, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {outfile.name}")

# ============================================================================
# Save summary statistics
# ============================================================================
print("\nCreating summary statistics file...")

summary = {
    'analysis_date': '2025-12-04',
    'ihs_analysis_date': 'October 2024',
    'total_ihs_regions': n_ihs_total,
    'ihs_regions_validated_by_xpehh': n_ihs_validated,
    'ihs_regions_not_validated': n_ihs_only,
    'validation_rate_percent': (n_ihs_validated / n_ihs_total * 100) if n_ihs_total > 0 else 0,
    'total_overlap_instances': len(high_conf),
    'xpehh_comparison': 'EUR vs EAS',
    'ihs_threshold': 2.0,
    'xpehh_threshold': 3.0,
    'mean_overlap_size_kb': high_conf['overlap_kb'].mean(),
    'median_overlap_size_kb': high_conf['overlap_kb'].median(),
    'mean_ihs_score_validated': ihs_scores.mean(),
    'mean_xpehh_score': xpehh_scores.mean(),
    'correlation_ihs_xpehh': correlation,
    'eur_specific_variants': len(eur_specific)
}

summary_df = pd.DataFrame([summary]).T
summary_df.columns = ['value']
summary_df.to_csv(output_dir / "iHS_validation_summary_statistics.tsv", sep='\t')
print(f"  Saved: iHS_validation_summary_statistics.tsv")

# Save validated regions list
print("\nSaving validated iHS regions...")
validated_regions = []
for key, region in ihs_regions_dict.items():
    validated_regions.append({
        'chr': 'X',
        'start': region['start'],
        'end': region['end'],
        'length_kb': (region['end'] - region['start']) / 1000,
        'ihs_max': region['score'],
        'ihs_peak': region['peak'],
        'validation_status': 'HIGH_CONFIDENCE',
        'evidence': 'iHS + XP-EHH'
    })

for _, region in ihs_only.iterrows():
    validated_regions.append({
        'chr': 'X',
        'start': region['start'],
        'end': region['end'],
        'length_kb': region['length_kb'],
        'ihs_max': region['max_ihs'],
        'ihs_peak': region['peak_position'],
        'validation_status': 'INCOMPLETE_OR_FALSE_POSITIVE',
        'evidence': 'iHS only'
    })

validated_df = pd.DataFrame(validated_regions)
validated_df = validated_df.sort_values('ihs_max', key=abs, ascending=False)
validated_df.to_csv(output_dir / "iHS_regions_with_validation_status.tsv",
                    sep='\t', index=False)
print(f"  Saved: iHS_regions_with_validation_status.tsv")

print("\n" + "=" * 80)
print("VALIDATION COMPLETE!")
print("=" * 80)
print(f"\nGenerated files:")
print(f"  1. iHS_validation_summary.png")
print(f"  2. iHS_validation_genome_wide.png")
print(f"  3. top_validated_regions_detail.png")
print(f"  4. iHS_XPEHH_correlation.png")
print(f"  5. iHS_validation_summary_statistics.tsv")
print(f"  6. iHS_regions_with_validation_status.tsv")
print(f"\n{'='*80}")
print("KEY FINDINGS")
print("=" * 80)
print(f"\nOctober iHS Analysis Results:")
print(f"  • Total iHS regions detected: {n_ihs_total}")
print(f"\nValidation with XP-EHH (EUR vs EAS):")
print(f"  • HIGH CONFIDENCE regions: {n_ihs_validated} ({n_ihs_validated/n_ihs_total*100:.1f}%)")
print(f"    → Both iHS and XP-EHH agree")
print(f"    → Strong evidence for positive selection")
print(f"\n  • INCOMPLETE/FALSE POSITIVE regions: {n_ihs_only} ({n_ihs_only/n_ihs_total*100:.1f}%)")
print(f"    → iHS signal without XP-EHH support")
print(f"    → May be incomplete sweeps or false positives")
print(f"\nSignal Properties:")
print(f"  • Correlation between iHS and XP-EHH: {correlation:.3f}")
print(f"  • Mean overlap size: {high_conf['overlap_kb'].mean():.1f} kb")
print(f"  • EUR-specific signals (positive in both EUR comparisons): {len(eur_specific)}")
print("\n" + "=" * 80)
