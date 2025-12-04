#!/usr/bin/env python3
"""
Comprehensive comparison of iHS and XP-EHH results
Following the Sabeti et al. framework for validating selection signals
"""

import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results")
PLOTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots")
OUTPUT_DIR = PLOTS_DIR / "ihs_xpehh_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Thresholds based on Sabeti et al.
IHS_THRESHOLD = 2.0  # |std iHS| > 2
XPEHH_STRONG = 3.0   # |std XP-EHH| > 3 for strong signal
XPEHH_SUGGESTIVE = 2.0  # |std XP-EHH| > 2 for suggestive signal

def parse_position(variant_id):
    """Extract position from variant ID (format: X:position:ref:alt)"""
    try:
        parts = variant_id.split(':')
        return int(parts[1])
    except:
        return None

def load_ihs_data(filepath):
    """Load iHS results"""
    print(f"Loading iHS data from {filepath.name}...")

    positions = []
    ihs_scores = []
    variant_ids = []

    with open(filepath, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variant_id = parts[0]
                pos = parse_position(variant_id)

                if pos is None:
                    continue

                try:
                    std_ihs = float(parts[4])
                    if np.isfinite(std_ihs):
                        positions.append(pos)
                        ihs_scores.append(std_ihs)
                        variant_ids.append(variant_id)
                except:
                    continue

    print(f"  Loaded {len(positions):,} variants with valid iHS scores")
    print(f"  Position range: {min(positions):,} - {max(positions):,}")
    print(f"  iHS range: {min(ihs_scores):.2f} - {max(ihs_scores):.2f}")

    return positions, ihs_scores, variant_ids

def load_xpehh_data(filepath):
    """Load XP-EHH results"""
    print(f"Loading XP-EHH data from {filepath.name}...")

    positions = []
    xpehh_scores = []
    variant_ids = []

    with open(filepath, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                variant_id = parts[1]
                pos = parse_position(variant_id)

                if pos is None:
                    continue

                try:
                    std_xpehh = float(parts[7])
                    if np.isfinite(std_xpehh):
                        positions.append(pos)
                        xpehh_scores.append(std_xpehh)
                        variant_ids.append(variant_id)
                except:
                    continue

    print(f"  Loaded {len(positions):,} variants with valid XP-EHH scores")
    print(f"  Position range: {min(positions):,} - {max(positions):,}")
    print(f"  XP-EHH range: {min(xpehh_scores):.2f} - {max(xpehh_scores):.2f}")

    return positions, xpehh_scores, variant_ids

def find_peaks(positions, scores, threshold, window_size=500000):
    """
    Identify genomic regions with clusters of high scores
    Returns list of (start, end, n_variants, max_score, peak_position)
    """
    # Find significant variants
    sig_indices = [i for i, score in enumerate(scores) if abs(score) > threshold]

    if len(sig_indices) == 0:
        return []

    # Sort by position
    sig_data = [(positions[i], scores[i]) for i in sig_indices]
    sig_data.sort(key=lambda x: x[0])

    # Group nearby variants into regions
    regions = []
    current_region = {
        'start': sig_data[0][0],
        'end': sig_data[0][0],
        'variants': [sig_data[0]],
        'max_score': sig_data[0][1],
        'peak_pos': sig_data[0][0]
    }

    for pos, score in sig_data[1:]:
        if pos - current_region['end'] <= window_size:
            # Extend current region
            current_region['end'] = pos
            current_region['variants'].append((pos, score))
            if abs(score) > abs(current_region['max_score']):
                current_region['max_score'] = score
                current_region['peak_pos'] = pos
        else:
            # Save current region and start new one
            regions.append(current_region)
            current_region = {
                'start': pos,
                'end': pos,
                'variants': [(pos, score)],
                'max_score': score,
                'peak_pos': pos
            }

    # Add last region
    regions.append(current_region)

    # Convert to simple tuples
    result = [(r['start'], r['end'], len(r['variants']), r['max_score'], r['peak_pos'])
              for r in regions]

    # Sort by absolute max score
    result.sort(key=lambda x: abs(x[3]), reverse=True)

    return result

def regions_overlap(r1_start, r1_end, r2_start, r2_end):
    """Check if two regions overlap"""
    return not (r1_end < r2_start or r2_end < r1_start)

def find_overlapping_regions(ihs_peaks, xpehh_peaks):
    """Find regions that show signals in both iHS and XP-EHH"""
    overlaps = []

    for ihs_region in ihs_peaks:
        ihs_start, ihs_end, ihs_n, ihs_max, ihs_peak = ihs_region

        for xpehh_region in xpehh_peaks:
            xpehh_start, xpehh_end, xpehh_n, xpehh_max, xpehh_peak = xpehh_region

            if regions_overlap(ihs_start, ihs_end, xpehh_start, xpehh_end):
                # Calculate overlap
                overlap_start = max(ihs_start, xpehh_start)
                overlap_end = min(ihs_end, xpehh_end)
                overlap_size = overlap_end - overlap_start

                overlaps.append({
                    'chr': 'X',
                    'overlap_start': overlap_start,
                    'overlap_end': overlap_end,
                    'overlap_kb': overlap_size / 1000,
                    'ihs_start': ihs_start,
                    'ihs_end': ihs_end,
                    'ihs_n': ihs_n,
                    'ihs_max': ihs_max,
                    'ihs_peak': ihs_peak,
                    'xpehh_start': xpehh_start,
                    'xpehh_end': xpehh_end,
                    'xpehh_n': xpehh_n,
                    'xpehh_max': xpehh_max,
                    'xpehh_peak': xpehh_peak
                })

    # Sort by iHS max score
    overlaps.sort(key=lambda x: abs(x['ihs_max']), reverse=True)

    return overlaps

def find_eur_specific_signals(eur_eas_data, eur_afr_data, threshold=3.0):
    """
    Find positions where BOTH EUR vs EAS and EUR vs AFR show strong positive XP-EHH
    This indicates EUR-specific selection (not just drift)
    """
    print(f"\nFinding EUR-specific selection signals (threshold: {threshold})...")

    eur_eas_pos, eur_eas_scores, _ = eur_eas_data
    eur_afr_pos, eur_afr_scores, _ = eur_afr_data

    # Create position-to-score mappings
    eur_eas_dict = {pos: score for pos, score in zip(eur_eas_pos, eur_eas_scores)}
    eur_afr_dict = {pos: score for pos, score in zip(eur_afr_pos, eur_afr_scores)}

    # Find positions present in both datasets
    common_positions = set(eur_eas_dict.keys()) & set(eur_afr_dict.keys())

    # Find EUR-specific signals (positive in both comparisons)
    eur_specific = []
    for pos in common_positions:
        eas_score = eur_eas_dict[pos]
        afr_score = eur_afr_dict[pos]

        if eas_score > threshold and afr_score > threshold:
            eur_specific.append({
                'position': pos,
                'eur_vs_eas': eas_score,
                'eur_vs_afr': afr_score,
                'mean_score': (eas_score + afr_score) / 2
            })

    # Sort by mean score
    eur_specific.sort(key=lambda x: x['mean_score'], reverse=True)

    print(f"  Found {len(eur_specific):,} EUR-specific variants")
    print(f"  (positive XP-EHH > {threshold} in BOTH EUR vs EAS and EUR vs AFR)")

    return eur_specific

def categorize_signals(ihs_peaks, xpehh_peaks, overlaps):
    """
    Categorize signals according to the framework:
    - High confidence: Both iHS and XP-EHH agree
    - XP-EHH only: Possible near-fixation sweep (iHS loses power)
    - iHS only: Possible incomplete sweep
    """

    # Get all iHS and XP-EHH peak positions
    ihs_regions = set()
    for start, end, _, _, _ in ihs_peaks:
        ihs_regions.add((start, end))

    xpehh_regions = set()
    for start, end, _, _, _ in xpehh_peaks:
        xpehh_regions.add((start, end))

    # Categorize
    high_conf = overlaps  # Regions in both

    # iHS only regions
    ihs_only = []
    for ihs_region in ihs_peaks:
        start, end, n, max_score, peak = ihs_region
        found_overlap = False
        for overlap in overlaps:
            if regions_overlap(start, end, overlap['xpehh_start'], overlap['xpehh_end']):
                found_overlap = True
                break
        if not found_overlap:
            ihs_only.append(ihs_region)

    # XP-EHH only regions
    xpehh_only = []
    for xpehh_region in xpehh_peaks:
        start, end, n, max_score, peak = xpehh_region
        found_overlap = False
        for overlap in overlaps:
            if regions_overlap(start, end, overlap['ihs_start'], overlap['ihs_end']):
                found_overlap = True
                break
        if not found_overlap:
            xpehh_only.append(xpehh_region)

    return {
        'high_confidence': high_conf,
        'ihs_only': ihs_only,
        'xpehh_only': xpehh_only
    }

def main():
    print("=" * 80)
    print("iHS vs XP-EHH Comparison Analysis")
    print("Following Sabeti et al. framework for validating selection signals")
    print("=" * 80)

    # Load iHS data
    print("\n" + "=" * 80)
    print("STEP 1: Loading iHS Results")
    print("=" * 80)
    ihs_file = RESULTS_DIR / "ihs" / "ALL.chrX.ihs.tsv"
    ihs_pos, ihs_scores, ihs_ids = load_ihs_data(ihs_file)

    # Load XP-EHH data
    print("\n" + "=" * 80)
    print("STEP 2: Loading XP-EHH Results")
    print("=" * 80)

    xpehh_dir = RESULTS_DIR / "xp_ehh"

    eur_eas_file = xpehh_dir / "chrX_EUR_vs_EAS"
    eur_afr_file = xpehh_dir / "chrX_EUR_vs_AFR"
    afr_eas_file = xpehh_dir / "chrX_AFR_vs_EAS"

    eur_eas_data = load_xpehh_data(eur_eas_file)
    eur_afr_data = load_xpehh_data(eur_afr_file)
    afr_eas_data = load_xpehh_data(afr_eas_file)

    # Find EUR-specific signals
    print("\n" + "=" * 80)
    print("STEP 3: Identifying EUR-Specific Selection Signals")
    print("=" * 80)
    eur_specific = find_eur_specific_signals(eur_eas_data, eur_afr_data, threshold=XPEHH_STRONG)

    # Save EUR-specific signals
    if len(eur_specific) > 0:
        eur_file = OUTPUT_DIR / "EUR_specific_signals.tsv"
        with open(eur_file, 'w') as f:
            f.write("chr\tposition\teur_vs_eas_xpehh\teur_vs_afr_xpehh\tmean_xpehh\n")
            for sig in eur_specific[:100]:  # Top 100
                f.write(f"X\t{sig['position']}\t{sig['eur_vs_eas']:.4f}\t"
                       f"{sig['eur_vs_afr']:.4f}\t{sig['mean_score']:.4f}\n")
        print(f"\nSaved top EUR-specific signals to: {eur_file}")

    # Find peaks in each dataset
    print("\n" + "=" * 80)
    print("STEP 4: Identifying Peak Regions")
    print("=" * 80)

    print("\niHS peaks (|std iHS| > 2.0):")
    ihs_peaks = find_peaks(ihs_pos, ihs_scores, IHS_THRESHOLD)
    print(f"  Found {len(ihs_peaks)} iHS peak regions")

    print("\nXP-EHH peaks (EUR vs EAS, |std XP-EHH| > 3.0):")
    eur_eas_peaks = find_peaks(eur_eas_data[0], eur_eas_data[1], XPEHH_STRONG)
    print(f"  Found {len(eur_eas_peaks)} EUR vs EAS peak regions")

    print("\nXP-EHH peaks (EUR vs AFR, |std XP-EHH| > 3.0):")
    eur_afr_peaks = find_peaks(eur_afr_data[0], eur_afr_data[1], XPEHH_STRONG)
    print(f"  Found {len(eur_afr_peaks)} EUR vs AFR peak regions")

    print("\nXP-EHH peaks (AFR vs EAS, |std XP-EHH| > 3.0):")
    afr_eas_peaks = find_peaks(afr_eas_data[0], afr_eas_data[1], XPEHH_STRONG)
    print(f"  Found {len(afr_eas_peaks)} AFR vs EAS peak regions")

    # Compare iHS with each XP-EHH comparison
    print("\n" + "=" * 80)
    print("STEP 5: Finding Overlapping Regions (iHS + XP-EHH)")
    print("=" * 80)

    print("\niHS vs EUR vs EAS XP-EHH:")
    overlaps_eur_eas = find_overlapping_regions(ihs_peaks, eur_eas_peaks)
    print(f"  Found {len(overlaps_eur_eas)} overlapping regions")

    print("\niHS vs EUR vs AFR XP-EHH:")
    overlaps_eur_afr = find_overlapping_regions(ihs_peaks, eur_afr_peaks)
    print(f"  Found {len(overlaps_eur_afr)} overlapping regions")

    print("\niHS vs AFR vs EAS XP-EHH:")
    overlaps_afr_eas = find_overlapping_regions(ihs_peaks, afr_eas_peaks)
    print(f"  Found {len(overlaps_afr_eas)} overlapping regions")

    # Categorize signals for EUR vs EAS
    print("\n" + "=" * 80)
    print("STEP 6: Categorizing Signals (EUR vs EAS comparison)")
    print("=" * 80)

    categories = categorize_signals(ihs_peaks, eur_eas_peaks, overlaps_eur_eas)

    print(f"\nHIGH CONFIDENCE (both iHS and XP-EHH): {len(categories['high_confidence'])} regions")
    print(f"  → Strong evidence for selection")

    print(f"\nXP-EHH ONLY: {len(categories['xpehh_only'])} regions")
    print(f"  → Possible near-fixation sweeps (iHS loses power)")

    print(f"\niHS ONLY: {len(categories['ihs_only'])} regions")
    print(f"  → Possible incomplete sweeps or iHS more sensitive")

    # Save results
    print("\n" + "=" * 80)
    print("STEP 7: Saving Results")
    print("=" * 80)

    # High confidence regions (EUR vs EAS)
    if len(categories['high_confidence']) > 0:
        outfile = OUTPUT_DIR / "high_confidence_regions_EUR_vs_EAS.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\toverlap_start\toverlap_end\toverlap_kb\t"
                   "ihs_max\tihs_peak_pos\txpehh_max\txpehh_peak_pos\t"
                   "ihs_region\txpehh_region\n")
            for region in categories['high_confidence']:
                f.write(f"{region['chr']}\t{region['overlap_start']}\t{region['overlap_end']}\t"
                       f"{region['overlap_kb']:.2f}\t{region['ihs_max']:.4f}\t{region['ihs_peak']}\t"
                       f"{region['xpehh_max']:.4f}\t{region['xpehh_peak']}\t"
                       f"{region['ihs_start']}-{region['ihs_end']}\t"
                       f"{region['xpehh_start']}-{region['xpehh_end']}\n")
        print(f"Saved high confidence regions: {outfile}")

    # iHS only regions
    if len(categories['ihs_only']) > 0:
        outfile = OUTPUT_DIR / "ihs_only_regions.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tstart\tend\tn_variants\tmax_ihs\tpeak_position\tlength_kb\n")
            for region in categories['ihs_only']:
                start, end, n, max_score, peak = region
                f.write(f"X\t{start}\t{end}\t{n}\t{max_score:.4f}\t{peak}\t"
                       f"{(end-start)/1000:.2f}\n")
        print(f"Saved iHS-only regions: {outfile}")

    # XP-EHH only regions
    if len(categories['xpehh_only']) > 0:
        outfile = OUTPUT_DIR / "xpehh_only_regions_EUR_vs_EAS.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tstart\tend\tn_variants\tmax_xpehh\tpeak_position\tlength_kb\n")
            for region in categories['xpehh_only']:
                start, end, n, max_score, peak = region
                f.write(f"X\t{start}\t{end}\t{n}\t{max_score:.4f}\t{peak}\t"
                       f"{(end-start)/1000:.2f}\n")
        print(f"Saved XP-EHH-only regions: {outfile}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)

    print("\nSignal Distribution:")
    print(f"  Total iHS peaks (|iHS| > 2): {len(ihs_peaks)}")
    print(f"  Total EUR vs EAS XP-EHH peaks (|XP-EHH| > 3): {len(eur_eas_peaks)}")
    print(f"  Total EUR vs AFR XP-EHH peaks (|XP-EHH| > 3): {len(eur_afr_peaks)}")
    print(f"  Total AFR vs EAS XP-EHH peaks (|XP-EHH| > 3): {len(afr_eas_peaks)}")

    print("\nValidation (iHS + EUR vs EAS XP-EHH):")
    print(f"  High confidence (both methods): {len(categories['high_confidence'])}")
    print(f"  iHS only (incomplete sweep?): {len(categories['ihs_only'])}")
    print(f"  XP-EHH only (near fixation?): {len(categories['xpehh_only'])}")

    if len(ihs_peaks) > 0:
        validation_rate = len(categories['high_confidence']) / len(ihs_peaks) * 100
        print(f"\n  iHS validation rate: {validation_rate:.1f}% of iHS peaks confirmed by XP-EHH")

    print("\nEUR-specific selection:")
    print(f"  Variants with strong signal in BOTH EUR comparisons: {len(eur_specific)}")

    # Print top high-confidence regions
    if len(categories['high_confidence']) > 0:
        print("\n" + "=" * 80)
        print("TOP 20 HIGH CONFIDENCE REGIONS (iHS + XP-EHH agree)")
        print("=" * 80)
        print(f"{'Rank':<5} {'Position Range (Mb)':<25} {'iHS':<10} {'XP-EHH':<10} {'Size (kb)':<12}")
        print("-" * 80)

        for i, region in enumerate(categories['high_confidence'][:20], 1):
            start_mb = region['overlap_start'] / 1e6
            end_mb = region['overlap_end'] / 1e6
            pos_range = f"{start_mb:.2f} - {end_mb:.2f}"
            print(f"{i:<5} {pos_range:<25} {region['ihs_max']:>9.2f} {region['xpehh_max']:>9.2f} "
                  f"{region['overlap_kb']:>11.1f}")

    print("\n" + "=" * 80)
    print("Analysis Complete!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
