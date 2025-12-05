#!/usr/bin/env python3
"""
Comprehensive EUR-specific selection analysis
Combines EUR iHS with XP-EHH (EUR vs EAS and EUR vs AFR)

Following Sabeti et al. framework:
- High confidence: EUR iHS + both XP-EHH comparisons agree
- Medium confidence: EUR iHS + one XP-EHH comparison
- EUR iHS only: Incomplete sweep
- XP-EHH only: Near-fixation sweep (iHS loses power)
"""

import sys
import math
from pathlib import Path
from collections import defaultdict

# Directories
RESULTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results")
PLOTS_DIR = Path("/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/plots")
OUTPUT_DIR = PLOTS_DIR / "EUR_selection_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Thresholds
IHS_THRESHOLD = 2.0
XPEHH_STRONG = 3.0
XPEHH_SUGGESTIVE = 2.0

def parse_position(variant_id):
    """Extract position from variant ID"""
    try:
        parts = variant_id.split(':')
        return int(parts[1])
    except:
        return None

def load_eur_ihs(filepath):
    """Load EUR-specific iHS results"""
    print(f"Loading EUR iHS from {filepath.name}...")

    data = {}
    with open(filepath, 'r') as f:
        header = f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variant_id = parts[0]
                pos = parse_position(variant_id)

                if pos is None:
                    continue

                try:
                    std_ihs = float(parts[4])
                    if not math.isnan(std_ihs) and not math.isinf(std_ihs):
                        data[pos] = {
                            'variant_id': variant_id,
                            'std_ihs': std_ihs
                        }
                except:
                    continue

    print(f"  Loaded {len(data):,} variants")
    print(f"  iHS range: {min(d['std_ihs'] for d in data.values()):.2f} - "
          f"{max(d['std_ihs'] for d in data.values()):.2f}")

    return data

def load_xpehh(filepath):
    """Load XP-EHH results"""
    print(f"Loading XP-EHH from {filepath.name}...")

    data = {}
    with open(filepath, 'r') as f:
        header = f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                variant_id = parts[1]
                pos = parse_position(variant_id)

                if pos is None:
                    continue

                try:
                    std_xpehh = float(parts[7])
                    if not math.isnan(std_xpehh) and not math.isinf(std_xpehh):
                        data[pos] = {
                            'variant_id': variant_id,
                            'std_xpehh': std_xpehh
                        }
                except:
                    continue

    print(f"  Loaded {len(data):,} variants")
    print(f"  XP-EHH range: {min(d['std_xpehh'] for d in data.values()):.2f} - "
          f"{max(d['std_xpehh'] for d in data.values()):.2f}")

    return data

def find_eur_specific_positions(eur_ihs, eur_eas, eur_afr):
    """
    Find positions with EUR-specific selection signals

    Categories:
    1. HIGH CONFIDENCE: EUR iHS > 2 AND EUR vs EAS > 3 AND EUR vs AFR > 3
    2. MEDIUM CONFIDENCE: EUR iHS > 2 AND (EUR vs EAS > 3 OR EUR vs AFR > 3)
    3. EUR iHS ONLY: EUR iHS > 2 but XP-EHH weak
    4. XP-EHH ONLY: Both XP-EHH > 3 but EUR iHS weak (near fixation)
    """

    # Get common positions across all datasets
    all_positions = set(eur_ihs.keys()) | set(eur_eas.keys()) | set(eur_afr.keys())

    high_confidence = []
    medium_confidence = []
    ihs_only = []
    xpehh_only = []

    for pos in all_positions:
        # Get scores (default to 0 if missing)
        ihs_score = eur_ihs.get(pos, {}).get('std_ihs', 0)
        eas_score = eur_eas.get(pos, {}).get('std_xpehh', 0)
        afr_score = eur_afr.get(pos, {}).get('std_xpehh', 0)

        # Check for EUR-specific signals (positive XP-EHH = EUR selection)
        strong_ihs = abs(ihs_score) > IHS_THRESHOLD
        strong_eas = eas_score > XPEHH_STRONG
        strong_afr = afr_score > XPEHH_STRONG

        # Categorize
        if strong_ihs and strong_eas and strong_afr:
            # All three methods agree - highest confidence
            high_confidence.append({
                'position': pos,
                'eur_ihs': ihs_score,
                'eur_vs_eas': eas_score,
                'eur_vs_afr': afr_score,
                'combined_score': abs(ihs_score) + eas_score + afr_score
            })

        elif strong_ihs and (strong_eas or strong_afr):
            # EUR iHS + one XP-EHH - medium confidence
            medium_confidence.append({
                'position': pos,
                'eur_ihs': ihs_score,
                'eur_vs_eas': eas_score,
                'eur_vs_afr': afr_score,
                'combined_score': abs(ihs_score) + max(eas_score, afr_score)
            })

        elif strong_ihs and not (strong_eas or strong_afr):
            # EUR iHS only - incomplete sweep
            ihs_only.append({
                'position': pos,
                'eur_ihs': ihs_score,
                'eur_vs_eas': eas_score,
                'eur_vs_afr': afr_score
            })

        elif (strong_eas and strong_afr) and not strong_ihs:
            # Both XP-EHH strong but iHS weak - near fixation
            xpehh_only.append({
                'position': pos,
                'eur_ihs': ihs_score,
                'eur_vs_eas': eas_score,
                'eur_vs_afr': afr_score,
                'combined_score': eas_score + afr_score
            })

    # Sort by combined score
    high_confidence.sort(key=lambda x: x['combined_score'], reverse=True)
    medium_confidence.sort(key=lambda x: x['combined_score'], reverse=True)
    xpehh_only.sort(key=lambda x: x['combined_score'], reverse=True)
    ihs_only.sort(key=lambda x: abs(x['eur_ihs']), reverse=True)

    return high_confidence, medium_confidence, ihs_only, xpehh_only

def cluster_into_regions(positions_data, window_size=500000):
    """Cluster nearby positions into genomic regions"""
    if len(positions_data) == 0:
        return []

    # Sort by position
    sorted_data = sorted(positions_data, key=lambda x: x['position'])

    regions = []
    current_region = {
        'start': sorted_data[0]['position'],
        'end': sorted_data[0]['position'],
        'positions': [sorted_data[0]],
        'max_combined': sorted_data[0].get('combined_score', abs(sorted_data[0]['eur_ihs']))
    }

    for data in sorted_data[1:]:
        pos = data['position']
        combined = data.get('combined_score', abs(data['eur_ihs']))

        if pos - current_region['end'] <= window_size:
            # Extend region
            current_region['end'] = pos
            current_region['positions'].append(data)
            if combined > current_region['max_combined']:
                current_region['max_combined'] = combined
        else:
            # Save and start new region
            regions.append(current_region)
            current_region = {
                'start': pos,
                'end': pos,
                'positions': [data],
                'max_combined': combined
            }

    # Add last region
    regions.append(current_region)

    # Sort by max combined score
    regions.sort(key=lambda x: x['max_combined'], reverse=True)

    return regions

def save_results(high_conf, medium_conf, ihs_only, xpehh_only, output_dir):
    """Save categorized results to files"""

    print("\nSaving results...")

    # High confidence variants
    if len(high_conf) > 0:
        outfile = output_dir / "HIGH_CONFIDENCE_EUR_selection_variants.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tposition\teur_ihs\teur_vs_eas_xpehh\teur_vs_afr_xpehh\tcombined_score\n")
            for v in high_conf[:1000]:  # Top 1000
                f.write(f"X\t{v['position']}\t{v['eur_ihs']:.4f}\t"
                       f"{v['eur_vs_eas']:.4f}\t{v['eur_vs_afr']:.4f}\t"
                       f"{v['combined_score']:.4f}\n")
        print(f"  Saved: {outfile}")

    # Medium confidence variants
    if len(medium_conf) > 0:
        outfile = output_dir / "MEDIUM_CONFIDENCE_EUR_selection_variants.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tposition\teur_ihs\teur_vs_eas_xpehh\teur_vs_afr_xpehh\tcombined_score\n")
            for v in medium_conf[:1000]:
                f.write(f"X\t{v['position']}\t{v['eur_ihs']:.4f}\t"
                       f"{v['eur_vs_eas']:.4f}\t{v['eur_vs_afr']:.4f}\t"
                       f"{v['combined_score']:.4f}\n")
        print(f"  Saved: {outfile}")

    # iHS only variants
    if len(ihs_only) > 0:
        outfile = output_dir / "iHS_ONLY_variants.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tposition\teur_ihs\teur_vs_eas_xpehh\teur_vs_afr_xpehh\n")
            for v in ihs_only[:1000]:
                f.write(f"X\t{v['position']}\t{v['eur_ihs']:.4f}\t"
                       f"{v['eur_vs_eas']:.4f}\t{v['eur_vs_afr']:.4f}\n")
        print(f"  Saved: {outfile}")

    # XP-EHH only variants
    if len(xpehh_only) > 0:
        outfile = output_dir / "XPEHH_ONLY_variants_near_fixation.tsv"
        with open(outfile, 'w') as f:
            f.write("chr\tposition\teur_ihs\teur_vs_eas_xpehh\teur_vs_afr_xpehh\tcombined_score\n")
            for v in xpehh_only[:1000]:
                f.write(f"X\t{v['position']}\t{v['eur_ihs']:.4f}\t"
                       f"{v['eur_vs_eas']:.4f}\t{v['eur_vs_afr']:.4f}\t"
                       f"{v['combined_score']:.4f}\n")
        print(f"  Saved: {outfile}")

def save_regions(regions, category, output_dir):
    """Save genomic regions to file"""
    if len(regions) == 0:
        return

    outfile = output_dir / f"{category}_EUR_selection_regions.tsv"
    with open(outfile, 'w') as f:
        f.write("chr\tstart\tend\tn_variants\tmax_combined_score\tlength_kb\t"
               "peak_position\tpeak_ihs\tpeak_eur_vs_eas\tpeak_eur_vs_afr\n")

        for region in regions:
            # Find peak variant
            peak = max(region['positions'], key=lambda x: x.get('combined_score', abs(x['eur_ihs'])))

            f.write(f"X\t{region['start']}\t{region['end']}\t{len(region['positions'])}\t"
                   f"{region['max_combined']:.4f}\t{(region['end']-region['start'])/1000:.2f}\t"
                   f"{peak['position']}\t{peak['eur_ihs']:.4f}\t"
                   f"{peak['eur_vs_eas']:.4f}\t{peak['eur_vs_afr']:.4f}\n")

    print(f"  Saved: {outfile}")

def main():
    print("=" * 80)
    print("EUR-SPECIFIC SELECTION ANALYSIS")
    print("Combining EUR iHS + XP-EHH (EUR vs EAS/AFR)")
    print("=" * 80)

    # Load data
    print("\n" + "=" * 80)
    print("STEP 1: Loading Data")
    print("=" * 80)

    eur_ihs = load_eur_ihs(RESULTS_DIR / "ihs_EUR" / "EUR.chrX.ihs.tsv")
    eur_eas = load_xpehh(RESULTS_DIR / "xp_ehh" / "chrX_EUR_vs_EAS")
    eur_afr = load_xpehh(RESULTS_DIR / "xp_ehh" / "chrX_EUR_vs_AFR")

    # Find EUR-specific positions
    print("\n" + "=" * 80)
    print("STEP 2: Identifying EUR-Specific Selection Signals")
    print("=" * 80)

    high_conf, medium_conf, ihs_only, xpehh_only = find_eur_specific_positions(
        eur_ihs, eur_eas, eur_afr
    )

    print(f"\nSignal Categories:")
    print(f"  HIGH CONFIDENCE (all 3 methods): {len(high_conf):,} variants")
    print(f"    → EUR iHS > 2 AND EUR vs EAS > 3 AND EUR vs AFR > 3")
    print(f"  MEDIUM CONFIDENCE (2 methods): {len(medium_conf):,} variants")
    print(f"    → EUR iHS > 2 AND (EUR vs EAS > 3 OR EUR vs AFR > 3)")
    print(f"  iHS ONLY (incomplete sweep): {len(ihs_only):,} variants")
    print(f"    → EUR iHS > 2 but weak XP-EHH")
    print(f"  XP-EHH ONLY (near fixation): {len(xpehh_only):,} variants")
    print(f"    → Both XP-EHH > 3 but weak iHS")

    # Cluster into regions
    print("\n" + "=" * 80)
    print("STEP 3: Clustering Signals into Genomic Regions")
    print("=" * 80)

    print("\nClustering variants within 500kb windows...")
    high_conf_regions = cluster_into_regions(high_conf)
    medium_conf_regions = cluster_into_regions(medium_conf)
    ihs_only_regions = cluster_into_regions(ihs_only)
    xpehh_only_regions = cluster_into_regions(xpehh_only)

    print(f"\nGenomic Regions:")
    print(f"  HIGH CONFIDENCE: {len(high_conf_regions)} regions")
    print(f"  MEDIUM CONFIDENCE: {len(medium_conf_regions)} regions")
    print(f"  iHS ONLY: {len(ihs_only_regions)} regions")
    print(f"  XP-EHH ONLY: {len(xpehh_only_regions)} regions")

    # Save results
    print("\n" + "=" * 80)
    print("STEP 4: Saving Results")
    print("=" * 80)

    save_results(high_conf, medium_conf, ihs_only, xpehh_only, OUTPUT_DIR)

    save_regions(high_conf_regions, "HIGH_CONFIDENCE", OUTPUT_DIR)
    save_regions(medium_conf_regions, "MEDIUM_CONFIDENCE", OUTPUT_DIR)
    save_regions(ihs_only_regions, "iHS_ONLY", OUTPUT_DIR)
    save_regions(xpehh_only_regions, "XPEHH_ONLY", OUTPUT_DIR)

    # Print top high-confidence regions
    if len(high_conf_regions) > 0:
        print("\n" + "=" * 80)
        print("TOP 20 HIGH CONFIDENCE EUR-SPECIFIC SELECTION REGIONS")
        print("All three methods agree: EUR iHS + EUR vs EAS + EUR vs AFR")
        print("=" * 80)
        print(f"{'Rank':<6} {'Region (Mb)':<20} {'Variants':<10} {'Combined':<10} {'Size (kb)':<10}")
        print("-" * 80)

        for i, region in enumerate(high_conf_regions[:20], 1):
            start_mb = region['start'] / 1e6
            end_mb = region['end'] / 1e6
            region_str = f"{start_mb:.2f}-{end_mb:.2f}"
            size_kb = (region['end'] - region['start']) / 1000

            print(f"{i:<6} {region_str:<20} {len(region['positions']):<10} "
                  f"{region['max_combined']:<10.2f} {size_kb:<10.1f}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("INTERPRETATION GUIDE")
    print("=" * 80)

    print("\nHIGH CONFIDENCE regions:")
    print("  → Strong EUR-specific selection (not just drift)")
    print("  → Top priority for functional follow-up")
    print("  → Look for genes related to:")
    print("     - High latitude adaptation (pigmentation, vitamin D)")
    print("     - European-specific pathogens")
    print("     - Metabolic adaptations")

    print("\nMEDIUM CONFIDENCE regions:")
    print("  → Likely EUR selection but less certain")
    print("  → Could be ongoing selection")
    print("  → Worth investigating if overlap with known candidates")

    print("\niHS ONLY regions:")
    print("  → Incomplete selective sweeps")
    print("  → Selection may be ongoing")
    print("  → Allele hasn't reached fixation yet")

    print("\nXP-EHH ONLY regions:")
    print("  → Near-fixation sweeps")
    print("  → Selection completed, allele at high frequency")
    print("  → iHS loses power at near-fixation")

    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("1. Annotate HIGH CONFIDENCE regions with genes")
    print("2. Check for known selection candidates:")
    print("   - SLC24A5, SLC45A2 (pigmentation)")
    print("   - HERC2/OCA2 (eye color)")
    print("   - LCT (lactase persistence)")
    print("   - EDAR (hair/tooth morphology - should be EAS-specific)")
    print("3. Look for enrichment in functional categories")
    print("4. Compare with Sabeti et al. Table 1 findings")

    print("\n" + "=" * 80)
    print("Analysis Complete!")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
