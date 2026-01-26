#!/usr/bin/env python3
"""
Generate candidate identification summary for X chromosome iHS analysis.
Shows total variants, candidates at Q99 threshold, and candidate rate per population.

Needed this for project report to summarize candidate counts using population-specific thresholds. 
It reads the Q99 thresholds from results/ihs_candidates_summary_all_populations.tsv
For each population, it count show many variants have |Std.iHS|≥Q99
Calculates candidate rate 
Outputs: results/candidate_identification_summary.tsv
"""

import sys
import os

def count_candidates_exceeding_threshold(ihs_file, q99_threshold):
    """Count variants where |Std iHS| >= Q99 threshold."""
    count = 0
    total = 0

    with open(ihs_file, 'r') as f:
        # Skip header
        header = f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            try:
                std_ihs = float(parts[4])  # Std iHS is the 5th column
                total += 1

                # Check if absolute value exceeds threshold
                if abs(std_ihs) >= q99_threshold:
                    count += 1
            except ValueError:
                continue

    return count, total

def main():
    # Input file path
    summary_file = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs_candidates_summary_all_populations.tsv"
    results_base = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results"

    # Read the summary file
    populations = []
    q99_thresholds = {}

    with open(summary_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue

            pop = parts[0]
            chrom = parts[1]
            q99 = float(parts[6])

            if chrom == 'X':
                populations.append(pop)
                q99_thresholds[pop] = q99

    # Sort populations
    populations.sort()

    # Generate summary
    results = []

    print("=" * 100)
    print("X Chromosome iHS Candidate Identification Summary")
    print("=" * 100)
    print()
    print(f"{'Population':<12} {'Total_Variants':>15} {'Q99_Threshold':>15} {'N_Candidates':>15} {'Candidate_Rate_%':>18}")
    print("-" * 100)

    for pop in populations:
        ihs_file = os.path.join(results_base, f"ihs_{pop}", f"{pop}.chrX.ihs.tsv")

        if not os.path.exists(ihs_file):
            print(f"{pop:<12} {'FILE NOT FOUND':>15}", file=sys.stderr)
            continue

        q99 = q99_thresholds[pop]
        n_candidates, total_variants = count_candidates_exceeding_threshold(ihs_file, q99)
        candidate_rate = (n_candidates / total_variants * 100) if total_variants > 0 else 0

        results.append({
            'population': pop,
            'total_variants': total_variants,
            'q99': q99,
            'n_candidates': n_candidates,
            'candidate_rate': candidate_rate
        })

        print(f"{pop:<12} {total_variants:>15,} {q99:>15.4f} {n_candidates:>15,} {candidate_rate:>17.4f}%")

    print("-" * 100)

    # Summary statistics
    if results:
        mean_rate = sum(r['candidate_rate'] for r in results) / len(results)
        sorted_rates = sorted(r['candidate_rate'] for r in results)
        median_rate = sorted_rates[len(sorted_rates)//2]
        min_rate = min(r['candidate_rate'] for r in results)
        max_rate = max(r['candidate_rate'] for r in results)

        print()
        print(f"Total populations analyzed: {len(results)}")
        print(f"Mean candidate rate: {mean_rate:.4f}%")
        print(f"Median candidate rate: {median_rate:.4f}%")
        print(f"Range: {min_rate:.4f}% - {max_rate:.4f}%")
        print("=" * 100)

    # Save to file
    output_file = "/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/candidate_identification_summary.tsv"
    with open(output_file, 'w') as f:
        f.write("Population\tTotal_Variants\tQ99_Threshold\tN_Candidates_Q99\tCandidate_Rate_%\n")
        for r in results:
            f.write(f"{r['population']}\t{r['total_variants']}\t{r['q99']:.6f}\t{r['n_candidates']}\t{r['candidate_rate']:.6f}\n")

    print(f"\nSummary saved to: {output_file}")

if __name__ == "__main__":
    main()
