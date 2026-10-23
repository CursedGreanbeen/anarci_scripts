#!/usr/bin/env python3
"""Filter BLAST hits by release date using SAbDab database."""

import csv
from datetime import datetime

# Input files
SABDAB_FILE = "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/sabdab_summary_all.tsv"
BLAST_HITS_FILE = "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/fasta-fv-only-levenshtein/combined_filtered_hits.csv"

# Output files
PRE_2024_OUTPUT = "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/fasta-fv-only-levenshtein/pre_2024_hits.csv"
SUMMARY_OUTPUT = "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/fasta-fv-only-levenshtein/date_filter_summary.txt"

# Cutoff: hits before this date are considered "pre-2024"
CUTOFF_DATE = datetime(2024, 1, 1)

def load_sabdab_dates():
    """Load PDB -> release date mapping from SAbDab TSV."""
    pdb_dates = {}
    with open(SABDAB_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_id = row['pdb'].upper()
            date_str = row['date']
            # Parse MM/DD/YY format
            try:
                release_date = datetime.strptime(date_str, '%m/%d/%y')
                pdb_dates[pdb_id] = release_date
            except ValueError:
                # Try alternative formats if needed
                pass
    return pdb_dates

def extract_pdb_id(hit_id):
    """Extract PDB ID from BLAST hit_id (format: pdb|XXXX|X)."""
    import re
    match = re.search(r'pdb\|([A-Z0-9]{4})\|', hit_id)
    return match.group(1).upper() if match else None

def main():
    print("Loading SAbDab dates...")
    pdb_dates = load_sabdab_dates()
    print(f"Loaded {len(pdb_dates)} PDB entries from SAbDab")

    print(f"\nFiltering BLAST hits (pre-2024 cutoff: {CUTOFF_DATE.date()})...")

    # Counters
    total_hits = 0
    pre_2024_hits = 0
    post_2024_hits = 0
    not_in_sabdab = 0

    # Field names for output
    fieldnames = ['source_file', 'chain_id', 'hit_id', 'hit_description',
                  'pct_identity', 'identity', 'align_length', 'evalue',
                  'bit_score', 'h3_covered', 'h3_identity', 'h3_align_length',
                  'h3_pct_identity', 'query_h3', 'hit_h3', 'release_date', 'pdb_id']

    pre_2024_rows = []
    unique_pre_2024_pdbs = set()
    unique_post_2024_pdbs = set()

    with open(BLAST_HITS_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_hits += 1
            pdb_id = extract_pdb_id(row['hit_id'])

            if pdb_id is None:
                continue

            if pdb_id not in pdb_dates:
                not_in_sabdab += 1
                continue

            release_date = pdb_dates[pdb_id]

            if release_date < CUTOFF_DATE:
                pre_2024_hits += 1
                unique_pre_2024_pdbs.add(pdb_id)
                pre_2024_rows.append({
                    **row,
                    'release_date': release_date.strftime('%Y-%m-%d'),
                    'pdb_id': pdb_id
                })
            else:
                post_2024_hits += 1
                unique_post_2024_pdbs.add(pdb_id)

    # Write pre-2024 hits
    with open(PRE_2024_OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(pre_2024_rows)

    # Write summary
    with open(SUMMARY_OUTPUT, 'w') as f:
        f.write("BLAST HITS DATE FILTER SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total BLAST hits: {total_hits:,}\n")
        f.write(f"Found in SAbDab: {total_hits - not_in_sabdab:,}\n")
        f.write(f"Not in SAbDab: {not_in_sabdab:,}\n\n")
        f.write(f"Pre-2024 hits: {pre_2024_hits:,} ({pre_2024_hits/total_hits*100:.1f}%)\n")
        f.write(f"2024+ hits: {post_2024_hits:,} ({post_2024_hits/total_hits*100:.1f}%)\n\n")
        f.write(f"Unique PDB IDs (pre-2024): {len(unique_pre_2024_pdbs)}\n")
        f.write(f"Unique PDB IDs (2024+): {len(unique_post_2024_pdbs)}\n\n")
        f.write(f"Output files:\n")
        f.write(f"  - Pre-2024 hits: {PRE_2024_OUTPUT}\n")
        f.write(f"  - Summary: {SUMMARY_OUTPUT}\n")

    print(f"\n{'=' * 50}")
    print(f"RESULTS:")
    print(f"  Total hits: {total_hits:,}")
    print(f"  Pre-2024 hits: {pre_2024_hits:,} ({pre_2024_hits/total_hits*100:.1f}%)")
    print(f"  2024+ hits: {post_2024_hits:,} ({post_2024_hits/total_hits*100:.1f}%)")
    print(f"  Not in SAbDab: {not_in_sabdab:,}")
    print(f"\n  Unique pre-2024 PDBs: {len(unique_pre_2024_pdbs)}")
    print(f"  Unique 2024+ PDBs: {len(unique_post_2024_pdbs)}")
    print(f"\nPre-2024 hits saved to: {PRE_2024_OUTPUT}")
    print(f"Summary saved to: {SUMMARY_OUTPUT}")

if __name__ == '__main__':
    main()
