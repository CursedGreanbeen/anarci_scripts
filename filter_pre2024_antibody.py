#!/usr/bin/env python3
"""Filter for pre-2024 antibody-only matches."""

import csv
from datetime import datetime

# Input files
SABDAB_FILE = "sabdab_summary_all.tsv"
BLAST_HITS_FILE = "results/fasta_fv_only_h3_summary/combined_filtered_hits.csv"

# Output files
OUTPUT_FILE = "results/fasta_fv_only_h3_summary/pre_2024_antibody_hits.csv"
SUMMARY_FILE = "results/fasta_fv_only_h3_summary/pre_2024_antibody_summary.txt"

# Cutoff date
CUTOFF_DATE = datetime(2024, 1, 1)

# Antibody keywords
ANTIBODY_KEYWORDS = [
    'light chain', 'heavy chain', 'fab', 'antibody', 'igg', 'igm', 'iga',
    'bcr', 'immunoglobulin', 'antigen-binding', 'complementarity', 'cdr',
    'variable domain', 'constant domain', 'scfv', 'nanobody', 'vhh',
    'single-domain', 'camelid', 'humanized', 'chimeric', 'monoclonal'
]

def is_antibody(description):
    """Check if hit description indicates an antibody chain."""
    desc = description.lower()
    return any(kw in desc for kw in ANTIBODY_KEYWORDS)

def load_sabdab_dates():
    """Load PDB -> release date mapping from SAbDab TSV."""
    pdb_dates = {}
    with open(SABDAB_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_id = row['pdb'].upper()
            date_str = row['date']
            try:
                release_date = datetime.strptime(date_str, '%m/%d/%y')
                pdb_dates[pdb_id] = release_date
            except ValueError:
                pass
    return pdb_dates

def extract_pdb_id(hit_id):
    """Extract PDB ID from BLAST hit_id."""
    import re
    match = re.search(r'pdb\|([A-Z0-9]{4})\|', hit_id)
    return match.group(1).upper() if match else None

def main():
    print("Loading SAbDab dates...")
    pdb_dates = load_sabdab_dates()
    print(f"Loaded {len(pdb_dates)} PDB entries from SAbDab")

    print("\nFiltering for pre-2024 antibody-only hits...")

    # Counters
    total_hits = 0
    pre_2024_antibody_hits = 0

    # Field names for output
    fieldnames = ['source_file', 'chain_id', 'hit_id', 'hit_description',
                  'pct_identity', 'identity', 'align_length', 'evalue',
                  'bit_score', 'h3_covered', 'h3_identity', 'h3_align_length',
                  'h3_pct_identity', 'release_date', 'pdb_id']

    output_rows = []
    unique_pdbs = set()
    unique_source_files = set()

    with open(BLAST_HITS_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_hits += 1
            pdb_id = extract_pdb_id(row['hit_id'])

            if pdb_id is None:
                continue

            if pdb_id not in pdb_dates:
                continue

            release_date = pdb_dates[pdb_id]

            # Skip if 2024+
            if release_date >= CUTOFF_DATE:
                continue

            # Skip if not antibody
            if not is_antibody(row['hit_description']):
                continue

            pre_2024_antibody_hits += 1
            unique_pdbs.add(pdb_id)
            unique_source_files.add(row['source_file'])

            output_rows.append({
                **row,
                'release_date': release_date.strftime('%Y-%m-%d'),
                'pdb_id': pdb_id
            })

    # Write output
    with open(OUTPUT_FILE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    # Write summary
    with open(SUMMARY_FILE, 'w') as f:
        f.write("PRE-2024 ANTIBODY HITS REPORT\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total BLAST hits scanned: {total_hits:,}\n")
        f.write(f"Pre-2024 antibody hits: {pre_2024_antibody_hits:,}\n\n")
        f.write(f"Unique PDB IDs: {len(unique_pdbs)}\n")
        f.write(f"Unique source files with matches: {len(unique_source_files)}\n\n")
        f.write(f"Output file: {OUTPUT_FILE}\n")

    print(f"\n{'=' * 50}")
    print(f"RESULTS:")
    print(f"  Pre-2024 antibody hits: {pre_2024_antibody_hits:,}")
    print(f"  Unique PDB IDs: {len(unique_pdbs)}")
    print(f"  Source files with matches: {len(unique_source_files)}")
    print(f"\nOutput saved to: {OUTPUT_FILE}")

if __name__ == '__main__':
    main()
