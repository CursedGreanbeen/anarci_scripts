#!/usr/bin/env python3
"""Filter for pre-2024 antibody-only matches - CLI version."""

import csv
import sys
import os
from datetime import datetime

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

def load_sabdab_dates(sabdab_file):
    """Load PDB -> release date mapping from SAbDab TSV."""
    pdb_dates = {}
    with open(sabdab_file, 'r') as f:
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
    if len(sys.argv) < 2:
        print("Usage: python filter_pre2024_antibody_cli.py <results_dir> [sabdab_file]")
        print("  results_dir: Directory containing combined_filtered_hits.csv")
        print("  sabdab_file: Path to sabdab_summary_all.tsv (default: ./sabdab_summary_all.tsv)")
        sys.exit(1)

    results_dir = sys.argv[1]
    sabdab_file = sys.argv[2] if len(sys.argv) > 2 else "sabdab_summary_all.tsv"

    # Input files
    blast_hits_file = os.path.join(results_dir, "combined_filtered_hits.csv")

    # Output files
    output_file = os.path.join(results_dir, "pre_2024_antibody_hits.csv")
    summary_file = os.path.join(results_dir, "pre_2024_antibody_summary.txt")

    # Cutoff date
    cutoff_date = datetime(2024, 1, 1)

    print(f"Processing: {results_dir}")
    print(f"Loading SAbDab dates from: {sabdab_file}")
    pdb_dates = load_sabdab_dates(sabdab_file)
    print(f"Loaded {len(pdb_dates)} PDB entries from SAbDab")

    print("\nFiltering for pre-2024 antibody-only hits...")

    # Counters
    total_hits = 0
    pre_2024_antibody_hits = 0

    # Get fieldnames from input file dynamically
    with open(blast_hits_file, 'r') as f:
        reader = csv.DictReader(f)
        input_fieldnames = reader.fieldnames

    # Output fieldnames: input fields + release_date + pdb_id
    output_fieldnames = input_fieldnames + ['release_date', 'pdb_id']

    output_rows = []
    unique_pdbs = set()
    unique_source_files = set()

    with open(blast_hits_file, 'r') as f:
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
            if release_date >= cutoff_date:
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
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    # Write summary
    with open(summary_file, 'w') as f:
        f.write("PRE-2024 ANTIBODY HITS REPORT\n")
        f.write(f"Input directory: {results_dir}\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total BLAST hits scanned: {total_hits:,}\n")
        f.write(f"Pre-2024 antibody hits: {pre_2024_antibody_hits:,}\n\n")
        f.write(f"Unique PDB IDs: {len(unique_pdbs)}\n")
        f.write(f"Unique source files with matches: {len(unique_source_files)}\n\n")
        f.write(f"Output file: {output_file}\n")

    print(f"\n{'=' * 50}")
    print(f"RESULTS:")
    print(f"  Pre-2024 antibody hits: {pre_2024_antibody_hits:,}")
    print(f"  Unique PDB IDs: {len(unique_pdbs)}")
    print(f"  Source files with matches: {len(unique_source_files)}")
    print(f"\nOutput saved to: {output_file}")

if __name__ == '__main__':
    main()
