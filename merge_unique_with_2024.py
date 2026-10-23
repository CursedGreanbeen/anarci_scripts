#!/usr/bin/env python3
"""
Merge pre_2024_unique.csv with queries_2024_only_details.csv from each subdirectory
and output to unique_fresh_queries.csv in each subdirectory
"""

import csv
from pathlib import Path

RESULTS_DIR = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results")

# Directories to process (as specified: 0-100, cdr-80-cutoff, levenshtein)
SUBDIRS = [
    "fasta-fv-only-0-100",
    "fasta-fv-only-cdr-80-cutoff",
    "fasta-fv-only-levenshtein"
]

def merge_subdirectory(subdir):
    """Merge pre_2024_unique.csv and queries_2024_only_details.csv in a single subdir."""

    unique_path = RESULTS_DIR / subdir / "pre_2024_unique.csv"
    details_path = RESULTS_DIR / subdir / "queries_2024_only_details.csv"
    output_path = RESULTS_DIR / subdir / "unique_fresh_queries.csv"

    all_rows = []
    fieldnames = None

    # Read pre_2024_unique.csv
    if unique_path.exists():
        print(f"Reading {unique_path}...")
        with open(unique_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            for row in reader:
                all_rows.append(row)
        print(f"  Found {len(all_rows)} pre_2024 unique rows")
    else:
        print(f"Warning: {unique_path} not found")
        return

    # Read queries_2024_only_details.csv
    if details_path.exists():
        print(f"Reading {details_path}...")
        with open(details_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            count_2024 = 0
            for row in reader:
                # Extract pdb_id from source_file (query), not hit_id
                source_file = row.get('source_file', '')
                row['pdb_id'] = source_file.replace('.fasta', '').upper() if source_file else ''
                all_rows.append(row)
                count_2024 += 1
        print(f"  Found {count_2024} 2024 rows")
    else:
        print(f"Warning: {details_path} not found")

    # Write combined output
    print(f"Writing {len(all_rows)} total rows to {output_path}...")

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"  Done!")

def merge_files():
    """Process each subdirectory independently."""

    for subdir in SUBDIRS:
        print(f"\nProcessing {subdir}...")
        merge_subdirectory(subdir)

    print(f"\nAll done!")

if __name__ == "__main__":
    merge_files()
