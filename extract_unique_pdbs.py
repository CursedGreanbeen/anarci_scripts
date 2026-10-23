#!/usr/bin/env python3
"""
Extract unique PDB IDs from unique_fresh_queries.csv in each subdirectory
and write to unique_pdbs.txt files
"""

import csv
from pathlib import Path

RESULTS_DIR = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results")

SUBDIRS = [
    "fasta-fv-only-0-100",
    "fasta-fv-only-cdr-80-cutoff",
    "fasta-fv-only-levenshtein"
]

def extract_pdbs(subdir):
    """Extract unique PDB IDs from unique_fresh_queries.csv."""

    input_path = RESULTS_DIR / subdir / "unique_fresh_queries.csv"
    output_path = RESULTS_DIR / subdir / "unique_pdbs.txt"

    if not input_path.exists():
        print(f"Warning: {input_path} not found")
        return

    print(f"Processing {input_path}...")

    pdb_ids = set()

    with open(input_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row.get('pdb_id', '').strip()
            if pdb_id:
                pdb_ids.add(pdb_id)

    # Sort and write to file
    sorted_pdbs = sorted(pdb_ids)

    with open(output_path, 'w') as f:
        for pdb_id in sorted_pdbs:
            f.write(f"{pdb_id}\n")

    print(f"  Wrote {len(sorted_pdbs)} unique PDB IDs to {output_path}")

def main():
    for subdir in SUBDIRS:
        extract_pdbs(subdir)
    print("\nDone!")

if __name__ == "__main__":
    main()
